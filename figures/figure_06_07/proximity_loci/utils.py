from scipy.stats import mannwhitneyu
from scipy.stats import genpareto
from pyliftover import LiftOver
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import polars as pl
import pickle as pk
import numpy as np
import pysam


COMPLEMENT = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]

def liftover_variant(chr:str, bp: int, ref:str, alt:str, fasta, lo):

    REF = ref.upper()
    ALT = alt.upper()

    result = lo.convert_coordinate(chr, bp - 1)

    if not result:
        print(chr, bp, 'failed')
        global fail
        fail = fail + 1
        return 000, REF, ALT

    new_chrom, new_pos0, strand, _ = result[0]
    new_pos = new_pos0 + 1  # back to 1-based

    # If the liftover flips strand, REF/ALT need to be reverse-complemented
    # before comparing against the target's forward-strand sequence.
    exp_ref, exp_alt = (REF, ALT) if strand == "+" else (revcomp(REF), revcomp(ALT))

    fetched_ref = fasta.fetch(new_chrom, new_pos0, new_pos0 + len(exp_ref)).upper()
    
    if fetched_ref == exp_ref:
        return new_pos, exp_ref, exp_alt
    else:
        return new_pos, exp_alt, exp_ref  # swap so new_ref matches the target


def match_distribution(reference_df, candidate_df, chrom_col="chr", maf_col="maf",
                        n_bins=10, replace=False, random_state=None, oversample_factor=100):
    """
    Sample rows from candidate_df so the joint (chromosome x MAF-bin)
    distribution of the sample matches reference_df.

    MAF bin edges are derived from reference_df (quantile bins, so each bin
    holds roughly equal numbers of reference variants), then applied to
    both reference_df and candidate_df so every variant falls into a
    comparable bin.
    """
    rng = np.random.default_rng(random_state)

    # Quantile bin edges from the reference set's MAF distribution
    _, bin_edges = pd.qcut(reference_df[maf_col], q=n_bins, retbins=True, duplicates="drop")

    ref = reference_df.copy()
    cand = candidate_df.copy()
    ref["maf_bin"] = pd.cut(ref[maf_col], bins=bin_edges, include_lowest=True)
    cand["maf_bin"] = pd.cut(cand[maf_col], bins=bin_edges, include_lowest=True)

    # How many reference variants fall in each (chrom, maf_bin) stratum
    target_counts = ref.groupby([chrom_col, "maf_bin"], observed=True).size()

    sampled_parts = []
    shortfalls = []

    for (chrom, mbin), n_needed in target_counts.items():
        pool = cand[(cand[chrom_col] == chrom) & (cand["maf_bin"] == mbin)]

        if len(pool) == 0:
            shortfalls.append((chrom, mbin, n_needed * oversample_factor, 0))
            continue

        n_draw = n_needed * oversample_factor
        if not replace and len(pool) < n_needed * oversample_factor:
            shortfalls.append((chrom, mbin, n_needed * oversample_factor, len(pool)))
            n_draw = len(pool)  # take everything available in this stratum

        idx = rng.choice(pool.index, size=n_draw, replace=replace)
        sampled_parts.append(cand.loc[idx])

    sampled = pd.concat(sampled_parts, ignore_index=True) if sampled_parts else cand.iloc[0:0].copy()

    if shortfalls:
        print(f"[match_distribution] {len(shortfalls)} stratum/strata under-supplied:")
        for chrom, mbin, needed, available in shortfalls[:20]:
            print(f"    chr={chrom} maf_bin={mbin}: needed {needed}, only {available} available")

    return sampled.drop(columns="maf_bin")


 
def build_locus_metrics(snps_df, peaks_df=None, locus_col="locus_id", variant_col="variant_id",
                         maf_col="MAF", chr_col="chr", start_col="start", stop_col="stop",
                         peak_chr_col="chr", peak_start_col="start", peak_end_col="end"):
    """
    Collapse a variant-level dataframe into per-locus metrics: snp_count,
    mean_maf, size (stop - start), and SNP_density (snp_count / size).
 
    If `peaks_df` (chr, start, end, ...) is given, also adds peak_count and
    peak_density (peak_count / size): the number of peaks overlapping each
    locus (interval overlap: peak.start <= locus.stop AND peak.end >=
    locus.start), normalized by locus size.
 
    chr/start/stop are constant within a locus, so they're pulled via
    drop_duplicates rather than aggregated with the count/mean.
    """
    metrics = (
        snps_df.groupby(locus_col)
        .agg(snp_count=(variant_col, "count"), mean_maf=(maf_col, "mean"))
        .reset_index()
    )
 
    bounds = snps_df[[locus_col, chr_col, start_col, stop_col]].drop_duplicates(subset=locus_col).copy()
    bounds = bounds.rename(columns={chr_col: "chr", start_col: "start", stop_col: "stop"})
    bounds["size"] = bounds["stop"] - bounds["start"]
 
    metrics = metrics.merge(bounds[[locus_col, "size"]], on=locus_col, how="left")
    metrics["SNP_density"] = metrics["snp_count"] / metrics["size"]
 
    if peaks_df is not None:
        peak_counts = _count_peaks_per_locus(bounds, peaks_df, locus_col,
                                              peak_chr_col, peak_start_col, peak_end_col)
        metrics = metrics.merge(peak_counts, on=locus_col, how="left")
        metrics["peak_count"] = metrics["peak_count"].fillna(0).astype(int)
        metrics["peak_density"] = metrics["peak_count"] / metrics["size"]
 
    return metrics
 
 
def _build_peak_chrom_index(peaks_df, peak_chr_col="chr", peak_start_col="start", peak_end_col="end"):
    """
    Sort peaks by start position once per chromosome. Assumes peaks within
    a chromosome are already non-overlapping (true for merged/'consensus
    peak' sets), so end positions are monotonic in start order too --
    that's what makes the binary search below valid.
    """
    idx = {}
    for chrom in peaks_df[peak_chr_col].unique():
        sub = peaks_df[peaks_df[peak_chr_col] == chrom].sort_values(peak_start_col)
        idx[chrom] = (sub[peak_start_col].to_numpy(), sub[peak_end_col].to_numpy())
    return idx
 
 
def _count_peaks_per_locus(bounds, peaks_df, locus_col, peak_chr_col="chr",
                            peak_start_col="start", peak_end_col="end"):
    """
    Number of peaks overlapping each locus in `bounds` (chr/start/stop),
    via binary search against a pre-sorted per-chromosome peak index --
    NOT a join. A join here builds a full cross-product of every locus
    against every peak sharing a chromosome before filtering, which is
    exactly the pattern that caused this to blow up in time/memory with a
    large peaks_df -- this replaces it with the same binary-search
    approach used for variant-in-locus lookups earlier in this pipeline.
 
    Assumes peaks within each chromosome don't overlap each other; if your
    peaks_df has raw overlapping peaks, merge them first (e.g. `bedtools
    merge`) before calling this.
    """
    chrom_index = _build_peak_chrom_index(peaks_df, peak_chr_col, peak_start_col, peak_end_col)
    counts = []
    for row in bounds.itertuples(index=False):
        peak_starts, peak_ends = chrom_index.get(row.chr, (np.array([]), np.array([])))
        lo = np.searchsorted(peak_ends, row.start, side="left")
        hi = np.searchsorted(peak_starts, row.stop, side="right")
        counts.append(max(0, hi - lo))
    return pd.DataFrame({locus_col: bounds[locus_col].to_numpy(), "peak_count": counts})
 
 
def _linear_bin_edges(values, n_bins):
    """Equal-width (linear) bin edges spanning the observed range of `values`."""
    values = np.asarray(values, dtype=float)
    vmin, vmax = values.min(), values.max()
    if vmax <= vmin:
        vmax = vmin + 1e-9
    return np.linspace(vmin, vmax, n_bins + 1)
 
 
def _exp_bin_edges(values, n_bins):
    """Exponentially (log-)spaced bin edges spanning the observed range of `values`."""
    values = np.asarray(values, dtype=float)
    vmin, vmax = values.min(), values.max()
    floor = vmin
    if floor <= 0:
        positive = values[values > 0]
        floor = positive.min() * 0.5 if positive.size else 1e-9
    if vmax <= floor:
        vmax = floor * 1.0001
    edges = np.geomspace(floor, vmax, n_bins + 1)
    # The floor above is only there to keep geomspace's math well-defined
    # (it can't start at 0). It must never exceed the true minimum, or
    # values at/below it (e.g. loci with 0 overlapping peaks) fall outside
    # every bin and get silently dropped as NaN by pd.cut.
    edges[0] = min(edges[0], vmin)
    return edges
 
 
def sample_null_regions(reference_metrics, candidate_metrics, locus_col="locus_id",
                         metric_cols=("mean_maf", "size", "SNP_density", "peak_density"),
                         linear_cols=("mean_maf",), n_bins=5, M=100,
                         random_state=None, verbose=True):
    """
    Draw M independent null sets of len(reference_metrics) loci from
    candidate_metrics, each stratified-sampled so its joint distribution
    over `metric_cols` matches reference_metrics. See build_locus_metrics
    for how metric_cols are computed.
 
    snp_count is deliberately left out of the default metric_cols since
    it's redundant with size + SNP_density together (SNP_density =
    snp_count / size) -- stratifying on all three at once would just
    multiply the number of strata for no extra matching power. Pass your
    own metric_cols if you want a different combination.
 
    Returns
    -------
    (null_draws, shortfalls, bin_edges)
    null_draws : list of length M, each a list of locus_ids.
    shortfalls : list of dicts describing any under-supplied strata.
    bin_edges : dict of {metric_col: edges_array} -- the actual bin edges
        used to stratify each metric, e.g. for plotting the same bins
        used for matching instead of independently recomputed ones.
    """
    rng = np.random.default_rng(random_state)
 
    ref = reference_metrics.copy()
    cand = candidate_metrics.copy()
 
    bin_cols = []
    bin_edges = {}
    for col in metric_cols:
        col_values = ref[col].to_numpy()
        finite_mask = np.isfinite(col_values)
        n_nonfinite = (~finite_mask).sum()
        if n_nonfinite:
            print(f"[sample_null_regions] {n_nonfinite} non-finite value(s) (NaN/inf) in "
                  f"reference '{col}' -- excluded when computing bin edges. This is often "
                  f"caused by a zero-size locus (size=0 -> density=count/0=inf); worth "
                  f"checking your locus bounds for that.")
        finite_values = col_values[finite_mask]
        true_min, true_max = finite_values.min(), finite_values.max()
 
        if col in linear_cols:
            edges = _linear_bin_edges(finite_values, n_bins)
        else:
            edges = _exp_bin_edges(finite_values, n_bins)
        # Hard guarantee regardless of which branch/degenerate case produced
        # `edges`: never extend beyond the true observed reference range.
        edges[0] = min(edges[0], true_min)
        edges[-1] = true_max
 
        bin_edges[col] = edges
        bin_col = f"{col}_bin"
        ref[bin_col] = pd.cut(ref[col], bins=edges, include_lowest=True)
        cand[bin_col] = pd.cut(cand[col], bins=edges, include_lowest=True)
        bin_cols.append(bin_col)
 
    strata = []
    for stratum_vals, ref_group in ref.groupby(bin_cols, observed=True):
        stratum_vals = stratum_vals if isinstance(stratum_vals, tuple) else (stratum_vals,)
        mask = pd.Series(True, index=cand.index)
        for col, val in zip(bin_cols, stratum_vals):
            mask &= (cand[col] == val)
        pool = cand[mask]
        strata.append({
            "stratum": stratum_vals,
            "reference_locus_ids": ref_group[locus_col].tolist(),
            "n_needed": len(ref_group),
            "n_available": len(pool),
            "pool_index": pool.index,
            "candidate_locus_ids_available": pool[locus_col].tolist(),
        })
 
    shortfalls = [s for s in strata if s["n_available"] < s["n_needed"]]
 
    if verbose and shortfalls:
        print(f"[sample_null_regions] {len(shortfalls)}/{len(strata)} strata are short on candidates:")
        for s in shortfalls:
            print(f"    stratum={s['stratum']}  needed={s['n_needed']}  available={s['n_available']}")
            print(f"        reference loci affected: {s['reference_locus_ids']}")
 
    null_draws = []
    for draw_i in range(M):
        sampled_parts = []
        for s in strata:
            n_draw = min(s["n_needed"], s["n_available"])
            if n_draw == 0:
                continue
            idx = rng.choice(s["pool_index"], size=n_draw, replace=False)
            sampled_parts.append(cand.loc[idx, locus_col])
        null_draws.append(pd.concat(sampled_parts).tolist() if sampled_parts else [])
 
    return null_draws, shortfalls, bin_edges


def plot_strat_density(lead_metrics, null_candidate_metrics, null_draws, metric, p_val, bin_edges=None):
    plt.figure(figsize=(12, 3))

    lo, hi = lead_metrics[metric].min(), lead_metrics[metric].max()
    bins = np.linspace(lo, hi, 30)

    for draw in null_draws:
        drawx = null_candidate_metrics[null_candidate_metrics['locus_id'].isin(draw)]
        sns.histplot(drawx[metric], label='null-set 0', bins=bins, element='step', fill=False, color='#aaa', alpha=0.05)

    sns.histplot(lead_metrics[metric], label='leads', bins=bins, element='step', fill=False, color='k')
    for edge in bin_edges[metric]:
        plt.axvline(edge, color='r', linestyle=':', linewidth=2, alpha=0.5)
    plt.title(f'{metric} of lead loci set vs all drawn sets [p different={p_val}]')
    plt.savefig(f"5_figs/strat_density_{metric}.png")


def pool_variants(variants: pl.DataFrame, loci: pl.DataFrame, value_col) -> pl.Series:
    """All variant values falling inside any locus in `loci`."""
    hits = variants.join(loci, on="chr", how="inner").filter(
        (pl.col("pos") >= pl.col("left_bound")) & (pl.col("pos") <= pl.col("right_bound"))
    )
    return hits[value_col]


def mwu_lead_vs_pooled_null(lead_metrics, null_candidate_metrics, null_draws, col="mean_maf"):
    pooled_null = pd.concat([
        null_candidate_metrics.loc[null_candidate_metrics["locus_id"].isin(draw), col]
        for draw in null_draws
    ])
    return mannwhitneyu(lead_metrics[col], pooled_null, alternative="two-sided")


def mwu_lead_vs_pooled_null_variants(variants, lead_loci, null_loci, null_draws, value_col):
    lead_values = pool_variants(variants, lead_loci, value_col)
    pooled_null = pl.concat([
        pool_variants(variants, null_loci.filter(pl.col("lead_id").is_in(lead_ids)), value_col)
        for lead_ids in null_draws.values()
    ])
    return mannwhitneyu(lead_values, pooled_null, alternative="two-sided")


def fit_gpd_shape(values, threshold=0.5, min_exceedances=10):
    """
    Fit a Generalized Pareto Distribution to the excesses of abs(values)
    over `threshold`, with loc fixed at 0 (excesses are threshold-centered
    by construction). Returns the shape parameter -- the tail "fatness":
    >0 heavy tail, 0 exponential decay, <0 bounded tail. Returns nan if
    there aren't enough exceedances to fit stably.
    """
    abs_values = np.abs(np.asarray(values))
    excesses = abs_values[abs_values > threshold] - threshold
    if len(excesses) < min_exceedances:
        return np.nan
    shape, loc, scale = genpareto.fit(excesses, floc=0)
    return shape
 
 
def gpd_tail_enrichment(variants, lead_loci, null_loci, null_draws, value_col,
                         threshold=0.5, min_exceedances=10):
    """
    Fit a GPD shape parameter to the tail (|value| > threshold) of the
    pooled lead-loci variants (observed), and to each bootstrap draw's
    pooled null-loci variants (M null shape parameters). Empirical
    two-sided p-value compares the observed shape against the spread of
    the M null shapes.
    """
    lead_values = pool_variants(variants, lead_loci, value_col)
    observed_shape = fit_gpd_shape(lead_values, threshold, min_exceedances)
 
    null_shapes = np.array([
        fit_gpd_shape(
            pool_variants(variants, null_loci.filter(pl.col("lead_id").is_in(lead_ids)), value_col),
            threshold, min_exceedances,
        )
        for lead_ids in null_draws.values()
    ])
 
    n_dropped = np.isnan(null_shapes).sum()
    if n_dropped:
        print(f" [enrich] {n_dropped}/{len(null_shapes)} draws had fewer than "
              f"{min_exceedances} exceedances above threshold={threshold} -- excluded from null")
    null_shapes = null_shapes[~np.isnan(null_shapes)]
 
    null_mean = null_shapes.mean()
    p_value = (1 + np.sum(np.abs(null_shapes - null_mean) >= np.abs(observed_shape - null_mean))) / (len(null_shapes) + 1)
 
    return {
        "observed_shape": observed_shape,
        "null_mean": null_mean,
        "null_std": null_shapes.std(),
        "null_shapes": null_shapes,
        "p_value": p_value,
    }
 

def parse_peaks(peak_strings) -> pl.DataFrame:
    """Parse 'chr:start-end' strings (e.g. 'chr2:202419829-202420329') into a chr/start/end DataFrame."""
    df = pl.DataFrame({"peak": list(peak_strings)})
    region = pl.col("peak").str.split(":").list.get(1)
    return df.with_columns([
        pl.col("peak").str.split(":").list.get(0).alias("chr"),
        region.str.split("-").list.get(0).cast(pl.Int64).alias("start"),
        region.str.split("-").list.get(1).cast(pl.Int64).alias("end"),
    ])
 

def per_locus_peak_counts(peaks: pl.DataFrame, loci: pl.DataFrame) -> np.ndarray:
    """Number of peaks overlapping each individual locus in `loci` (one value per locus, including 0)."""
    hits = peaks.join(loci, on="chr", how="inner").filter(
        (pl.col("start") <= pl.col("right_bound")) & (pl.col("end") >= pl.col("left_bound"))
    )
    counts = hits.group_by("lead_id").len().rename({"len": "peak_count"})
    per_locus = (
        loci.select("lead_id")
        .join(counts, on="lead_id", how="left")
        .with_columns(pl.col("peak_count").fill_null(0))
    )
    return per_locus["peak_count"].to_numpy()
 
 
def fit_nb_dispersion(counts, min_n=10):
    """
    Method-of-moments estimate of the Negative Binomial (NB2) dispersion
    parameter alpha, where variance = mean + alpha * mean^2. alpha=0 is
    Poisson; alpha>0 is overdispersion (a fatter tail than Poisson would
    predict) -- this is the "fatness" parameter for count data, playing
    the same role the GPD shape parameter played for the continuous logfc
    tail earlier.
 
    Method-of-moments (not MLE) is used deliberately: it's closed-form, so
    it can't fail to converge across many bootstrap draws the way an MLE
    optimizer sometimes can on small/sparse samples. If you want the MLE
    version instead, statsmodels' NegativeBinomial (intercept-only) gives
    mu/alpha directly, at the cost of needing to handle non-convergence.
 
    Returns nan if there are too few loci to estimate variance reliably.
    """
    counts = np.asarray(counts, dtype=float)
    if len(counts) < min_n:
        return np.nan
    mean = counts.mean()
    if mean == 0:
        return np.nan
    var = counts.var(ddof=1)
    alpha = (var - mean) / mean ** 2
    return max(alpha, 0.0)  # NB alpha can't be negative; clip sample noise that dips below 0
 
 
def nb_dispersion_enrichment(peaks, lead_loci, null_loci, null_draws, min_n=10):
    """
    Fit the NB dispersion parameter once to the real lead loci's per-locus
    peak counts (observed), and once to each bootstrap draw's null loci
    (M null values). Empirical two-sided p-value compares the observed
    dispersion against the spread of the M null dispersions.
    """
    lead_counts = per_locus_peak_counts(peaks, lead_loci)
    observed_alpha = fit_nb_dispersion(lead_counts, min_n)
 
    null_alphas = np.array([
        fit_nb_dispersion(per_locus_peak_counts(peaks, null_loci.filter(pl.col("lead_id").is_in(lead_ids))), min_n)
        for lead_ids in null_draws.values()
    ])
 
    n_dropped = np.isnan(null_alphas).sum()
    if n_dropped:
        print(f"   [nb_dispersion_enrichment] {n_dropped}/{len(null_alphas)} draws had fewer than "
              f"{min_n} loci -- excluded from null")
    null_alphas = null_alphas[~np.isnan(null_alphas)]
 
    null_mean = null_alphas.mean()
    p_value = (1 + np.sum(np.abs(null_alphas - null_mean) >= np.abs(observed_alpha - null_mean))) / (len(null_alphas) + 1)
 
    return {
        "observed_alpha": observed_alpha,
        "null_mean": null_mean,
        "null_std": null_alphas.std(),
        "null_alphas": null_alphas,
        "p_value": p_value,
    }
 
 
