#!/usr/bin/env python3

"""Variant scoring script."""

import gc
import json
import h5py
import argparse
import numpy as np
import pandas as pd
from pathlib import Path
from keras.models import load_model

from variant_utils import *


def maybe_str_or_int(arg):
    """
    Convert argument to int or return 'all'.

    Args:
        arg (str): Argument to convert.

    Returns:
        int or str: Converted value.

    Raises:
        argparse.ArgumentTypeError: If the argument cannot be converted to int or is not 'all'.
    """
    try:
        return int(arg)
    except ValueError:
        pass
    if arg == "all":
        return arg
    raise argparse.ArgumentTypeError("Number must be an int or 'all'")


def fetch_variant_args():
    """
    Fetch and parse command line arguments for interpretations.

    Returns:
        argparse.Namespace: Parsed command line arguments.
    """
    parser = argparse.ArgumentParser(description="Perform variant scoring with a trained model")

    parser.add_argument("-v", "--variants", type=str, required=True, help="a TSV file containing a list of variants to score")
    parser.add_argument("-sc", "--schema", type=str, choices=['snp', 'peak', 'indel'], default='peak', help="Format for the input variants list")
    parser.add_argument("-m", "--model_path", type=str, required=True, help="Path to trained model for variant scoring")
    parser.add_argument("-g", "--genome", type=str, required=True, help="Genome fasta")
    parser.add_argument("-chrs", "--chrom_sizes", type=str, required=True, help="Path to TSV file with chromosome sizes")
    parser.add_argument("-o", "--output_path", type=str, required=True, help="Output path")

    parser.add_argument("-p", "--peaks", type=str, required=False, nargs='?', const='', default=None, help="3-column BED file containing peak regions")
    parser.add_argument("-pg", "--peak_genome", type=str, required=False, help="Genome fasta for peaks")
    parser.add_argument("-ps", "--peak_chrom_sizes", type=str, required=False, help="Path to TSV file with chromosome sizes for peak genome")

    parser.add_argument("-me", "--method", type=str, choices=['shap', 'ism'], default='ism', help="Method to get variant scores")
    parser.add_argument("-al", "--attribution_length", type=int, default=0, help="Window around SNP to compute attribution scores (for shap)")
    parser.add_argument("-fo", "--forward_only", action='store_true', help="Run variant scoring only on forward sequence (for ism)")
    parser.add_argument("-norm", "--normalize", action="store_true", help="Gradient normalization (for shap)")

    parser.add_argument("-ns", "--num_shuffles", type=int, default=0, help="Number of permuted scores per SNP")
    parser.add_argument("-np", "--num_peaks", type=maybe_str_or_int, default="all", help="Number of peaks to use for peak percentile calculation")

    parser.add_argument("-b", "--batch_size", type=int, default=1024, help="Batch size")
    parser.add_argument("-s", "--seed", type=int, default=42, help="Seed for reproducibility when sampling")

    parser.add_argument("-chr", "--chrom", type=str, help="Only score SNPs in selected chromosome")
    parser.add_argument("--hdf5", action='store_true', help="Save detailed predictions in hdf5 file")

    args = parser.parse_args()
    return args


def main(args):
    if args.method.lower() == 'shap':
        tf.compat.v1.disable_v2_behavior()
        tf.compat.v1.disable_eager_execution()

    output_path = Path(args.output_path)
    Path(output_path).mkdir(parents=True, exist_ok=True)

    # Write all the command line arguments to a json file
    with open(output_path.joinpath(Path("variant_params.json")), "w") as f:
        json.dump(vars(args), f, ensure_ascii=False, indent=4)

    np.random.seed(args.seed)
    if args.forward_only:
        print("Running variant scoring only for forward sequences")

    # Load model
    model = load_model(args.model_path, compile=False)
    print("Got the model")
    model.summary()

    # Infer input length
    input_length = model.input_shape[1]
    print("Inferred model input length:", input_length)

    # Load variants
    variants_table = load_variant_table(args.variants, args.schema)
    variants_table = variants_table.fillna('-')
    print("Original variants table shape:", variants_table.shape)

    chrom_sizes = pd.read_csv(args.chrom_sizes, header=None, sep='\t', names=['chrom', 'size'])
    chrom_sizes_dict = chrom_sizes.set_index('chrom')['size'].to_dict()

    # Check only specific chromosome
    if args.chrom:
        variants_table = variants_table.loc[variants_table['chr'] == args.chrom]
        print("Chromosome variants table shape:", variants_table.shape)

    # Get valid variants
    if args.schema == 'snp':
        print("Running variant scoring centered on SNPs")
        variants_table = variants_table.loc[variants_table.apply(lambda x: get_variants(x.chr, x.pos, input_length, chrom_sizes_dict), axis=1)]
    elif args.schema == 'peak':
        print("Running variant scoring on the same peaks")
        variants_table = variants_table.loc[variants_table.apply(lambda x: get_peaks(x.chr, x.start, x.end, input_length, chrom_sizes_dict), axis=1)]
    elif args.schema == 'indel':
        print("Running indel scoring on the same peaks")
        variants_table = variants_table.loc[variants_table.apply(lambda x: get_peaks(x.chr, x.start, x.end, input_length, chrom_sizes_dict), axis=1)]
    else:
        raise NameError("Schema should be either 'snp' or 'peak' or 'indel'.")
    variants_table.reset_index(drop=True, inplace=True)
    print("Final variants table shape:", variants_table.shape)

    # Permute other variants to test the significance of the variants across contexts
    shuf_variants_table = create_shuffle_table(variants_table, args.seed, args.num_shuffles)
    print("Shuffled variants table shape:", shuf_variants_table.shape)

    if len(shuf_variants_table) > 0:
        if args.method.lower() == 'shap':
            shuf_variant_ids, shuf_ref_pred, shuf_alt_pred = fetch_variant_shap(
                model,
                shuf_variants_table,
                input_length,
                args.genome,
                args.batch_size,
                shuffle=True,
                schema=args.schema,
                attribution_length=args.attribution_length,
                normalize=args.normalize
            )
        elif args.method.lower() == 'ism':
            shuf_variant_ids, shuf_ref_pred, shuf_alt_pred = fetch_variant_predictions(
                model,
                shuf_variants_table,
                input_length,
                args.genome,
                args.batch_size,
                shuffle=True,
                forward_only=args.forward_only,
                schema=args.schema
            )
        else:
            raise NameError("Method should be either 'shap' or 'ism'.")

        print('Got shuffled variant predictions!')

    if args.peaks:
        if args.peak_chrom_sizes is None:
            args.peak_chrom_sizes = args.chrom_sizes
        if args.peak_genome is None:
            args.peak_genome = args.genome

        peak_chrom_sizes = pd.read_csv(args.peak_chrom_sizes, header=None, sep='\t', names=['chrom', 'size'])
        peak_chrom_sizes_dict = peak_chrom_sizes.set_index('chrom')['size'].to_dict()

        peaks = pd.read_csv(args.peaks, header=None, sep='\t', usecols=[0, 1, 2], names=['chr', 'start', 'end'])
        print("Original peak table shape:", peaks.shape)

        peaks.sort_values(by=['chr', 'start', 'end'], ascending=[True, True, True], inplace=True)
        peaks.drop_duplicates(subset=['chr', 'start', 'end'], inplace=True)
        peaks = peaks.loc[peaks.apply(lambda x: get_peaks(x.chr, x.start, x.end, input_length, peak_chrom_sizes_dict), axis=1)]
        peaks.reset_index(drop=True, inplace=True)
        print("De-duplicated peak table shape:", peaks.shape)

        if isinstance(args.num_peaks, int):
            num = min(args.num_peaks, len(peaks))
            peaks = peaks.sample(n=num, random_state=args.seed, ignore_index=True)
            print("Subsampled peak table shape:", peaks.shape)

        if args.method.lower() == 'shap':
            preds = fetch_peak_shap(
                model,
                peaks,
                input_length,
                args.peak_genome,
                args.batch_size,
                normalize=args.normalize
            )
        elif args.method.lower() == 'ism':
            preds = fetch_peak_predictions(
                model,
                peaks,
                input_length,
                args.peak_genome,
                args.batch_size,
                forward_only=args.forward_only
            )
        else:
            raise NameError("Method should be either 'shap' or 'ism'.")

        print('Got peak predictions!')

        if len(shuf_variants_table) > 0:
            shuf_logfc, shuf_diff, shuf_ref_percentile, shuf_alt_percentile = get_variant_scores_with_peaks(
                shuf_ref_pred,
                shuf_alt_pred,
                preds
            )

            shuf_max_percentile = np.maximum(shuf_ref_percentile, shuf_alt_percentile)
            shuf_percentile_change = shuf_alt_percentile - shuf_ref_percentile
            shuf_abs_percentile_change = np.abs(shuf_percentile_change)
            shuf_abs_logfc = np.abs(shuf_logfc)
            shuf_abs_diff = np.abs(shuf_diff)
            shuf_logfc_max_percentile = shuf_logfc * shuf_max_percentile
            shuf_abs_logfc_max_percentile = shuf_abs_logfc * shuf_max_percentile
            shuf_diff_max_percentile = shuf_diff * shuf_max_percentile
            shuf_abs_diff_max_percentile = shuf_abs_diff * shuf_max_percentile

            assert shuf_abs_logfc.shape == shuf_logfc.shape
            assert shuf_abs_logfc_max_percentile.shape == shuf_abs_logfc.shape
            assert shuf_abs_logfc_max_percentile.shape == shuf_max_percentile.shape
            assert shuf_max_percentile.shape == shuf_ref_percentile.shape
            assert shuf_max_percentile.shape == shuf_alt_percentile.shape
            assert shuf_max_percentile.shape == shuf_percentile_change.shape
            assert shuf_abs_percentile_change.shape == shuf_percentile_change.shape
            assert shuf_abs_diff.shape == shuf_diff.shape
            assert shuf_abs_diff_max_percentile.shape == shuf_abs_diff.shape
            assert shuf_abs_diff_max_percentile.shape == shuf_max_percentile.shape

    else:
        if len(shuf_variants_table) > 0:
            shuf_logfc, shuf_diff = get_variant_scores(
                shuf_ref_pred,
                shuf_alt_pred)

            shuf_abs_logfc = np.squeeze(np.abs(shuf_logfc))
            shuf_abs_diff = np.squeeze(np.abs(shuf_diff))

    # Get model prediction for variants
    ref_peak_pred, alt_peak_pred = None, None
    if args.method.lower() == 'shap':
        variant_ids, (ref_pred, alt_pred), (ref_peak_pred, alt_peak_pred), \
            ((ref_shap, ref_seqs), (alt_shap, alt_seqs)) = fetch_variant_shap(
            model,
            variants_table,
            input_length,
            args.genome,
            args.batch_size,
            shuffle=False,
            schema=args.schema,
            attribution_length=args.attribution_length,
            normalize=args.normalize
        )

        print("Saving 'contribution' scores")
        np.savez(output_path.joinpath(Path('ref_shap.npz')), ref_shap)
        np.savez(output_path.joinpath(Path('ref_ohe.npz')), ref_seqs)

        np.savez(output_path.joinpath(Path('alt_shap.npz')), alt_shap)
        np.savez(output_path.joinpath(Path('alt_ohe.npz')), alt_seqs)

        del ref_shap, alt_shap, ref_seqs, alt_seqs
        gc.collect()

    elif args.method.lower() == 'ism':
        variant_ids, ref_pred, alt_pred = fetch_variant_predictions(
            model,
            variants_table,
            input_length,
            args.genome,
            args.batch_size,
            shuffle=False,
            forward_only=args.forward_only,
            schema=args.schema
        )
    else:
        raise NameError("Method should be either 'shap' or 'ism'.")

    print('Got variant predictions!')

    if args.peaks:
        (logfc, diff, ref_percentile, alt_percentile), \
            (logfc_peak, diff_peak, ref_peak_percentile, alt_peak_percentile) = get_variant_scores_with_peaks(
            ref_pred,
            alt_pred,
            preds,
            ref_peak_pred,
            alt_peak_pred
        )

    else:
        (logfc, diff), (logfc_peak, diff_peak) = get_variant_scores(
            ref_pred,
            alt_pred,
            ref_peak_pred,
            alt_peak_pred
        )

    assert np.array_equal(variants_table["variant_id"].tolist(), variant_ids)
    max_variants_table = variants_table.copy()
    max_variants_table['max_delta_ref'] = diff.min(axis=1)
    max_variants_table['max_delta_alt'] = diff.max(axis=1)
    max_variants_table['max_delta'] = np.abs(diff).max(axis=1)
    max_variants_table['max_df_topic_ref'] = diff.argmin(axis=1) + 1
    max_variants_table['max_df_topic_alt'] = diff.argmax(axis=1) + 1
    max_variants_table['max_df_topic'] = np.abs(diff).argmax(axis=1) + 1

    if args.method.lower() == 'shap':
        max_variants_table['max_delta_peak_ref'] = diff_peak.min(axis=1)
        max_variants_table['max_delta_peak_alt'] = diff_peak.max(axis=1)
        max_variants_table['max_delta_peak'] = np.abs(diff_peak).max(axis=1)
        max_variants_table['max_df_peak_topic_ref'] = diff_peak.argmin(axis=1) + 1
        max_variants_table['max_df_peak_topic_alt'] = diff_peak.argmax(axis=1) + 1
        max_variants_table['max_df_peak_topic'] = np.abs(diff_peak).argmax(axis=1) + 1

    max_variants_table['max_fc_ref'] = logfc.min(axis=1)
    max_variants_table['max_fc_alt'] = logfc.max(axis=1)
    max_variants_table['max_fc'] = np.abs(logfc).max(axis=1)
    max_variants_table['max_fc_topic_ref'] = logfc.argmin(axis=1) + 1
    max_variants_table['max_fc_topic_alt'] = logfc.argmax(axis=1) + 1
    max_variants_table['max_fc_topic'] = np.abs(logfc).argmax(axis=1) + 1

    if args.method.lower() == 'shap':
        max_variants_table['max_fc_peak_ref'] = logfc_peak.min(axis=1)
        max_variants_table['max_fc_peak_alt'] = logfc_peak.max(axis=1)
        max_variants_table['max_fc_peak'] = np.abs(logfc_peak).max(axis=1)
        max_variants_table['max_fc_peak_topic_ref'] = logfc_peak.argmin(axis=1) + 1
        max_variants_table['max_fc_peak_topic_alt'] = logfc_peak.argmax(axis=1) + 1
        max_variants_table['max_fc_peak_topic'] = np.abs(logfc_peak).argmax(axis=1) + 1

    variants_table["ref_pred"] = ref_pred.tolist()
    variants_table["alt_pred"] = alt_pred.tolist()
    variants_table["diff"] = diff.tolist()
    variants_table["logfc"] = logfc.tolist()

    if args.method.lower() == 'shap':
        variants_table["ref_peak_pred"] = ref_peak_pred.tolist()
        variants_table["alt_peak_pred"] = alt_peak_pred.tolist()
        variants_table["diff_peak"] = diff_peak.tolist()
        variants_table["logfc_peak"] = logfc_peak.tolist()

    if len(shuf_variants_table) > 0:
        variants_table["logfc.pval"] = get_pvals(variants_table["logfc"].tolist(), shuf_logfc, tail="both")
        variants_table["diff.pval"] = get_pvals(variants_table["diff"].tolist(), shuf_diff, tail="both")
        variants_table["abs_logfc.pval"] = get_pvals(variants_table["abs_logfc"].tolist(), shuf_abs_logfc, tail="right")
        variants_table["abs_diff.pval"] = get_pvals(variants_table["abs_diff"].tolist(), shuf_abs_diff, tail="right")
    if args.peaks:
        variants_table["ref_percentile"] = ref_percentile.tolist()
        variants_table["alt_percentile"] = alt_percentile.tolist()
        variants_table["max_percentile"] = np.maximum(ref_percentile, alt_percentile).tolist()
        variants_table["percentile_change"] = (alt_percentile - ref_percentile).tolist()
        variants_table["diff_x_max_percentile"] = (diff * np.maximum(ref_percentile, alt_percentile)).tolist()
        variants_table["diff_x_percentile_change"] = (diff * np.abs(alt_percentile - ref_percentile)).tolist()
        variants_table["logfc_x_max_percentile"] = (logfc * np.maximum(ref_percentile, alt_percentile)).tolist()
        variants_table["logfc_x_percentile_change"] = (logfc * np.abs(alt_percentile - ref_percentile)).tolist()

        max_variants_table['max_percdf_topic_ref'] = (alt_percentile - ref_percentile).argmin(axis=1) + 1
        max_variants_table['max_percdf_topic_alt'] = (alt_percentile - ref_percentile).argmax(axis=1) + 1
        max_variants_table['max_percdf_topic'] = np.abs(alt_percentile - ref_percentile).argmax(axis=1) + 1
        max_variants_table['max_df_x_perc_topic_ref'] = (diff * np.maximum(ref_percentile, alt_percentile)).argmin(axis=1) + 1
        max_variants_table['max_df_x_perc_topic_alt'] = (diff * np.maximum(ref_percentile, alt_percentile)).argmax(axis=1) + 1
        max_variants_table['max_df_x_perc_topic'] = (np.abs(diff) * np.maximum(ref_percentile, alt_percentile)).argmax(axis=1) + 1
        max_variants_table['max_df_x_percdf_topic_ref'] = (diff * np.abs(alt_percentile - ref_percentile)).argmin(axis=1) + 1
        max_variants_table['max_df_x_percdf_topic_alt'] = (diff * np.abs(alt_percentile - ref_percentile)).argmax(axis=1) + 1
        max_variants_table['max_df_x_percdf_topic'] = (np.abs(diff) * np.abs(alt_percentile - ref_percentile)).argmax(axis=1) + 1

        if len(shuf_variants_table) > 0:
            variants_table["max_percentile.pval"] = get_pvals(
                variants_table["max_percentile"].tolist(),
                shuf_max_percentile, tail="right"
            )
            variants_table['percentile_change.pval'] = get_pvals(
                variants_table["percentile_change"].tolist(),
                shuf_percentile_change, tail="both"
            )
            variants_table["abs_percentile_change.pval"] = get_pvals(
                variants_table["abs_percentile_change"].tolist(),
                shuf_abs_percentile_change, tail="right"
            )
            variants_table["logfc_x_max_percentile.pval"] = get_pvals(
                variants_table["logfc_x_max_percentile"].tolist(),
                shuf_logfc_max_percentile, tail="both"
            )
            variants_table["abs_logfc_x_max_percentile.pval"] = get_pvals(
                variants_table["abs_logfc_x_max_percentile"].tolist(),
                shuf_abs_logfc_max_percentile, tail="right"
            )
            variants_table["diff_x_max_percentile.pval"] = get_pvals(
                variants_table["diff_x_max_percentile"].tolist(),
                shuf_diff_max_percentile, tail="right"
            )
            variants_table["abs_diff_x_max_percentile.pval"] = get_pvals(
                variants_table["abs_diff_x_max_percentile"].tolist(),
                shuf_abs_diff_max_percentile,
                tail="right"
            )

    print()
    print(variants_table.head())
    print("Output score table shape:", variants_table.shape)
    print()
    print(max_variants_table.head())
    print("Output max scores table shape:", max_variants_table.shape)
    print()

    variants_table.to_parquet(output_path.joinpath(Path("variant_scores.parquet")), engine="pyarrow", compression="snappy")
    max_variants_table.to_parquet(output_path.joinpath(Path("max_variant_scores.parquet")), engine="pyarrow", compression="snappy")
    print("Results saved in:", output_path)

    # Store predictions at variants
    if args.hdf5:
        with h5py.File(output_path.joinpath(Path("variant_predictions.h5")), 'w') as f:
            observed = f.create_group('observed')
            observed.create_dataset('ref_pred', data=ref_pred, compression='gzip', compression_opts=9)
            observed.create_dataset('alt_pred', data=alt_pred, compression='gzip', compression_opts=9)
            if len(shuf_variants_table) > 0:
                shuffled = f.create_group('shuffled')
                shuffled.create_dataset('shuf_ref_pred', data=shuf_ref_pred, compression='gzip', compression_opts=9)
                shuffled.create_dataset('shuf_alt_pred', data=shuf_alt_pred, compression='gzip', compression_opts=9)
                shuffled.create_dataset('shuf_logfc', data=shuf_logfc, compression='gzip', compression_opts=9)
                shuffled.create_dataset('shuf_abs_logfc', data=shuf_abs_logfc, compression='gzip', compression_opts=9)
                shuffled.create_dataset('shuf_diff', data=shuf_diff, compression='gzip', compression_opts=9)
                shuffled.create_dataset('shuf_abs_diff', data=shuf_abs_diff, compression='gzip', compression_opts=9)
                if args.peaks:
                    shuffled.create_dataset('shuf_max_percentile', data=shuf_max_percentile, compression='gzip', compression_opts=9)
                    shuffled.create_dataset('shuf_percentile_change', data=shuf_percentile_change, compression='gzip', compression_opts=9)
                    shuffled.create_dataset('shuf_abs_percentile_change', data=shuf_abs_percentile_change, compression='gzip', compression_opts=9)
                    shuffled.create_dataset('shuf_logfc_x_max_percentile', data=shuf_logfc_max_percentile, compression='gzip', compression_opts=9)
                    shuffled.create_dataset('shuf_abs_logfc_max_percentile', data=shuf_abs_logfc_max_percentile, compression='gzip', compression_opts=9)
                    shuffled.create_dataset('shuf_diff_x_max_percentile', data=shuf_diff_max_percentile, compression='gzip', compression_opts=9)
                    shuffled.create_dataset('shuf_abs_diff_x_max_percentile', data=shuf_abs_diff_max_percentile, compression='gzip', compression_opts=9)

    print("DONE")
    print()


if __name__ == "__main__":
    args = fetch_variant_args()
    main(args)
