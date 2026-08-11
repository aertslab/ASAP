#!/usr/bin/env python3
"""
define_blocks.py

Defines LD-based "linkage regions" around GP2 GWAS lead variants, using
plink2 --r2-phased output against 1000G, then lifts the regions to
hs1/T2T/chm13 and (optionally) intersects them with donor variants.
"""

import argparse
import gzip
import shutil
import subprocess
import sys
import zipfile
import urllib.request
from collections import Counter
from pathlib import Path
 
import numpy as np
import pandas as pd
from pyliftover import LiftOver


# --------------------------------------------------------------------------
# Region construction
# --------------------------------------------------------------------------

def load_leads(leads_path: Path) -> pd.DataFrame:
    leads = pd.read_csv(leads_path, sep="\t")
    print(f" [leads] Lead variants: {leads.shape[0]}")
    return leads


def load_r2(r2_path: Path) -> pd.DataFrame:
    r2 = pd.read_csv(r2_path, sep="\t")
    print(f" [r2] {r2.shape[0]} r2 matches between leads and 1KG variants")
    return r2


def build_raw_regions(leads: pd.DataFrame, r2: pd.DataFrame, r2_threshold: float,  save_path:str=None) -> pd.DataFrame:
    """
    Demarcate lead-variant centered regions for a single r2 threshold.
    Algorithm: Find the furthest variant 'left' and 'right' of the lead with an r2 to that lead
    being higher than the threshold.
    r2_threshold: float, e.g. 0.6 (rerun the script with a different value to compare thresholds).
    """
    ids_a = set(r2["ID_A"])
    linked = r2[r2["PHASED_R2"] > r2_threshold]
    regions = []

    leads_in_1kg = pd.read_csv('3_variants/leads_in_1kg.txt', header=None, names=['ID'])

    for lead_id in set(leads["ID"]):

        # If lead id NOT in subset of 1000G
        if lead_id not in list(leads_in_1kg['ID']):
            chr_ = leads[leads["ID"] == lead_id]["CHR"].values[0]
            lead = leads[leads["ID"] == lead_id]["BP"].values[0]
            regions.append({
                "chr": chr_,
                "lead_id": lead_id,
                "lead_bp": lead,
                "left_raw": lead,
                "right_raw": lead,
                "left_dist": 0,
                "right_dist": 0,
                "size": 0,
                "r2_threshold": r2_threshold,
                "1kg": False,
            })
            continue

        r2_block = linked[linked["ID_A"] == lead_id]
        linked_ids = set(r2_block["ID_A"])

        # Has no r2 with any variant at this threshold
        if lead_id not in linked_ids:
            chr_ = leads[leads["ID"] == lead_id]["CHR"].values[0]
            lead = leads[leads["ID"] == lead_id]["BP"].values[0]
            regions.append({
                "chr": chr_,
                "lead_id": lead_id,
                "lead_bp": lead,
                "left_raw": lead,
                "right_raw": lead,
                "left_dist": 0,
                "right_dist": 0,
                "size": 0,
                "r2_threshold": r2_threshold,
                "1kg": True,
            })
            continue

        # Has r2 with at leats one variant at this threshold
        chr_ = r2_block["#CHROM_A"].values[0]
        lead = r2_block["POS_A"].values[0]
        left_border = min([lead, min(r2_block["POS_B"])])
        right_border = max([lead, max(r2_block["POS_B"])])
        ld_left = max([0, r2_block["POS_A"].values[0] - left_border])
        ld_right = max([0, right_border - r2_block["POS_A"].values[0]])
        regions.append({
            "chr": chr_,
            "lead_id": lead_id,
            "lead_bp": lead,
            "left_raw": left_border,
            "right_raw": right_border,
            "left_dist": ld_left,
            "right_dist": ld_right,
            "size": ld_left + ld_right,
            "r2_threshold": r2_threshold,
            "1kg": True,
        })

    regions_raw = pd.DataFrame(regions)
    regions_raw["hg38id_short"] = (
        ["chr"] * len(regions_raw)
        + regions_raw["chr"].astype(str)
        + ["_"] * len(regions_raw)
        + regions_raw["lead_bp"].astype(str)
    )

    print(f" [regions] {regions_raw[regions_raw['1kg'] == False].shape[0]} leads not in 1KG")
    print(f" [regions] {regions_raw[regions_raw['size'] == 0].shape[0]} leads have no r2 at threshold {r2_threshold}")

    if save_path:
        print(f" [regions] Exporting {regions_raw.shape[0]} regions to {save_path}")
        regions_raw.to_csv(save_path, sep="\t", index=False)
    
    return regions_raw


def apply_buffer(regions_raw: pd.DataFrame, buffer_bp: int) -> pd.DataFrame:
    """Apply a lower bound of `buffer_bp` at each side of the lead variant."""
    print(f" [sizing] Applying a minimal half-size of {buffer_bp}...")
    regions_raw['left_bound'] = regions_raw.apply(lambda r: r['lead_bp'] - buffer_bp if r['lead_bp'] - r['left_raw'] < buffer_bp else r['left_raw'], axis=1)
    regions_raw['right_bound'] = regions_raw.apply(lambda r: r['lead_bp'] + buffer_bp if r['right_raw'] - r['lead_bp'] < buffer_bp else r['right_raw'], axis=1)
    regions_raw['left_dist'] = regions_raw.apply(lambda r: r['lead_bp'] - r['left_bound'], axis=1)
    regions_raw['right_dist'] = regions_raw.apply(lambda r: r['right_bound'] - r['lead_bp'], axis=1)
    regions_raw['size'] = regions_raw.apply(lambda r: r['right_bound'] - r['left_bound'], axis=1)
    return regions_raw


def add_nearest_gene(regions: pd.DataFrame, leads: pd.DataFrame) -> pd.DataFrame:
    print(f" [gene] Adding nearest gene according to GP2...")
    out = regions_raw_gene = regions.merge(leads[['ID','Nearest Gene']], left_on='lead_id', right_on='ID'
             ).rename(columns={'Nearest Gene':'nearest_gene'}).drop(columns=['r2_threshold','ID'])
    return out


# --------------------------------------------------------------------------
# Main
# --------------------------------------------------------------------------

def parse_args():
    parser = argparse.ArgumentParser(
        description="Define LD-based linkage regions around GWAS lead variants.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--leads", type=Path, default=Path("2_source/GP2_leads.tsv"),
                         help="Path to GP2 leads TSV")
    parser.add_argument("--r2", type=Path, default=Path("1_plink/linkage.vcor"),
                         help="Path to plink2 --r2-phased .vcor output")
    parser.add_argument("--out-dir", type=Path, default=Path("4_regions"),
                         help="Directory for region output files")
    parser.add_argument("--variants-dir", type=Path, default=Path("3_variants"),
                         help="Directory for donor variant output files")
    parser.add_argument("--donor-bcf", type=Path, default=Path("2_source/donor_vars.bcf"),
                         help="Path to donor variant BCF")
    parser.add_argument("--r2-threshold", type=float, default=0.6,
                         help="R-squared threshold used to define linkage regions "
                              "(e.g. 0.6 keeps links with PHASED_R2 > 0.6). Regions are "
                              "built for this one threshold; rerun with a different value "
                              "to compare thresholds.")
    parser.add_argument("--buffer", type=int, default=1500,
                         help="Lower bound (bp) applied at each side of the lead variant")
    parser.add_argument("--project-root", type=Path, default=Path("."),
                         help="Root directory under which 1_plink/2_source/3_variants/"
                              "4_regions/5_figs live (used with --setup)")
    return parser.parse_args()


def main():
    args = parse_args()

    print(args.leads)

    leads = load_leads(args.project_root / args.leads)
    r2 = load_r2(args.project_root / args.r2)

    regions_raw = build_raw_regions(leads, r2, args.r2_threshold, save_path = args.project_root / args.out_dir / "regions_0_raw.tsv")

    regions_buffered = apply_buffer(regions_raw, args.buffer)
    regions_with_gene = add_nearest_gene(regions_buffered, leads)
    regions_with_gene['chr'] = [f"chr{chr}" for chr in regions_with_gene['chr']]
    regions_with_gene.to_csv(args.project_root / args.out_dir / "regions.tsv", sep="\t", index=False)


if __name__ == "__main__":
    main()
