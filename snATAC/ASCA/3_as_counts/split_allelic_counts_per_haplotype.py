import pandas as pd
import argparse

def main():
    # Parse command-line arguments
    parser = argparse.ArgumentParser(description="Split counts file into haplotype-specific counts.")
    parser.add_argument('--counts', type=str, required=True, help='Path to counts file (input).')
    parser.add_argument('--hap1', type=str, required=True, help='Path to output file for haplotype 1.')
    parser.add_argument('--hap2', type=str, required=True, help='Path to output file for haplotype 2.')
    args = parser.parse_args()

    counts_path = args.counts
    output_path1 = args.hap1
    output_path2 = args.hap2

    # Read counts file (first three rows are comments)
    counts = pd.read_csv(counts_path, sep='\t', skiprows=3)

    # Filter for heterozygous sites and split haplotype column
    counts = counts[(counts['haplotype'] == "0|1") | (counts['haplotype'] == "1|0")]
    counts[['hap1', 'hap2']] = counts['haplotype'].str.split('|', expand=True)

    # add end postion (always 1 bp) and sort
    counts["end"] = counts.position + 1
    counts = counts.sort_values(by=["#chrom", "position", "end"])

    # Compute counts for each haplotype
    counts["hap1_count"] = counts.apply(
        lambda row: row["ref_as_count"] if row["hap1"] == "0" else row["alt_as_count"], axis=1
    )
    counts["hap2_count"] = counts.apply(
        lambda row: row["ref_as_count"] if row["hap2"] == "0" else row["alt_as_count"], axis=1
    )

    # Save haplotype-specific counts to separate files
    counts[["#chrom", "position", "end", "hap1_count"]].to_csv(output_path1, sep='\t', index=False, header=False)
    counts[["#chrom", "position", "end", "hap2_count"]].to_csv(output_path2, sep='\t', index=False, header=False)

if __name__ == "__main__":
    main()