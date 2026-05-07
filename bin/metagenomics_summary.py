#!/usr/bin/env python3


import pandas as pd
import sys


def metagenomics_summary(
    bracken_tsv,
    blast_tsv,
    sample_name,
    blast_best_hits_out=None,
    summary_out=None,
):
    # Build default output names using the sample prefix when paths are not provided.
    if blast_best_hits_out is None:
        blast_best_hits_out = f"{sample_name}_blast_best_hits.tsv"
    if summary_out is None:
        summary_out = f"{sample_name}_kracken_blast_summary.tsv"

    # Load Bracken and BLAST reports.
    bracken_df = pd.read_csv(bracken_tsv, sep="\t")
    blast_df = pd.read_csv(blast_tsv, sep="\t")

    # Compute relative abundance per taxon from Bracken columns.
    bracken_df['Abundance (%)'] = (
        bracken_df['kraken_assigned_reads']
        / (bracken_df['new_est_reads'] / bracken_df['fraction_total_reads'])
    )

    # Extract taxonomy id and sample id encoded in BLAST query id.
    blast_df['kraken_taxid'] = blast_df['queryid'].str.split('_').str[1].astype(int)
    blast_df['sampleid'] = blast_df['queryid'].str.split('_').str[0].astype(str)

    # Keep only high-confidence BLAST alignments and de-duplicate top hits.
    blast_best_hits = (
        blast_df[
            (blast_df['alignment length'] > 100)
            & (blast_df['%identity'] > 80)
            & (blast_df['evalue'] < 1e-5)
            & (blast_df['query_coverage'] > 50)
        ]
        .sort_values(by='bitscore', ascending=False)
        .drop_duplicates(subset=['staxids'])
        .drop_duplicates(subset=['queryid'])
    )
    blast_best_hits.sort_values(by='stitle', inplace=True)
    blast_best_hits.to_csv(blast_best_hits_out, sep="\t", index=False)

    # Join Bracken abundance with BLAST best hits using taxonomy id.
    bracken_select = bracken_df[['name', 'taxonomy_id', 'kraken_assigned_reads', 'Abundance (%)']]
    merged_df = pd.merge(
        bracken_select,
        blast_best_hits,
        left_on='taxonomy_id',
        right_on='kraken_taxid',
        how='inner',
    )
    # Trim species title to genus + species for cleaner report output.
    merged_df['stitle'] = merged_df['stitle'].str.split(' ').str[0:2].str.join(' ')


    print(merged_df)

    # Aggregate species hits per sample/genus and format final summary columns.
    summary = (
        merged_df.groupby(['sampleid', 'name', 'kraken_assigned_reads', 'Abundance (%)'])
        .agg(species_list=('stitle', ', '.join))
        .reset_index()
    )
    summary['Abundance (%)'] = summary['Abundance (%)'].apply(lambda x: f"{x*100:.4f}%")
    summary = summary.rename(
        columns={
            'sampleid': 'SampleID',
            'name': 'kraken2-Genus_level',
            'kraken_assigned_reads': 'Read_Count',
            'species_list': 'BLAST_hits_Species',
        }
    )
    summary = summary.sort_values(by='Read_Count', ascending=False)
   
    print(summary)
    summary.to_csv(summary_out, sep="\t", index=False)




def main():
    if len(sys.argv) < 4:
        print("Usage: metagenomics_summary.py <bracken_tsv> <blast_tsv> <sample_name>")
        sys.exit(1)

    # Parse positional CLI arguments and run the report generator.
    bracken_tsv = sys.argv[1]
    blast_tsv = sys.argv[2]
    sample_name = sys.argv[3]
    metagenomics_summary(bracken_tsv, blast_tsv, sample_name)

if __name__ == "__main__":
    main()