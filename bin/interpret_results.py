#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to interpret pathogen detection results.

This script processes outputs from Abricate and mapping results from samtools to determine
the presence or absence of pathogens in a given sample. It generates a summary
TSV report including QC status, read counts and abundance for each pathogen.
"""

import csv
import sys


def interpret_results(samplename, pathogen_target_map, abricate_file, mappedreads_file, output_file):
    """
    Interpret results from pathogen detection pipeline.

    Args:
        samplename (str): Name of the sample being processed.
        pathogen_target_map (str): Path to TSV file mapping pathogens to their PCR target markers.
        abricate_file (str): Path to Abricate output file (tab-delimited).
        mappedreads_file (str): Path to file containing mapped read counts (space-delimited).
        output_file (str): Path where the output TSV will be written.
    """
    # ---- SECTION 1: Read pathogen target map ----
    # Build mappings: pathogen -> PCR targets, pathogen -> LIMS display name,
    # and identify LIMS names shared by multiple pathogens (e.g. "Leptospira species")
    with open(pathogen_target_map, 'r') as pathogen_map_file:
        next(pathogen_map_file)
        pathogen_dict = {}           # {pathogen_name: [PCR targets]}
        lims_pathogen_dict = {}      # {pathogen_name: LIMS display name}
        monospp_lims_names = set()
        multispp_lims_names = set()  # LIMS names shared by multiple pathogens
        for line in pathogen_map_file:
            line = line.strip().split('\t')
            pathogen = line[0]
            pcr_targets = []
            for t in line[1:-1]:
                if t != 'NA':
                    pcr_targets.append(t)
            pathogen_dict[pathogen] = pcr_targets
            lims_name = line[-1]
            lims_pathogen_dict[pathogen] = lims_name
            if lims_name in monospp_lims_names:
                multispp_lims_names.add(lims_name)
            else:
                monospp_lims_names.add(lims_name)

    # ---- SECTION 2: Read abricate results ----
    # Collect unique pathogens detected via the RESISTANCE column
    with open(abricate_file, 'r') as abricate_f:
        reader = csv.DictReader(abricate_f, delimiter='\t')
        pathogens_in_sample = []
        for col in reader:
            if col['RESISTANCE'] not in pathogens_in_sample:
                pathogens_in_sample.append(col['RESISTANCE'])

    # ---- SECTION 3: Read mapped reads ----
    # Format: Amplicon_Name Size ReadCount Relative_Abundance (space-delimited)
    with open(mappedreads_file, 'r') as mappedreads_f:
        next(mappedreads_f)
        read_counts = {}  # {pcr_target: (count, abundance)}
        for line in mappedreads_f:
            line = line.strip().split(' ')
            read_counts[line[0]] = (int(line[2]), float(line[3]))

        # If CFV_parA is present, remove generic C. fetus so subspecies are reported separately
        if 'CFV_parA' in read_counts:
            pathogens_in_sample.remove('Campylobacter_fetus')

    # ---- SECTION 4: Check for PCR inhibition ----
    # Flag if no abricate pathogens matched the target map (including internal control)
    detected_map_pathogens = []
    for p in pathogens_in_sample:
        if p in pathogen_dict:
            detected_map_pathogens.append(p)

    if detected_map_pathogens:
        global_qc_status = "OK"
    else:
        global_qc_status = "Check for PCR inhibition"

    # ---- SECTION 5: Write output file ----
    with open(output_file, 'w') as results_file:
        results_file.write('Sample\tPathogen\tQC_Status\tReadCount\tAbundance(%)\tResult\tComment\n')

        lims_results = {}  # One output line per LIMS name; Positive overwrites Negative
        for pathogen, targets in pathogen_dict.items():
            lims_pathogen_name = lims_pathogen_dict.get(pathogen, pathogen)
            if pathogen in pathogens_in_sample:
                total_count = 0
                total_abundance = 0
                for target in targets:
                    count, abund = read_counts.get(target, (0, 0))
                    total_count += count
                    total_abundance += abund

                # When multiple pathogens share a LIMS name, add the specific name as a comment
                if lims_pathogen_name in multispp_lims_names:
                    comment = pathogen
                else:
                    comment = ''
                output_line = f'{samplename}\t{lims_pathogen_name}\t{global_qc_status}\t{total_count}\t{total_abundance}\tPositive\t{comment}\n'
                lims_results[lims_pathogen_name] = output_line
            else:
                output_line = f'{samplename}\t{lims_pathogen_name}\t{global_qc_status}\t0\t0\tNegative\t\n'
                if lims_pathogen_name not in lims_results:
                    lims_results[lims_pathogen_name] = output_line

        for output_line in lims_results.values():
            results_file.write(output_line)


if __name__ == "__main__":
    if len(sys.argv) < 6:
        print("Usage: python interpret_results.py <samplename> <pathogen_target_map> <abricate_file> <mappedreads_file> <output_file>")
        sys.exit(1)

    samplename = sys.argv[1]
    pathogen_target_map = sys.argv[2]
    abricate_file = sys.argv[3]
    mappedreads_file = sys.argv[4]
    output_file = sys.argv[5]
    interpret_results(samplename, pathogen_target_map, abricate_file, mappedreads_file, output_file)
