#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to interpret pathogen detection results.

This script processes outputs from Abricate and mapping results to determine
the presence or absence of pathogens in a given sample. It generates a summary
CSV report including QC status and read counts.
"""

import csv
import sys

def interpret_results(samplename, pathogen_target_map, abricate_file, mappedreads_file, output_file):
    """
    Interpret results from pathogen detection pipeline.

    Args:
        samplename (str): Name of the sample being processed.
        pathogen_target_map (str): Path to TSV file mapping pathogens to their target markers.
        abricate_file (str): Path to Abricate output file (tab-delimited).
        mappedreads_file (str): Path to file containing mapped read counts.
        output_file (str): Path where the output CSV will be written.
    """

    # Read the pathogen to target mapping file
    # Defines which targets correspond to which pathogen
    with open(pathogen_target_map, 'r') as pathogen_map_file:
        next(pathogen_map_file)  # Skip header line
        pathogen_dict = {}
        for line in pathogen_map_file:
            line = line.strip().split('\t')
            pathogen = line[0]
            targets = line[1:]
            pathogen_dict[pathogen] = targets
        
    # Read Abricate output to find detected resistance genes/pathogens
    with open(abricate_file, 'r') as abricate_f:
        reader = csv.DictReader(abricate_f, delimiter='\t')
        pathogens_in_sample = []
        for col in reader:
            if col['RESISTANCE'] not in pathogens_in_sample:
                pathogens_in_sample.append(col['RESISTANCE'])
        print(f"Pathogens in sample: {pathogens_in_sample}")

    # Read mapped read counts for targets
    with open(mappedreads_file, 'r') as mappedreads_f:
        next(mappedreads_f)  # Skip header line
        read_counts = {}
        for line in mappedreads_f:
            line = line.strip().split(' ')
            # line[0] is target name, line[2] is count
            if len(line) > 2:
                read_counts[line[0]] = int(line[2])
        print(f"Read counts: {read_counts}")

    # Check for PCR inhibition
    # If none of the mapped pathogens are found in the sample, flag for inhibition check
    detected_map_pathogens = [p for p in pathogens_in_sample if p in pathogen_dict]
    if not detected_map_pathogens:
        global_qc_status = "Check for PCR inhibition"
    else:
        global_qc_status = "OK"

    # Write output to CSV file
    with open(output_file, 'w') as results_file:
        results_file.write('SampleName\tPathogen\tQC_Status\tReadCount\tResult\n')
        total_reads = {}
        
        for pathogen, targets in pathogen_dict.items():
            if pathogen in pathogens_in_sample:
                # Calculate total reads for the first two targets of the pathogen
                PCR_target1 = pathogen_dict[pathogen][0]
                PCR_target2 = pathogen_dict[pathogen][1]
                
                count1 = read_counts.get(PCR_target1, 0)
                count2 = read_counts.get(PCR_target2, 0)
                total_count = count1 + count2
                
                total_reads[pathogen] = total_count
                results_file.write(f'{samplename}\t{pathogen}\t{global_qc_status}\t{total_count}\tPositive\n')
            else:
                # Pathogen not found
                total_reads[pathogen] = 0
                results_file.write(f'{samplename}\t{pathogen}\t{global_qc_status}\t0\tNegative\n')
        
    # File is automatically closed by context manager, but keeping explicit close if desired (removed here as repetitive)

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
