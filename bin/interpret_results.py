#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Script to interpret pathogen detection results.

This script processes outputs from Abricate and mapping results from samtools to determine
the presence or absence of pathogens in a given sample. It generates a summary
CSV report including QC status, read counts and abundance for each pathogen.
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
    # map pathogens to pcr targets by creating dictionary with pathogen name as key and targets as values
    with open(pathogen_target_map, 'r') as pathogen_map_file:
        next(pathogen_map_file)  # Skip header line
        pathogen_dict = {}
        for line in pathogen_map_file:
            line = line.strip().split('\t')
            pathogen = line[0]
            targets = line[1:]
            pathogen_dict[pathogen] = targets
        
   # creates a list of pathogens in the panel by looking at the resistance column of the abricate file
    with open(abricate_file, 'r') as abricate_f:
        reader = csv.DictReader(abricate_f, delimiter='\t')
        pathogens_in_sample = []
        for col in reader:
            # Use the 'RESISTANCE' column as the marker for detected
            # elements. Duplicates are skipped so we only have a unique list.
            if col['RESISTANCE'] not in pathogens_in_sample:
                pathogens_in_sample.append(col['RESISTANCE'])
        print(f"Pathogens in sample: {pathogens_in_sample}")

   # Converts the mapped reads file into a dictionary with targets as key and the read counts as values
    # New format includes: Amplicon_Name Size ReadCount Relative_Abundance
    with open(mappedreads_file, 'r') as mappedreads_f:
        next(mappedreads_f)  # Skip header line
        read_counts = {}  # Store {target: (count, abundance)}
        for line in mappedreads_f:
            line = line.strip().split(' ')
            # line[0] is target name, line[2] is read count, line[3] is relative abundance
            if len(line) > 3:
                # Store both count and relative abundance as a tuple
                read_counts[line[0]] = (int(line[2]), float(line[3]))
        print(f"Read counts: {read_counts}")

    # Check for PCR inhibition
    # If Abricate did not report any of the pathogens that we have targets
    # for (i.e. none of the pathogens in the sample appear in the mapping),
    # flag the sample to 'Check for PCR inhibition'. Otherwise QC is OK.
    detected_map_pathogens = [p for p in pathogens_in_sample if p in pathogen_dict]
    if not detected_map_pathogens:
        global_qc_status = "Check for PCR inhibition"
    else:
        global_qc_status = "OK"

    # Write output to CSV file
    with open(output_file, 'w') as results_file:
        # Write header for output CSV-like TSV that downstream tools expect
        results_file.write('Sample\tPathogen\tQC_Status\tReadCount\tAbundance(%)\tResult\n')
        total_reads = {}
        
        for pathogen, targets in pathogen_dict.items():
            if pathogen in pathogens_in_sample:
                # Calculate total reads for the two targets of the pathogen
                PCR_target1 = pathogen_dict[pathogen][0]
                PCR_target2 = pathogen_dict[pathogen][1]

                # Lookup read counts and abundance; defaulting to 0 if a target wasn't found
                # read_counts values are tuples of (count, abundance)
                count1, abund1 = read_counts.get(PCR_target1, (0, 0))
                count2, abund2 = read_counts.get(PCR_target2, (0, 0))
                # adding the counts and abundace for each amplicon
                total_count = count1 + count2
                total_abundance = abund1 + abund2

                total_reads[pathogen] = total_count
                results_file.write(f'{samplename}\t{pathogen}\t{global_qc_status}\t{total_count}\t{total_abundance}\tPositive\n')
            else:
                # Pathogen not found by Abricate; report as Negative with 0 reads
                total_reads[pathogen] = 0
                results_file.write(f'{samplename}\t{pathogen}\t{global_qc_status}\t0\t0\tNegative\n')
        
    
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
