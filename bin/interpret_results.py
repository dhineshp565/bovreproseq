#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import csv

def interpret_results(samplename,pathogen_target_map,abricate_file,mappedreads_file,output_file):
    with open (pathogen_target_map, 'r') as pathogen_map_file:
        next(pathogen_map_file)  # Skip header line
        pathogen_dict = {}
        for line in pathogen_map_file:
            line = line.strip().split('\t')
            pathogen = line[0]
            targets = line[1:]
            pathogen_dict[pathogen] = targets
        

    with open (abricate_file, 'r') as abricate_file:
            reader = csv.DictReader(abricate_file, delimiter='\t')
            pathogens_in_sample=[]
            for col in reader:
                if col['RESISTANCE'] not in pathogens_in_sample:
                        pathogens_in_sample.append(col['RESISTANCE'])
            print (pathogens_in_sample)

    with open (mappedreads_file, 'r') as mappedreads_file:
        next(mappedreads_file)  # Skip header line
        read_counts = {}
        for line in mappedreads_file:
            line = line.strip().split(' ')
            read_counts[line[0]]=int(line[2])
        print (read_counts)

    # Check for PCR inhibition
    detected_map_pathogens = [p for p in pathogens_in_sample if p in pathogen_dict]
    if not detected_map_pathogens:
        global_qc_status = "Check for PCR inhibition"
    else:
        global_qc_status = "OK"

    with open (output_file, 'w') as results_file:
        results_file.write('SampleName\tPathogen\tQC_Status\tReadCount\tResult\n')
        total_reads={}
        for pathogen, targets in pathogen_dict.items():
            if pathogen in pathogens_in_sample:
                PCR_target1=pathogen_dict[pathogen][0]
                PCR_target2=pathogen_dict[pathogen][1]
                total_count=read_counts.get(PCR_target1,0)+read_counts.get(PCR_target2,0)
                total_reads[pathogen]=total_count
                results_file.write(f'{samplename}\t{pathogen}\t{global_qc_status}\t{total_count}\tPositive\n')
            else :
                total_reads[pathogen]=0
                results_file.write(f'{samplename}\t{pathogen}\t{global_qc_status}\t0\tNegative\n')
        
    results_file.close()

if __name__ == "__main__":
    import sys
    samplename = sys.argv[1]
    pathogen_target_map = sys.argv[2]
    abricate_file = sys.argv[3]
    mappedreads_file = sys.argv[4]
    output_file = sys.argv[5]
    interpret_results(samplename,pathogen_target_map,abricate_file,mappedreads_file,output_file)