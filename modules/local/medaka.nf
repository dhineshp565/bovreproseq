#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process medaka {
	publishDir "${params.out_dir}/medaka",mode:"copy"
	label "high"
	input:
	tuple val(SampleName),path(SamplePath),path(consensus)
	output:
	tuple val(SampleName),path("${SampleName}_medaka_consensus.fasta"),emit:consensus
	path("${SampleName}_medaka_consensus.fasta"),emit: cons_only
	script:
	"""
	
	if [ \$(wc -l < "${consensus}" ) -gt 1 ]
		then
		medaka_consensus -i ${SamplePath} -d ${consensus} -o ${SampleName}_medaka_consensus
		mv ${SampleName}_medaka_consensus/consensus.fasta ${SampleName}_medaka_consensus.fasta
	else 
		echo ">${SampleName} No consensus sequences to polish" > ${SampleName}_medaka_consensus.fasta
	fi
	"""

}
