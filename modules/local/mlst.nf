#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process mlst {
	publishDir "${params.out_dir}/mlst/",mode:"copy"
	label "low"
	input:
	tuple val(SampleName),path(consensus)
	output:
	path("${SampleName}_MLST_results.csv")
	script:
	"""
	mlst --legacy --scheme campylobacter_nonjejuni_9 ${consensus} > ${SampleName}_MLST.csv
	campmlst.sh ${SampleName}_MLST.csv
	"""
}
