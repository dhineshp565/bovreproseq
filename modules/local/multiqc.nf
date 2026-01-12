#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process multiqc {
	publishDir "${params.out_dir}/",mode:"copy"
	label "low"
	input:
	path '*'
	output:
	file ("multiqc_report.html")
	file ("multiqc_data")
	script:
	"""
	multiqc .
	"""
}
