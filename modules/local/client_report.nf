#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process client_report {
	publishDir "${params.out_dir}/client_report/",mode:"copy"
	label "low"
	input:
	path (abricate)
	output:
	path("Bovreproseq_results_report.html")
	script:
	"""
	"""

}
