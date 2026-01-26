#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process abricate{
	publishDir "${params.out_dir}/abricate/",mode:"copy"
	label "low"
	input:
	tuple val(SampleName),path(consensus)
	path(dbdir)
	path (targetlist)
	output:
	path("${SampleName}_abricate.csv"),emit:abricate
	path("${SampleName}_results.csv"),emit:results
	script:
	"""
	abricate --datadir ${dbdir} --db Bovreproseq -minid 80  -mincov 85 --quiet ${consensus} 1> ${SampleName}_abricate.csv
	sed -i "s/_consensus//g" "${SampleName}_abricate.csv"
	interpret_results.sh ${SampleName} ${targetlist}
	"""
	
}
