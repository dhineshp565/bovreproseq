#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process splitbam {
	publishDir "${params.out_dir}/splitbam",mode:"copy"
	label "medium"
	input:
	tuple val(SampleName),path(samfile)
    path (primerbed)
    val(readcount)
    val(consensus_mode)
    val (qscore)

	output:
	val(SampleName),emit:SampleName
	tuple val(SampleName), path("${SampleName}_mappedreads.txt"),emit:mapped
	path("${SampleName}_idxstats.txt"),emit:idxstats
	tuple val(SampleName),path("${SampleName}_consensus.fasta"),emit:consensus
	path("${SampleName}_consensus.fasta"),emit:(cons_only)
	path("${SampleName}_unfilt_stats.txt"),emit:unfilt_stats
	path ("${SampleName}_full_length_mappedreads.txt"),emit:full_reads
	tuple val(SampleName),path("${SampleName}.bam"),path("${SampleName}.bai"),emit:target_bam
	
	script:
	"""
	splitbam.sh ${SampleName} ${samfile} ${primerbed} ${readcount} ${consensus_mode} ${qscore}

	"""
}
