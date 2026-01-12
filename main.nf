#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include local modules
include { make_csv } from './modules/local/make_csv.nf'
include { merge_fastq } from './modules/local/merge_fastq.nf'
include { nanoplot } from './modules/local/nanoplot.nf'
include { porechop } from './modules/local/porechop.nf'
include { minimap2 } from './modules/local/minimap2.nf'
include { splitbam } from './modules/local/splitbam.nf'
include { medaka } from './modules/local/medaka.nf'
include { multiqc } from './modules/local/multiqc.nf'
include { kraken2 } from './modules/local/kraken2.nf'
include { krona_kraken } from './modules/local/krona_kraken.nf'
include { make_report } from './modules/local/make_report.nf'
include { abricate } from './modules/local/abricate.nf'
include { mlst } from './modules/local/mlst.nf'
include { client_report } from './modules/local/client_report.nf'
include { bracken } from './modules/local/bracken.nf'

// Include dehost subworkflow
include { DEHOST } from './subworkflows/dehost.nf'


workflow {
	data=Channel
	.fromPath(params.input)
	merge_fastq(make_csv(data).splitCsv(header:true).map { row-> tuple(row.SampleName,row.SamplePath)})
	reference=file("${baseDir}/Bovreproseq_reference.fasta")
	primerbed=file("${baseDir}/Bovreproseq_primer.bed")
	//read statistics
	nanoplot(merge_fastq.out.unfilt_reads)	
		//trim barcodes and adapter sequences
	if (params.trim_barcodes){
		porechop(merge_fastq.out.reads)
		
		// Dehost after porechop (if enabled)
		if (params.dehost) {
			host_ref = file(params.host_db)
			DEHOST(porechop.out, host_ref)
			reads_for_alignment = DEHOST.out.dehosted_reads
		} else {
			reads_for_alignment = porechop.out
		}
		
		minimap2(reference, reads_for_alignment)
		 
	} else {
		// Dehost after merge_fastq
		if (params.dehost) {
			host_ref = file(params.host_db)
			DEHOST(merge_fastq.out.reads, host_ref)
			reads_for_alignment = DEHOST.out.dehosted_reads
		} else {
			reads_for_alignment = merge_fastq.out.reads
		}
		
		minimap2(reference, reads_for_alignment)
	}
	

	

	kraken2(reads_for_alignment, params.kraken_db)
	bracken(kraken2.out, params.kraken_db)
	
	
	// create consensus
	splitbam(minimap2.out,primerbed,params.read_count_threshold,params.consensus_mode,params.qscore)


	//medaka polishing

		// Pair fastq file and consensus files by sample name
	paired_fastq_consensus = reads_for_alignment.join(splitbam.out.consensus)
    	.map { sample, fastq, consensus -> tuple(sample, fastq, consensus)}
	
	
	medaka (paired_fastq_consensus)


	krona_kraken(kraken2.out.kraken2_raw.collect())
	
	// qc report using split bam out put
	stats=splitbam.out.unfilt_stats
	idxstats=splitbam.out.idxstats
	nanoqc=nanoplot.out.stats_ufilt	
	dehosted_stats=DEHOST.out.stats
	multiqc(nanoqc.mix(stats,idxstats,dehosted_stats).collect())
	
	// abricate 
	dbdir=("${baseDir}/Bovreproseq_db")
	targetlist=("${baseDir}/Bovreproseq_targetlist.txt")
	abricate(medaka.out.consensus,dbdir,targetlist)
	
	//tax=("${baseDir}/taxdb")
	//blast_cons(splitbam.out.consensus,tax,db1)

	mlst(medaka.out.consensus)
	//generate report
	rmd_file=file("${baseDir}/Bovreproseq_tabbed.Rmd")
	make_report(make_csv.out,krona_kraken.out.raw,splitbam.out.mapped.collect(),abricate.out.abricate.collect(),abricate.out.results.collect(),rmd_file,mlst.out.collect(),medaka.out.cons_only.collect(),bracken.out.collect())
	
	
}
