#!/usr/bin/env nextflow
nextflow.enable.dsl=2

// Include local modules
include { make_csv } from './modules/local/make_csv.nf'
include { merge_fastq } from './modules/local/merge_fastq.nf'
include { nanoplot } from './modules/local/nanoplot.nf'
include { porechop } from './modules/local/porechop.nf'
//include { minimap2 } from './modules/local/minimap2.nf'
//include { splitbam } from './modules/local/splitbam.nf'
//include { medaka } from './modules/local/medaka.nf'
include { multiqc } from './modules/local/multiqc.nf'
include { kraken2 } from './modules/local/kraken2.nf'
include { krona_kraken } from './modules/local/krona_kraken.nf'
include { make_report } from './modules/local/make_report.nf'
include { make_metagenomics_report } from './modules/local/make_metagenomics_report.nf'
include { abricate } from './modules/local/abricate.nf'
include { mlst } from './modules/local/mlst.nf'
include { client_report } from './modules/local/client_report.nf'
include { bracken } from './modules/local/bracken.nf'

// Include dehost subworkflow
include { DEHOST } from './subworkflows/dehost.nf'
include { AMPLICONS } from './subworkflows/amplicons.nf'
include { METAGENOMICS } from './subworkflows/metagenomics.nf'
include { QCREADS } from './subworkflows/qcreads.nf'


workflow {
	// Merge fastq files and qc filtering
    QCREADS(params.input, params.qscore, params.trim_barcodes)

    // Dehost logic
    if (params.dehost) {
        host_ref = file(params.host_db)
        DEHOST(QCREADS.out.reads, host_ref)
        reads_for_alignment = DEHOST.out.dehosted_reads
    } else {
        reads_for_alignment = QCREADS.out.reads
    }

	reference=file("${baseDir}/Bovreproseq_reference.fasta")
	primerbed=file("${baseDir}/Bovreproseq_primer.bed")
    
    // Bovreproseq path
    AMPLICONS(reads_for_alignment,reference,primerbed,params.read_count_threshold,params.consensus_mode,params.qscore)
	
	if (params.dehost) {
		stats=DEHOST.out.stats
	} else {
		stats=AMPLICONS.out.stats
	} 
	//kraken2(reads_for_alignment, params.kraken_db)
	//bracken(kraken2.out, params.kraken_db)
	//krona_kraken(kraken2.out.kraken2_raw.collect())
	
	// qc report using split bam out put
	idxstats=AMPLICONS.out.idxstats
	nanoqc=QCREADS.out.read_stats	
	multiqc(nanoqc.mix(idxstats,stats).collect())
	
	// abricate 
	dbdir=("${baseDir}/Bovreproseq_db")
	targetlist=("${baseDir}/Bovreproseq_targetlist.txt")
	abricate(AMPLICONS.out.consensus,dbdir,targetlist)
	
	//tax=("${baseDir}/taxdb")
	//blast_cons(splitbam.out.consensus,tax,db1)

	mlst(AMPLICONS.out.consensus)
	IGVREPORTS(reference,AMPLICONS.out.bam.collect())

	METAGENOMICS (reads_for_alignment,params.kraken_db,params.blastdb_path,params.blastdb_name)
	//generate report
	rmd_file=file("${baseDir}/Bovreproseq_tabbed.Rmd")
	meta_bracken=METAGENOMICS.out.bracken_output.map{ sample, file -> file }.collect()
	meta_blast= METAGENOMICS.out.blast_best.map {sample, file -> file}.collect()
	make_report(QCREADS.out.csv,METAGENOMICS.out.krona_output,AMPLICONS.out.mapped.collect(),abricate.out.abricate.collect(),abricate.out.results.collect(),rmd_file,mlst.out.collect(),AMPLICONS.out.cons_only.collect(),meta_bracken,meta_blast)
	
	//generate metagenomics report
	rmd_meta_file=file("${baseDir}/ont_metagenomics.Rmd")
    make_metagenomics_report(QCREADS.out.csv,METAGENOMICS.out.krona_output,rmd_meta_file,meta_bracken,meta_blast)
	
	
	
}
