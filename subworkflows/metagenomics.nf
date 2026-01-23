#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { kraken2 } from '../modules/local/kraken2.nf'
include { bracken } from '../modules/local/bracken.nf'
include { krona_kraken } from '../modules/local/krona_kraken.nf'
include { medaka } from '../modules/local/medaka.nf'


process extract_reads {
    publishDir "${params.out_dir}/binned_reads", mode: 'copy'
    label "high"

    input:
    tuple val(SampleName), path(fastq), path(bracken), path(kraken_output), path(kraken_report)

    output:
    tuple val(SampleName), path("*.fastq"), path("${SampleName}_taxids.txt"), path("${SampleName}_taxinfo.txt")

    script:
    """
    # Get TaxIDs and names from Bracken output (skip header, get columns 2 and 1)
    tail -n +2 ${bracken} | awk -F'\t' '{print \$2"\t"\$1}' > ${SampleName}_taxinfo.txt
    
    # Extract just taxIDs for the loop
    tail -n +2 ${bracken} | cut -f 2 > ${SampleName}_taxids.txt
    
    while read -r taxid; do
        extract_kraken_reads.py \\
            -k ${kraken_output} \\
            -s ${fastq} \\
            -o ${SampleName}_\${taxid}.fastq \\
            -t \${taxid} \\
            --include-children \\
            --report ${kraken_report}
    done < ${SampleName}_taxids.txt
    """
}

process megahit {
    publishDir "${params.out_dir}/megahit_metagenomes/", mode: "copy"
    label "high"

    input:
    tuple val(SampleName), path(fastqs), path(taxids), path(taxinfo)

    output:
    tuple val(SampleName), path("${SampleName}_metagenomes.fasta")

    script:
    """
    touch ${SampleName}_metagenomes.fasta
    
    # Loop through taxids using cat - create combined assembly file
    for taxid in \$(cat ${taxids}); do
      
        megahit -r "${SampleName}_\${taxid}.fastq" -o ${SampleName}_\${taxid}_assembly --k-list 41,61,81,99 --prune-level 1 --min-contig-len 200

        # Process output - append to combined assembly file
        if [ -s "${SampleName}_\${taxid}_assembly/final.contigs.fa" ]; then
            # Rename headers to include taxID and append
             sed -E "s/>k[0-9]+_/>${SampleName}_\${taxid}_contig_/g" "${SampleName}_\${taxid}_assembly/final.contigs.fa" >> ${SampleName}_metagenomes.fasta
        fi
    done
    
    # If no contigs were assembled, add placeholder
    if [ ! -s ${SampleName}_metagenomes.fasta ]; then
        echo ">${SampleName}_no_assembly" > ${SampleName}_metagenomes.fasta
        echo "NNNN" >> ${SampleName}_metagenomes.fasta
    fi
    """
}

process refseq_masher {
    publishDir "${params.out_dir}/refseq_masher/", mode: "copy"
    label "high"

    input:
    tuple val(SampleName), path(assembly)

    output:
    path("${SampleName}_refseqmasher.tsv")

    script:
    """
    # Check if assembly file has content and valid FASTA records
    if [ -s ${assembly} ] && grep -q "^>" ${assembly}; then
        refseq_masher -vv matches --top-n-results 50 ${assembly} > ${SampleName}_refseqmasher.tsv
    else
        # Create empty result file with header if no assembly
        echo -e "sample\ttop_taxonomy_name\tdistance\tpvalue\tmatching\tfull_taxonomy\ttaxonomic_subspecies\ttaxonomic_species\ttaxonomic_genus\ttaxonomic_family\ttaxonomic_order\ttaxonomic_class\ttaxonomic_phylum\ttaxonomic_superkingdom\tsubspecies\tspecies\tgenus\tfamily\torder\tclass\tphylum\tsuperkingdom" > ${SampleName}_refseqmasher.tsv
        echo -e "${SampleName}\tNO_ASSEMBLY\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A\tN/A" >> ${SampleName}_refseqmasher.tsv
    fi
    """
}

//performs blast of the consensus sequences
process blast_cons {
	publishDir "${params.out_dir}/blast/",mode:"copy"
	containerOptions "-v ${params.blastdb_path}:${params.blastdb_path}"

	label "high"
	input:
	tuple val(SampleName),path(consensus)
	path (blastdb_path)
	val(blastdb_name)

	output:
	tuple val(SampleName), path("${SampleName}_report_blast.tsv"), emit: blast_report
	tuple val(SampleName), path("${SampleName}_blast.tsv"), emit: blast_raw
	tuple val(SampleName), path ("${SampleName}_report_blast_best.tsv"),emit: blast_best
	script:
	"""
	

	# Run BLAST
	blastn -task megablast -perc_identity 95 -db ${blastdb_path}/${blastdb_name} -query ${consensus} -out ${SampleName}_blast.tsv -outfmt "7 qseqid sseqid length qcovs pident evalue staxids ssciname scomnames stitle" -max_target_seqs 5

	# Build report files from BLAST output
	make_blast_report.sh "${SampleName}"

    """

}

workflow METAGENOMICS {
    take:
	reads // Channel: tuple val(SampleName), path(fastq)
	kraken_db  // Path: kraken database
    blastdb_path // Path to blast database
    blastdb_name // Name of blast database
	
	main:
	// Run kraken2
	kraken2(reads, kraken_db)
	
	// Prepare bracken input by joining kraken outputs
	bracken_input = kraken2.out.kraken_output
	    .join(kraken2.out.kraken_report)
	    .map { sample, krakenraw, krakenReport -> tuple(sample, krakenraw, krakenReport) }
	
	// Run bracken
	bracken(bracken_input, kraken_db)
	
	// Prepare extract_reads input by joining all required channels
	extract_reads_input = reads
	    .join(kraken2.out.kraken_output)
	    .join(kraken2.out.kraken_report)
	    .join(bracken.out)
	    .map { sample, fastq, krakenCsv, krakenReport, brackenTsv -> 
	        tuple(sample, fastq, brackenTsv, krakenCsv, krakenReport) 
	    }
	
	extract_reads(extract_reads_input)
	megahit(extract_reads.out)
    //medaka_input = reads.join(megahit.out).map{sample,fastq,draft_fasta -> tuple (sample,fastq,draft_fasta)}
   // medaka(medaka_input)
    krona_kraken(kraken2.out.kraken_report.map{ sample, file -> file }.collect())
    //refseq_masher(megahit.out)
    blast_cons(megahit.out,blastdb_path,blastdb_name)
	
	emit:
	metagenomes = megahit.out.map{sample,consensus -> consensus}
	kraken_output = kraken2.out.kraken_output
	kraken_report = kraken2.out.kraken_report
	bracken_output = bracken.out
    krona_output = krona_kraken.out
    //refmash=refseq_masher.out
    blast_report = blast_cons.out.blast_report
    blast_raw = blast_cons.out.blast_raw
    blast_best=blast_cons.out.blast_best
	
}