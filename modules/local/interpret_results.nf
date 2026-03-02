#!/usr/bin/env nextflow

process interpret_results {
    publishDir "${params.out_dir}/bovreproseq_lims/", mode: "copy"
    label "low"

    input:
    tuple val(SampleName), path(mapped_reads),path(abricate_results)
    path(pathogentargetmap)

    output:
    path("${SampleName}_bovreproseq_results.tsv")

    script:
    """
    interpret_results.py ${SampleName} ${pathogentargetmap} ${abricate_results} ${mapped_reads} ${SampleName}_bovreproseq_results.tsv
    
    """
}

