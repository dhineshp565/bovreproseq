#!/usr/bin/env nextflow

process bovreproseq_lims {
    publishDir "${params.out_dir}/bovreproseq_lims/", mode: "copy"
    label "high"

    input:
    tuple val(SampleName), path(mapped_reads),path(abricate_results)
    path(pathogentargetmap)

    output:
    path("${SampleName}_bovreproseq_lims_report.tsv")

    script:
    """
    interpret_results.py ${SampleName} ${pathogentargetmap} ${abricate_results} ${mapped_reads} ${SampleName}_bovreproseq_lims_report.tsv
    """
}