#!/usr/bin/env nextflow


process bovreproseq_lims {
     publishDir "${params.out_dir}/bovreproseq_lims/", mode: "copy"
    label "low"

    input:
    path(software_version)
    path(result_files)

    output:
    path("bovreproseq_LIMSfile_*.tsv")

    script:
    """
    date=\$(date '+%Y-%m-%d_%H-%M-%S')
    
    awk 'FNR==1 && NR!=1 { while (/^Sample/) getline; } 1 {print}' ${result_files} > bovreproseq_LIMSfile.tsv

    sed -i 's/_/ /g' bovreproseq_LIMSfile.tsv

    cat ${software_version} "bovreproseq_LIMSfile.tsv" > bovreproseq_LIMSfile_\${date}.tsv
    
    """
}
