#!/usr/bin/env nextflow

process bracken {
publishDir "${params.out_dir}/bracken/",mode:"copy"
label "high"

input:
tuple val(SampleName), path(kraken_output),path(kraken_report)
path(kraken_db)

output:
tuple val(SampleName), path("${SampleName}_bracken.tsv")

script:
"""
set +e  # Don't exit on error, we'll handle it manually

    # 1. Efficiently find the deepest taxonomic level with >= 10 reads
    # Prioritizes Species (S), then Genus (G), down to Domain (D)
    BRACKEN_LEVEL=\$(awk '\$2 >= 10 && \$3 ~ /^[SKOFPCKPD]\$/ {levels[\$3]++} 
        END { 
            split("S G F O C P K D", order, " "); 
            for (i in order) if (levels[order[i]]) { print order[i]; exit } 
            print "S" 
        }' ${kraken_report})

    # 2. Run Bracken; use || to handle failure gracefully
    bracken -d ${kraken_db} -i ${kraken_report} -o bracken_raw.tsv -r 150 -l \$BRACKEN_LEVEL || \
    echo -e "name\\ttaxonomy_id\\ttaxonomy_lvl\\tkraken_assigned_reads\\tadded_reads\\tnew_est_reads\\tfraction_total_reads\\nNo microbial reads\\t0\\t\$BRACKEN_LEVEL\\t0\\t0\\t0\\t0.00" > bracken_raw.tsv

    # 3. Sort by relative abundance (column 7) and keep header
    head -n 1 bracken_raw.tsv > ${SampleName}_bracken.tsv
    tail -n +2 bracken_raw.tsv | sort -t\$'\t' -k7,7nr >> ${SampleName}_bracken.tsv
    """
}