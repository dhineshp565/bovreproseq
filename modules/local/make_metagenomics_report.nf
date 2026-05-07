#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process make_metagenomics_report {
        publishDir "${params.out_dir}/",mode:"copy"

        label "high"
        input:
        path (csv)
        path(krona_reports_raw)
        path(rmdfile)
        path (meta_summary)
        path(blast_hits)
        output:
        path("*.html")
        script:
        """

        cp ${csv} samples.csv
        cp ${krona_reports_raw} rawreads.html
        cp ${rmdfile} report.Rmd


        Rscript -e 'rmarkdown::render(input="report.Rmd",params=list(csv="samples.csv",krona="rawreads.html",meta_summary="${meta_summary}",blast_hits="${blast_hits}"),output_file=paste0("ont_metagenomics_results_report_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"), ".html"))'
        """

}
