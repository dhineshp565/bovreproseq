#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process make_metagenomics_report {
        publishDir "${params.out_dir}/",mode:"copy"

        label "low"
        input:
        path (csv)
        path(krona_reports_raw)
        path(rmdfile)
        path (bracken)
        path(blast)
        output:
        path("*.html")
        script:
        """

        cp ${csv} samples.csv
        cp ${krona_reports_raw} rawreads.html
        cp ${rmdfile} report.Rmd


        Rscript -e 'rmarkdown::render(input="report.Rmd",params=list(csv="samples.csv",krona="rawreads.html"),output_file=paste0("ont_metageomics_results_report_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"), ".html"))'
        """

}
