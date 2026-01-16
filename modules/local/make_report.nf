#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process make_report {
	publishDir "${params.out_dir}/",mode:"copy"
	
	label "low"
	input:
	path (csv)
	path(krona_reports_raw)
	path(mappedreads)
	path(abricate)
	path(results)
	path(rmdfile)
	path(mlst)
	path(cons)
    path (bracken)
	path(blast)
	output:
	path("*.html")
	script:
	"""
	
	cp ${csv} samples.csv
	cp ${krona_reports_raw} rawreads.html
	cp ${rmdfile} report.Rmd
	

	Rscript -e 'rmarkdown::render(input="report.Rmd",params=list(csv="samples.csv",krona="rawreads.html"),output_file=paste0("Bovreproseq_results_report_", Sys.Date(), "_", format(Sys.time(), "%H-%M-%S"), ".html"))'
	"""

}
