# Bovreproseq
Bioinformatics pipeline for targeted syndromic next-generation sequencing panel for simultaneous detection of pathogens associated with bovine reproductive failure\
For more info please see https://doi.org/10.1128/jcm.01433-24 \
Outputs html files with results, consensus sequences, metagenomics reports (Kraken2/Bracken/Blast), Krona plots, IGV reports and MultiQC quality reports.

## Usage
```bash
nextflow run main.nf \
    --input <path_to_fastq_directory> \
    --out_dir <output_directory> \
    --kraken_db <path_to_kraken2_db> \
    --blastdb_path <path_to_blast_db> \
    --blastdb_name <blast_db_name>
```

## Parameters

### Required
* `--input`: Path to input directory containing FASTQ files (can be gzipped).
* `--out_dir`: Directory where results will be saved.
* `--kraken_db`: Path to the Kraken2 database directory.
* `--blastdb_path`: Path to the directory containing BLAST databases.
* `--blastdb_name`: Name of the BLAST database to use.

### Optional
* `--dehost`: Enable host read filtering (default: true).
* `--host_db`: Path to kraken2 database with Bos taurus reference genome FASTA (required if `--dehost` is true).
* `--qscore`: Minimum mean Q-score filter for reads (default: 10).
* `--trim_barcodes`: Enable barcode and adapter trimming using Porechop (default: false).
* `--read_count_threshold`: Threshold for minimum read count for consensus generation amplicons (default: 10).
* `--consensus_mode`: Mode for consensus generation with samtools('simple' or 'bayesian').

## Workflow Description
1. **QC & Preprocessing**: Merges FASTQ files, filters by Q-score (Nanoq), and optionally trims barcodes (Porechop). Generates QC metrics (NanoPlot).
2. **Dehosting (Optional)**: Uses kraken2 filter out Bos taurus reads
3. **Targeted Amplicon Analysis** : 
    - Aligns reads to the `Bovreproseq_reference.fasta`.
    - Generates consensus sequences using minimap2, samtools and medaka.
    - Identifies pathogens in the sample using Abricate.
    - Performs MLST typing to get subspecies ID of Campylobacter fetus.
    - Generates IGV HTML reports for visual inspection of alignments.
4. **Metagenomics**: 
    - Classifies reads using Kraken2 and estimates abundance with Bracken.
    - Visualizes taxonomy with Krona.
    - Validation of hits using BLAST against specified databases.
5. **Reporting**: 
    - Consolidates results into an interactive HTML report.
    - Generates a separate Metagenomics report.
    - MultiQC report for run statistics.

## Dependencies
* nextflow
* docker
* wsl2 (if on Windows)

## Software Used
* **[nanoq](https://github.com/esteinig/nanoq)**: Steinig, E., & Coin, L. (2022). Nanoq: ultra-fast quality control for nanopore reads. *Journal of Open Source Software*, 7(69), 2991. https://doi.org/10.21105/joss.02991
* **[NanoPlot](https://github.com/wdecoster/NanoPlot)**: De Coster, W., D'Hert, S., & Schultz, D. T. (2018). NanoPack: visualizing and processing long-read sequencing data. *Bioinformatics*, 34(12), 2066–2068. https://doi.org/10.1093/bioinformatics/bty149
* **[porechop](https://github.com/rrwick/Porechop)**: Wick, R. R. (2017). Porechop: adapter trimmer for Oxford Nanopore reads. https://github.com/rrwick/Porechop
* **[minimap2](https://github.com/lh3/minimap2)**: Li, H. (2018). Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*, 34(18), 3094–3100. https://doi.org/10.1093/bioinformatics/bty191
* **[kraken2](https://ccb.jhu.edu/software/kraken2/)**: Wood, D. E., Lu, J., & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257. https://doi.org/10.1186/s13059-019-1891-0
* **[bracken](https://ccb.jhu.edu/software/bracken/)**: Lu, J., Breitwieser, F. P., Thielen, P., & Salzberg, S. L. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104. https://doi.org/10.7717/peerj-cs.104
* **[KrakenTools](https://github.com/jenniferlu717/KrakenTools)**: Lu, J. (2020). KrakenTools. https://github.com/jenniferlu717/KrakenTools
* **[krona](https://github.com/marbl/Krona/wiki)**: Ondov, B. D., Bergman, N. H., & Phillippy, A. M. (2011). Interactive metagenomic visualization in a Web browser. *BMC Bioinformatics*, 12(1), 385. https://doi.org/10.1186/1471-2105-12-385
* **[megahit](https://github.com/voutcn/megahit)**: Li, D., Liu, C.-M., Luo, R., Sadakane, K., & Lam, T.-W. (2015). MEGAHIT: an ultra-fast single-node solution for large and complex metagenomics assembly via succinct de Bruijn graph. *Bioinformatics*, 31(10), 1674–1676. https://doi.org/10.1093/bioinformatics/btv033
* **[samtools](http://www.htslib.org/)**: Danecek, P., Bonfield, J. K., Liddle, J., Marshall, J., Ohan, V., Pollard, M. O., ... & Li, H. (2021). Twelve years of SAMtools and BCFtools. *GigaScience*, 10(2), giab008. https://doi.org/10.1093/gigascience/giab008
* **[bedtools](https://bedtools.readthedocs.io/en/latest/)**: Quinlan, A. R., & Hall, I. M. (2010). BEDTools: a flexible suite of utilities for comparing genomic features. *Bioinformatics*, 26(6), 841–842. https://doi.org/10.1093/bioinformatics/btq033
* **[abricate](https://github.com/tseemann/abricate)**: Seemann, T. (2016). Abricate: mass screening of contigs for antimicrobial resistance or virulence genes. https://github.com/tseemann/abricate
* **[medaka](https://github.com/nanoporetech/medaka)**: Oxford Nanopore Technologies. Medaka: Sequence correction provided by Oxford Nanopore Technologies. https://github.com/nanoporetech/medaka
* **[mlst](https://github.com/tseemann/mlst)**: Seemann, T. (2016). mlst: scan contig files against PubMLST typing schemes. https://github.com/tseemann/mlst (PubMLST: Jolley, K. A., & Maiden, M. C. J. (2010). BIGSdb: Scalable analysis of bacterial genome variation at the population level. *BMC Bioinformatics*, 11, 595. https://doi.org/10.1186/1471-2105-11-595)
* **[blast+](https://blast.ncbi.nlm.nih.gov/Blast.cgi)**: Camacho, C., Coulouris, G., Avagyan, V., Ma, N., Papadopoulos, J., Bealer, K., & Madden, T. L. (2009). BLAST+: architecture and applications. *BMC Bioinformatics*, 10, 421. https://doi.org/10.1186/1471-2105-10-421
* **[igv-reports](https://github.com/igvteam/igv-reports)**: Robinson, J. T., Thorvaldsdóttir, H., Winckler, W., Guttman, M., Lander, E. S., & Getz, G. (2011). Integrative genomics viewer. *Nature Biotechnology*, 29(1), 24–26. https://doi.org/10.1038/nbt.1754
* **[multiqc](https://multiqc.info/)**: Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. *Bioinformatics*, 32(19), 3047–3048. https://doi.org/10.1093/bioinformatics/btw354
* **[rmarkdown](https://rmarkdown.rstudio.com/)**: Allaire, J., Xie, Y., McPherson, J., Luraschi, J., Ushey, K., Atkins, A., Wickham, Cheng, J., Chang, W., & Iannone, R. (2021). rmarkdown: Dynamic Documents for R. R package version 2.7. https://rmarkdown.rstudio.com/
