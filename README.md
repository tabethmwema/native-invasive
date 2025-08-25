I study how biological invasions shape genome–wide diversity by comparing native vs. invasive populations across multiple species. Using public reference genomes and whole-genome variant data, I quantify:

_ **Window-based nucleotide diversity (π)** in 10 kb windows with weighted means,

_ **Regional contrasts along chromosomes (telomere vs. centromere),

_ **Observed heterozygosity and per-species Δ (invasive − native),

_ **Pooled and per-species statistics (Wilcoxon tests, Cliff’s delta with CIs), with balanced down-sampling where needed.

The repo contains R scripts for summaries, effect sizes, and publication-ready figures that mirror the analyses reported in the manuscript.


### Data sources & extraction

_ **NCBI Datasets for reference genomes and assemblies (downloaded via datasets/wget).



### Alignment & variant calling toolchain

_ **minimap2 - whole-genome alignment/mapping to the reference

_ **samtools - BAM conversion, sorting, and indexing

_ **bcftools - mpileup + call for SNP/indel VCFs, basic filtering

_ **VCFtools - windowed π calculations (10 kb) and diversity summaries

### R packages
Post-processing, visualization, and statistics are done in R (dplyr, tidyr, ggplot2, gghalves, rstatix, effsize, readr, stringr).
