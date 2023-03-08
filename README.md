# WGBS_downstream
Collection of scripts for downstream analysis of methylation bedgraph files produced by the EpiDiverse [WGBS pipeline](https://github.com/EpiDiverse/wgbs)

The EpiDiverse [WGBS pipeline](https://github.com/EpiDiverse/wgbs) is a great tool performing mapping and methylation calling of Bisulfite Sequencing datasets from non-model species. This repository contains scripts for downstream analysis of the [WGBS pipeline](https://github.com/EpiDiverse/wgbs) output, including coverage filtering, merging samples, NA filtering and region-specific average methylation calculation. Methylation of regions can be for both individual regions (eg. individual genes) or whole genomic features (eg. average methylation of all genes).
<br/> 
<br/> 
The repo is composed of two different workflows:

[WGBS_simpleworkflow](https://github.com/Dario-Galanti/WGBS_downstream/tree/main/WGBS_simpleworkflow)<br/>
Contains scripts for downstream analysis of the single-sample bedgraph files produced by the [WGBS pipeline](https://github.com/EpiDiverse/wgbs), including coverage filtering, merging samples into uninbed files, NA filtering and filtering of non-variable positions.
In addition, it is possible to calculate simple average methylation of regions (mean methylation in [Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0)).<br/>
The workhorse here are unionbed files (positions as rows and samples as columns) with the computed fraction of methylated/total reads for each sample and position.<br/>
Eg.<br/>
chrom start end Sample1  Sample2  Sample3 ...<br/>
Chr1  2234 2235 13.63 14.29 90.90
<br/> 
<br/> 

[WGBS_completeworkflow](https://github.com/Dario-Galanti/WGBS_downstream/tree/main/WGBS_completeworkflow)<br/>
Contains scripts for downstream analysis of the single-sample bedgraph files produced by the [WGBS pipeline](https://github.com/EpiDiverse/wgbs), including coverage filtering, merging samples into unincount files and NA filtering.
In addition, it is possible to:<br/>
1) calculate both mean and weighted methylation of regions according to [Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0), optionally correcting for non-conversion rates.<br/>
2) to classify genes methylation status using a binomial model adapted from [Takuno and Gaut 2013](https://www.pnas.org/doi/abs/10.1073/pnas.1215380110) and [Niedethuth et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1059-0).<br/>

The workhorse here are unioncount files (positions as rows and samples as columns) with the methylated/total read counts for each sample and position.<br/>
Eg.<br/>
chrom start end Sample1  Sample2  Sample3 ...<br/>
Chr1  2234 2235 3/22 4/28 10/11

<br/>

PUBLICATIONS: <br/>
[Genetic and environmental drivers of large-scale epigenetic variation in Thlaspi arvense](https://journals.plos.org/plosgenetics/article?id=10.1371/journal.pgen.1010452)
