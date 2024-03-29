# WGBS_simpleworkflow
Collection of scripts for downstream analysis of methylation bedgraph files produced by the EpiDiverse [WGBS pipeline](https://github.com/EpiDiverse/wgbs)

The EpiDiverse [WGBS pipeline](https://github.com/EpiDiverse/wgbs) is a great tool performing mapping and methylation calling of Bisulfite Sequencing datasets from non-model species. This folder contains scripts for initial downstream analysis of the single-sample bedgraph files produced by the [WGBS pipeline](https://github.com/EpiDiverse/wgbs), including coverage filtering, merging samples into uninbed files, NA filtering and filtering of positions based on variability between samples.
<br/> 
<br/> 

SCRIPTS DESCRIPTION: <br/>
[1_filter_bedGraphs.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_simpleworkflow/1_filter_bedGraphs.sh)<br/>
Coverage filtering of single-sample bedGraph files and sample specific QC. Positions with cov < cov_flt are discarded.

[2_bedGraphs_to_unionbed.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_simpleworkflow/2_bedGraphs_to_unionbed.sh)<br/>
Combining coverage-filtered bedGraph files in a unionbed file, applying maximum NAs filtering and extracting post-filtering average methylation.

[3_filter_unionbed.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_simpleworkflow/3_filter_unionbed.sh)<br/>
Filtering unionbed file for (1) maximum number of NAs (as in previous script); (2) Variability of positions: This method is similar to Minor Allele Frequency (MAF) filtering commonly used for genetic data and is based on defining a Minor Epiallele Frequency/Count (MEF/MEC), which is the proportion/number of samples which need to have differential methylation from all others and an Epiallele Difference (ED), which is the minimum methylation difference to define differential methylation between two samples. It can be used for example to exclude non variable positions before running clustering or PCA analysis.

[region_meth](https://github.com/Dario-Galanti/WGBS_downstream/tree/main/WGBS_simpleworkflow/region_meth)<br/>
Contains scripts to extract methylation of genomic features (genes, promoters, TEs) from multisample unionbed files.
