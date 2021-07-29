# EpiWGBS_downstream
Collection of scripts for downstream analysis of methylation bedgraph files produced by the EpiDiverse [WGBS pipeline](https://github.com/EpiDiverse/wgbs)

The EpiDiverse [WGBS pipeline](https://github.com/EpiDiverse/wgbs) is a great tool performing mapping and methylation calling of Bisulfite Sequencing datasets from non-model species. This repository contains scripts for downstream analysis of the single-sample bedgraph files produced by the pipeline, including coverage filtering, samples merging, NA filtering and region specific methylation calculation.
<br/> 
<br/> 

SCRIPTS DESCRIPTION: <br/>
[1_filter_bedGraphs.sh](https://github.com/Dario-Galanti/EpiWGBS_downstream/blob/main/1_filter_bedGraphs.sh)<br/>
Coverage filtering of single-sample bedGraph files and sample specific QC. Positions with cov < cov_flt are discarded.

[2_bedGraphs_to_unionbed.sh](https://github.com/Dario-Galanti/EpiWGBS_downstream/blob/main/2_bedGraphs_to_unionbed.sh)<br/>
Combining coverage-filtered bedGraph files in a unionbed file, applying maxNA filtering and extracting post-filtering average methylation.

[average_over_bed.py](https://github.com/Dario-Galanti/EpiWGBS_downstream/blob/main/average_over_bed.py)<br/> 
Python script to extract average methylation of distinct non-overlapping regions (stored in a bed file), from a multisample position unionbed file. It outputs a unionbed file with regions as rows instead of positions.

