# WGBS_downstream
Collection of scripts for downstream analysis of methylation bedgraph files produced by the EpiDiverse [WGBS pipeline](https://github.com/EpiDiverse/wgbs)

The EpiDiverse [WGBS pipeline](https://github.com/EpiDiverse/wgbs) is a great tool performing mapping and methylation calling of Bisulfite Sequencing datasets from non-model species. This repository contains scripts for downstream analysis of the [WGBS pipeline](https://github.com/EpiDiverse/wgbs) output, including coverage filtering, samples merging, NA filtering and region-specific average methylation calculation.
<br/> 
<br/> 
The repo is composed of two main folders:

[EpiWGBS_downstream](https://github.com/Dario-Galanti/WGBS_downstream/tree/main/EpiWGBS_downstream)<br/>
Contains scripts for initial downstream analysis of the single-sample bedgraph files produced by the [WGBS pipeline](https://github.com/EpiDiverse/wgbs), including coverage filtering, merging samples into uninbed files, NA filtering and filtering of position based on variability between samples.

[region_meth](https://github.com/Dario-Galanti/WGBS_downstream/tree/main/region_meth)<br/>
Contains scripts to extract methylation of genomic features (genes, promoters, TEs) from multisample unionbed files.

