# region_meth

Collection of scripts to extract 1) mean methylation, 2) weighted methylation and 3) fraction of methylated cytosines (all methylation stats described in [Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0)) of either:<br/>
1) whole genomic features (eg. all genes, CDS ...)<br/>
2) individual regions (eg. individual genes, CDS, TEs, DMRs...).
<br/>

[1_Features_to_bed.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/region_meth/1_Features_to_bed.sh)<br/>
Extract individual genes, CDSs, introns, promoters and TEs to bed files, strating from .gff annotation files of Thlaspi [Nunn et al. 2022](https://onlinelibrary.wiley.com/doi/full/10.1111/pbi.13775).
<br/>

[2_Features_methylation_Epi.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/region_meth/1_Features_to_bed.sh)<br/>
Extract average methylation statistics of whole genomic features (eg. all genes, CDSs, promoters, TEs ...). First we intersect unioncount file with genomic features and then we extract the statistics. We parallelize each genomic feature to save memory and time.
<br/>

[3_Regions_methylation.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/region_meth/3_Regions_methylation.sh)<br/>
Extract average methylation statistics of individual regions (eg. all genes, CDSs, promoters ...). Starting from intersected unioncount files produced by [2_Features_methylation_Epi.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/region_meth/1_Features_to_bed.sh), and two bed files for Forward and Reverse regions, it can extract methylation statistics for each individual region.<br/>
NB: Whether 1) mean methylation, 2) weighted methylation or 3) fraction of methylated cytosines depends on the python script defined at the beginning of the file.<br/>
[unioncount_reg_avgmet.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/region_meth/unioncount_reg_avgmet.sh) calculates simple average methylation (mean methylation) <br/>
[unioncount_reg_fracMetCs.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/region_meth/unioncount_reg_fracMetCs.sh) calculates the fraction of methylated cytosines per region<br/>
[unioncount_reg_weighmet.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/region_meth/unioncount_reg_weighmet.sh) calculates weighted methylation <br/>

