# WGBS_completeworkflow

Scripts for downstream analysis of the single-sample bedgraph files produced by the [WGBS pipeline](https://github.com/EpiDiverse/wgbs), including coverage filtering, merging samples into unincount files and NA filtering.<br/>
In addition, it is possible to:<br/>
1) calculate both mean and weighted methylation of regions according to [Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0), optionally correcting for non-conversion rates.<br/>
2) classify genes methylation status using a binomial model adapted from [Takuno and Gaut 2013](https://www.pnas.org/doi/abs/10.1073/pnas.1215380110) and [Niedethuth et al. 2016](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-1059-0).<br/>

The workflow is based on "unioncount" files (positions as rows and samples as columns) with the methylated/total read counts for each sample and position, which allow to store both methylation and coverage info.<br/>
Eg.<br/>
chrom start end Sample1  Sample2  Sample3 ...<br/>
Chr1  2234 2235 3/22 4/28 10/11
<br/>

[1_bedGraphs_to_unioncount.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/1_bedGraphs_to_unioncount.sh)<br/>
Combine single-sample bedGraph files from the [WGBS pipeline](https://github.com/EpiDiverse/wgbs) in a unioncount file, apply coverage filtering, maxNA filtering and extract post-filtering average and weighted methylation according to [Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0).
<br/>

[2a_unioncount_binom_byscaff.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/2a_unioncount_binom_byscaff.sh)<br/>
Wrapper script to run [unioncount_binom.R](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/unioncount_binom.R) parallelizing by scaffold and perform non-conversion rates correction ([Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0)). The parallelization uses a queuing system to parallelize by scaffold and save time and memory.
<br/>

[unioncount_binom.R](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/unioncount_binom.R)<br/>
Correct for non-conversion rates as described in [Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0). Perform a binomial test on each C to test if methylation is higher than background frequency (non-conversion rate) and convert methylated counts to zero if binomial test is not significant (P > 0.01).
<br/>

[2b_comb_unioncount.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/2b_comb_unioncount.sh)<br/>
Combine individual scaffold results from [2a_unioncount_binom_byscaff.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/2a_unioncount_binom_byscaff.sh) and extract average and weighted methylation [Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0).
<br/>

[region_meth](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/WGBS_completeworkflow/region_meth)<br/>
Collection of scripts to extract mean and weighted methylation or fraction of methylated cytosines (all methylation stats described in [Schultz et al 2012](https://www.cell.com/trends/genetics/fulltext/S0168-9525(12)00171-0)) of either:<br/>
1) whole genomic features (eg. all genes, CDS ...)<br/>
2) individual regions (eg. individual genes, CDS, TEs, DMRs...).
<br/>Gene Body Methylation classification is also included
<br/>

