# region_meth
Scripts to extract methylation of genomic features (genes, promoters, TEs) from multisample unionbed files.


SCRIPTS DESCRIPTION: <br/>
[genes_to_beds.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/region_meth/genes_to_beds.sh)<br/>
Simple script to extract bed files for forward (Fw) and reverse (Rv) genes from an annotation.

[TEs_to_beds.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/region_meth/TEs_to_beds.sh)<br/>
Script to extract bed files for forward (Fw) and reverse (Rv) TEs from the Thlaspi TE annotation. Overlapping TEs on the same strand are merged and superfamilies composing each merged region are recorded.

[get_promoters.py](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/region_meth/get_promoters.py)<br/>
Python script using the output of [genes_to_beds.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/region_meth/genes_to_beds.sh) to extract promoter regions. Promoter size is user defined (2kb default) or shorter in case the previous gene on the same strand is closer than the defined distance.

[Regions_methylation.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/region_meth/Regions_methylation.sh)<br/>
Wrapper script extracting methylation of two sets of regions stored in two bed files (it can use the output from any of [genes_to_beds.sh](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/region_meth/genes_to_beds.sh), [get_promoters.py](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/region_meth/get_promoters.py) or [TEs_to_beds.sh](https://github.com/Dario-GalantiWGBS_downstream/blob/main/region_meth/TEs_to_beds.sh)) by running [average_over_bed.py](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/average_over_bed.py). Regions within each of the two input files have to be non-overlapping, but regions in the two different files can overlap (eg. genes on different strands). It extracts i) Average methylation of all regions (eg. all genes), ii) Methylation of individual regions (average meth of each individual gene) and iii) Average methylation of all regions with methylation higher then a user defined cutoff (eg. meth of all genes with av. meth > 5% across all accessions).

[unionbed_reg_AVG_allcont.py](https://github.com/Dario-Galanti/WGBS_downstream/blob/main/region_meth/unionbed_reg_AVG_allcont.py)<br/>
Python script to plot the average (and optionally st.dev) distributions of regions in all contexts at the same time.

