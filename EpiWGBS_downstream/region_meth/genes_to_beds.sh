#!/bin/bash

### Aim: Using a gff annotation file, extract gene regions to 2 bed files (Fw and Rv genes)
### NB: This script might need to be adapted for other annotations using different separators
### Author: Dario Galanti Nov 2020
### Input: gene annotation file (/scr/epi/genomes/thlaspi_arvense/archive/thlaspi_clr/annotation/thlaspi.filtered.gff)
### Output: 2 bed files (Fw and Rv genes) with no headers. Columns are "chr" "start" "end" "strand" "geneID"
### Run: bash genes_to_beds.sh thlaspi.gff outdir
### Run: bash genes_to_beds.sh /scr/epi/genomes/thlaspi_arvense/archive/thlaspi_clr/annotation/thlaspi.filtered.gff /scr/episan/RP07/region_meth/gene_meth/gene_meth

### IMPORTANT TO KEEP IN MIND !!!
## Start position in gff is 1 based, in bed is 0 based!!!
## There are overlapping genes on the two different strands, but usually not on the same strand!

## DEFINE INPUT FILE
reg_gff=$1
## DEFINE OUTPUT DIRECTORY
outDir=$2

for i in "Fw" "Rv";
do
	if [ $i == "Fw" ];then strand="+"; else strand="-";fi
	
	# DEFINE OUTPUT FILE NAME
	reg_bed=${outDir}/${i}_genes.bed
##  Split gff by strand, subset genes and convert to bed
	awk -v str="$strand" -F"\t|;" 'OFS="\t"{if ($3=="gene" && $7==str) {split($9,id,"=");print $1,($4-1),$5,$7,id[2]}}' $reg_gff | sort -k1,1V -k2,2n > $reg_bed	# v3 and v5 annotation
	overlap_bp=$(bedtools merge -i $reg_bed -d -1 -c 3 -o count | awk '$4>1' | wc -l)
	if [ $overlap_bp == 0 ]; then echo No overlapping genes on $i strand, all good
	else echo WARNING: There are overlapping genes in the $i strand. The total overlap is $overlap_bp bp
	fi
	
done
