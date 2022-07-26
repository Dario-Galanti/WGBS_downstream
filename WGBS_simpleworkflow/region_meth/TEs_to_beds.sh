#!/bin/bash

### Aim: Using a TE gff annotation file, merge strand-specific overlapping TEs and extract resulting regions to 2 bed files (Fw and Rv TEs)
### NB: This script is specific for the TE Thlaspi annotation, it might need to be adapted for other annotations
### Author: Dario Galanti Nov 2020
### Input: TE annotation file (/scr/epi/genomes/thlaspi_arvense/archive/thlaspi_clr/annotation/tes/thlaspi.TE.gff)
### Output: 2 bed files (Fw and Rv TEs) with no headers. Columns are "chr" "start" "end" "strand" "TEfamilies"
### Run: bash TEs_to_beds.sh thlaspi.TE.gff
### Run: bash TEs_to_beds.sh /scr/epi/genomes/thlaspi_arvense/archive/thlaspi_clr/annotation/tes/thlaspi.TE.gff

### PRINCIPLE
## There are overlapping and nested TEs even on the same strand.
## Here we split the two strands to diminish the overlapping problem and merge all strand-specific overlapping TEs.
## We record the superfamilies present in each merged region

### IMPORTANT TO KEEP IN MIND !!!
## Start position in gff is 1 based, in bed is 0 based!!!


## DEFINE INPUT FILES
TE_gff=$1
## DEFINE OUTPUT DIRECTORY
outDir=/scr/episan/RP07/region_meth/TE_meth/TE_meth

for i in "Fw" "Rv";
do
	if [ $i == "Fw" ];then str="+"; else str="-";fi
	
	# DEFINE OUTPUT FILE NAME
	reg_bed=${outDir}/${i}_mergedTEs.bed
##  Split gff by strand, convert to bed and merge TEs
	awk -v s="$str" '{if($7 == s){OFS="\t";fam=substr($9,8,3);print $1,($4-1),$5,$7,fam}}' $TE_gff | bedtools merge -i stdin -d -1 -c 4,5 -o distinct,distinct | sort -k1,1V -k2,2n > $reg_bed

done
