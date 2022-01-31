#!/bin/bash
## Author: Dario Galanti Feb 2021
## Aim: Prepare unionbed file to run PCA. Filter for NAs and unvariable positions
## Input: Unionbed file
## Run: sh 3_filter_unionbed.sh /scr/episan/RP07/bam_aligned/WGBS/unionbed_v5_3cov_25NAs/CHH_unionbed_v5_25NAs.bed /scr/episan/RP07/bam_aligned/WGBS/unionbed_v5_3cov_5NAs_5MEC_ED/CHH_unionbed_v5_5NAs.bed /scr/episan/RP07/bam_aligned/WGBS/unionbed_v5_3cov_5NAs_5MEC_ED/CHH_unionbed_v5_5NAs_2MEC_5ED.bed
## Run: sbatch --partition test --cpus-per-task 4 --mem 20G --time 20:00:00 --wrap "sh 3_filter_unionbed.sh"

## DEFINE INPUT AND OUTPUT
fin=$1
fout1=$2
fout2=$3

cont=$(basename $unionbed | grep -o 'C[pH][HG]')	# Extract context

## DEFINE FILTERS
NAfilt=0.02
# IMPORTANT: Define MEF and ED to further filter DMRs
#MEF=0.05	# Define Minor Epiallele Frequency (proportion of samples which need to have differential methylation from the others)
MEC=2	# Define Minor Epiallele Count (number of samples which need to have differential methylation from the others)
# Define Epiallele Difference (ED): Minimum methylation difference to define different epialleles.
#if [ ${cont} == "CpG" ];then ED=20;else ED=15;fi
ED=5

## 1) FILTER NAs
awk -v NAfilt="$NAfilt" 'OFS="\t"{if(NR==1){print}else{NA=0;for(i=4;i<=NF;i++){if($i=="NA"){NA++}};if((NA/(NF-3))<=NAfilt){print}}}' $fin > $fout1

## 2) FILTER NON-VARIABLE POSITIONS
## Using MEF
#awk -v MEF=$MEF -v ED=$ED 'OFS="\t"{if(NR==1){print}else{for(i=4;i<=NF;i++){if($i!="NA"){list[i]=$i}};asort(list);
#spls=length(list); MEC=int(MEF*spls+0.999); MECdiff=(list[(spls-MEC+1)]-list[MEC]);
#if(MECdiff>ED){print}; delete list}}' $fout1 > $fout2

## Using MEC
awk -v MEC=$MEC -v ED=$ED 'OFS="\t"{if(NR==1){print}else{for(i=4;i<=NF;i++){if($i!="NA"){list[i]=$i}};asort(list);
spls=length(list); MECdiff=(list[(spls-MEC+1)]-list[MEC]);
if(MECdiff>ED){print}; delete list}}' $fout1 > $fout2

rm $fout1

