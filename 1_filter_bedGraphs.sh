#!/bin/bash
## Author: Dario Galanti Feb 2021
## Aim: Filter single-sample bedGraph files (EpiDiverse WGBS output) by coverage and output some stats (coverage and n° positions covered).
## Aim: Positions with cov < cov_flt are discarded.
## Input: EpiDiverse wgbs output Directory with dir structure: wgbs_output/${context}/${sample}.bedGraph
## Output: Coverage-filtered bedgraph files with same dir structure as input (wgbs_filt_output/${context}/${sample}.bedGraph)
## Run: bash 1_filter_bedGraphs.sh input_Dir context
## Run: sbatch --partition test --cpus-per-task 2 --mem 40G --time 20:00:00 --wrap "bash 1_filter_bedGraphs.sh wgbs_v3 CpG"

## NB: The QC file produced by this script should be used to check samples coverage.
## NB: Outlier samples with very low coverage should be removed by manually deleting the corresponding bedgraph files in the output of this script.
## I REMOVED "TA_SE_09_10_F1_HC0_M1_1" MANUALLY

## FILES DIRECTORIES AND COV FILTERING:
cov_flt=3			# Position coverage: If p_cov > this parameter, position is kept for further analysis
inDir=$1
cont=$2						#NB: This needs to be set in format (CpG, CHG, CHH)
outDir=bedGraphs_v3_${cov_flt}cov
av_meth_fout=$outDir/${cont}_av_meth.txt

## PREPARE OUT DIR AND QC FILE
mkdir -p $outDir/${cont}
#echo -e Sample"\t"Av_cov_${cont}"\t"Pos_${cov_flt}cov_filt"\t"Av_${cont}_meth > $av_meth_fout				#Stats QC file header (1)
echo -e Sample"\t"Pos_${cov_flt}cov_filt"\t"Av_${cont}_meth > $av_meth_fout							#Stats QC file header (2) Slightly quicker!

for f in $inDir/${cont}/*.bedGraph;
do
	spl=$(basename $f .bedGraph)
	fout=$outDir/$cont/$(basename $f)
	awk -v filt="$cov_flt" '{if(NR==1){print}else{cov=($5+$6);if(cov>filt){OFS="\t"; print}}}' $f > $fout
	#tail -n+2 $f | awk -v filt="$cov_flt" -v spl="$spl" '{c+=($5+$6);if(cov>filt){n++;m+=$4}} END {OFS="\t"; print spl,(c/NR),n,(m/n)}' >> $av_meth_fout	#Stats QC file (1)
	tail -n+2 $fout | awk -v spl="$spl" '{m+=$4} END {OFS="\t"; print spl,NR,(m/NR)}' >> $av_meth_fout		#Stats QC file (2) Slightly quicker!
done



