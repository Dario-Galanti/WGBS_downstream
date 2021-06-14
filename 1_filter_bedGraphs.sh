#!/bin/bash
## Author: Dario Galanti Feb 2021
## Aim: Filter bedGraph files by coverage
## Input: wgbs outDir with dir tree: wgbs_output/${context}/${sample}.bedGraph
## Output: multisample unionbed file harbouring methylation info for all cytosines and all samples
## Run: bash 1_filter_bedGraphs.sh input_Dir context
## Run: sbatch --partition test --cpus-per-task 2 --mem 40G --time 20:00:00 --wrap "bash 1_filter_bedGraphs.sh wgbs_v3 CpG"

## TO IMPLEMENT: Add an average sample coverage filtering to remove samples with very low coverage!!
## IMPORTANT: For now sample filtering has to be done manually looking at the QC files
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



