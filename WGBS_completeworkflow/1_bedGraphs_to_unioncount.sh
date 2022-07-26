#!/bin/bash
## Author: Dario Galanti May 2022
## Aim: Combine bedGraph files in a unioncount file, storing coverage info, apply coverage filtering, maxNA filtering and extract post-filtering average and weighted methylation
## NB: There is the option to scale all high coverage positions to a maximum treshold. Explanation is below before the command
## Documentation (EpiDiverse WGBS pipeline): https://github.com/EpiDiverse/wgbs
## Input: coverage-filtered bedgraph files with dir structure: wgbs_filt_output/${context}/${sample}.bedGraph (produced by 1_filter_bedGraphs.sh)
## Output: multisample simil-unionbed file with methylated and unmethylated counts (possibly with NAs depending on filtNAs parameter)
## Run: bash 1_bedGraphs_to_unioncount_flt.sh filtered_bedGraphs_Dir context
## Run: sbatch --partition test --cpus-per-task 4 --mem 40G --time 24:00:00 --wrap "bash 1_bedGraphs_to_unioncount.sh bedGraphs_v5_3cov CpG"
## Dependencies: bedtools (https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)
## Resources: 80G and 48h for CHH

## NB: The previous script (1_filter_bedGraphs.sh) produces QC files that should be used to check samples coverage.
## NB: Outlier samples with very low coverage should be removed before running this script by manually deleting the corresponding bedgraph files.
## I REMOVED "TA_SE_09_10_F1_HC0_M1_1" MANUALLY

## FILE DIRECTORIES AND NA FILTERING:
inDir=$1
cont=$2							#NB: This needs to be set in format (CpG, CHG, CHH)
min_cov=3						# Positions with lower coverage will be converted to "NA"
filtNAs=0.25					# Position with more NAs will be discarded
#max_cov=200					#OPTIONAL: Positions with higher coverage will be scaled to this value
outDir=unioncounts_v5_${min_cov}cov_${filtNAs}NAs
fout=$outDir/${cont}_unioncount_v5_${filtNAs}NAs.bed
#fout2=$outDir/${cont}_unioncount_v5_3to200cov_${filtNAs}NAs.bed
av_meth_fout=$outDir/${cont}_av_meth_${filtNAs}NAs.txt

## PREPARE OUT DIR
mkdir -p $outDir/
mkdir -p $outDir/input_bedgraphs/${cont}

## 	FORMAT AND FILTER BEDGRAPHS (Remove header and combine coverage columns to meth/total counts)
for f in ${inDir}/${cont}/*.bedGraph;
do
	temp_bedGr=$outDir/input_bedgraphs/${cont}/$(basename $f)
	awk min_cov=$min_cov 'NR>1{OFS="\t";tot=($5+$6);if(tot>min_cov){print $1,$2,$3,$5"/"tot}}' $f > $temp_bedGr
done

## MAKE UNIONBED (containing meth/total counts for each position) AND FILTER NAs
spl_fins=$(ls ${outDir}/input_bedgraphs/${cont}/*.bedGraph | tr "\n" " ")
spls=$(ls ${outDir}/input_bedgraphs/${cont}/*.bedGraph | rev | cut -d"/" -f1 | cut -d"." -f2- | rev | tr "\n" " ")
bedtools unionbedg -i ${spl_fins} -filler NA -header -names ${spls} | awk -v filtNAs="$filtNAs" '{OFS="\t";if(NR==1){print}else{NA=0;for(i=4;i<=NF;i++){if($i=="NA"){NA++}};if((NA/(NF-3))<=filtNAs){print}}}' > $fout

## REMOVE FORMATTED BEDGRAPHS
rm -r $outDir/input_bedgraphs/${cont}

## OPTIONAL: SCALE HIGH COVERAGE POSITIONS TO MAX THRESHOLD
## When applying weighted methylation, positions with very high coverage will have a much higher impact, and this is influenced by structural variation
## We limit this effect by twisting all positions with coverage higher than a treshold, to that treshold.
## Eg. Threshold=300. 100/500 -> (300*100/500)/300 -> 60/300
#awk -v MC=$max_cov '{if(NR==1){print}else{for(i=1;i<=NF;i++){i==NF?ORS="\n":ORS="\t";if(i<4 || $i=="NA"){print $i}else{split($i,a,"/");if(a[2]>MC){m=int((a[1]*MC/a[2])+0.499);print m"/"MC}else{print $i}}}}}' $fout > $fout2


## CALCULATE AVG AND WEIGHTED METHYLATION AFTER NA FILTERING
# We output also met and total counts so that weighted methylation can be calculated also accross contexts
head -1 $fout | cut -f4- | tr "\t" "\n" > $outDir/${cont}_tmp1.txt										#Print sample names in a column
#tail -n+2 $fout | awk '{for(i=4;i<=NF;i++){if($i!="NA"){meth[i]+=$i;c[i]++}}}END{for(i=4;i<=NF;i++){printf "%.2f\n",(meth[i]/c[i])}}' > $outDir/${cont}_tmp2.txt
awk 'NR>1{for(i=4;i<=NF;i++){if($i!="NA"){split($i,a,"/");metsum[i]+=a[1];totsum[i]+=a[2];AVGsum[i]+=(a[1]/a[2]);c[i]++}}} END{for(i=4;i<=NF;i++){printf "%.2f\t%s\t%s\t%.2f\n",(100*AVGsum[i]/c[i]),metsum[i],totsum[i],(100*metsum[i]/totsum[i])}}' $fout > $outDir/${cont}_tmp2.txt
echo -e Sample"\t"Avg_m${cont}_${filtNAs}NAs"\t"Met_${cont}counts"\t"Tot_${cont}counts"\t"Weight_m${cont}_${filtNAs}NAs > $av_meth_fout
paste $outDir/${cont}_tmp1.txt $outDir/${cont}_tmp2.txt >> $av_meth_fout
rm $outDir/${cont}_tmp*.txt

