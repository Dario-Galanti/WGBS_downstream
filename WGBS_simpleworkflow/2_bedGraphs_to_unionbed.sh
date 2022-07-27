#!/bin/bash
## Author: Dario Galanti Feb 2021
## Aim: Combine bedGraph files in a unionbed applying maxNA filtering and extracting post-filtering average methylation
## Documentation (EpiDiverse WGBS pipeline): https://github.com/EpiDiverse/wgbs
## Input: coverage-filtered bedgraph files with dir structure: wgbs_filt_output/${context}/${sample}.bedGraph (produced by 1_filter_bedGraphs.sh)
## Output: multisample unionbed file with methylation calls (possibly with NAs depending on filtNAs parameter)
## Run: bash 2_bedGraphs_to_unionbed.sh filtered_bedGraphs_Dir context
## Run: sbatch --partition test --cpus-per-task 2 --mem 40G --time 20:00:00 --wrap "bash 2_bedGraphs_to_unionbed.sh bedGraphs_v3_3cov CpG"
## Dependencies: bedtools (https://bedtools.readthedocs.io/en/latest/content/bedtools-suite.html)

## NB: The previous script (1_filter_bedGraphs.sh) produces QC files that should be used to check samples coverage.
## NB: Outlier samples with very low coverage should be removed before running this script by manually deleting the corresponding filtered bedgraph files.
## I REMOVED "TA_SE_09_10_F1_HC0_M1_1" MANUALLY

## FILE DIRECTORIES AND NA FILTERING:
inDir=$1
cont=$2						#NB: This needs to be set in format (CpG, CHG, CHH)
filtNAs=0.25
cov_flt=$(echo $inDir | rev | cut -d"_" -f1 | rev)
outDir=unionbed_v3_${cov_flt}_${filtNAs}NAs
tmp_fout=$outDir/${cont}_unionbed_v3.bed
fout=$outDir/${cont}_unionbed_v3_${filtNAs}NAs.bed
av_meth_fout=$outDir/${cont}_av_meth_${filtNAs}NAs.txt

## PREPARE OUT DIR
mkdir -p $outDir/
mkdir -p $outDir/input_bedgraphs/${cont}

## 	FORMAT BEDGRAPHS (Remove header and coverage columns)
for f in ${inDir}/${cont}/*.bedGraph;
do
	temp_bedGr=$outDir/input_bedgraphs/${cont}/$(basename $f)
	#tail -n+2 $f | cut -f-4 > $temp_bedGr					# This bedGraph column is rounded down > underestimation of methylation levels
	tail -n+2 $f | awk '{met=100*$5/($5+$6);printf "%s\t%s\t%s\t%.2f\n", $1,$2,$3,met}' > $temp_bedGr		# This is more precise
done

## MAKE UNIONBED
spl_fins=$(ls ${outDir}/input_bedgraphs/${cont}/*.bedGraph | tr "\n" " ")
spls=$(ls ${outDir}/input_bedgraphs/${cont}/*.bedGraph | rev | cut -d"/" -f1 | cut -d"." -f2- | rev | tr "\n" " ")
bedtools unionbedg -i ${spl_fins} -filler NA -header -names ${spls} > $tmp_fout

## NA FILTERING
head -1 $tmp_fout > $fout
tail -n+2 $tmp_fout | awk -v filtNAs="$filtNAs" '{NA=0;for(i=4;i<=NF;i++){if($i=="NA"){NA++}};if((NA/(NF-3))<=filtNAs){print}}' >> $fout

## REMOVE FORMATTED BEDGRAPHS AND UNFILTERED UNIONBEDG
rm -r $outDir/input_bedgraphs/${cont}
#rm $tmp_fout

## CALCULATE AV METH AFTER NA FILTERING
head -1 $fout | cut -f4- | tr "\t" "\n" > $outDir/${cont}_tmp1.txt																			#Print sample names in a column
tail -n+2 $fout | awk '{for(i=4;i<=NF;i++){if($i!="NA"){meth[i]+=$i;c[i]++}}}END{for(i=4;i<=NF;i++){printf "%.2f\n",(meth[i]/c[i])}}' > $outDir/${cont}_tmp2.txt
echo -e Sample"\t" ${cont}_avmeth_${filtNAs}NAs > $av_meth_fout
paste $outDir/${cont}_tmp1.txt $outDir/${cont}_tmp2.txt >> $av_meth_fout
rm $outDir/${cont}_tmp*.txt

