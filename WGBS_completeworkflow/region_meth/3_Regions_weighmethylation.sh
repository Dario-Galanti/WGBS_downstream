#!/bin/bash

### Author: Dario Galanti Feb 2021
### Aim: Extract methylation statistics of individual regions (genes, CDSs, TEs or promoters) from a unioncount file.
### Aim: depending on the python script used to iterate through regions it can produce three different stats for each region:
### 	1): unioncount_reg_avgmet.py	-> calculates simple average methylation (mean metylation in Schultz et al. 2012).
### 	2): unioncount_reg_weighmet.py	-> calculates weighted methylation (see Schultz et al. 2012).
### 	3): unioncount_reg_fracMetCs.py -> calculates the fraction of methylated cytosines (see Schultz et al. 2012).
### Explanation: genes/promoters on different strands often overlap and unioncount_reg_...py need discrete non-overlapping regions, so we run it two times, once for each strand.
### Input 1): position-sorted ("sort -k1,1V -k2,2n") unioncount file with meth/total read counts, already intersected with regions
### Input 2): Two bed files for Forward and Reverse regions (genes, TEs or promoters). Columns 4 and 5 should contain "strand" and "geneID" info (extra columns will be ignored)
### Run: bash 3_Regions_methylation.sh ${context}_${feature}_unionbed.bed regions_${feature}.bed 
### Run: sbatch --partition test --cpus-per-task 4 --mem 20G --time 18:00:00 --wrap "bash 3_Regions_weighmethylation.sh /scr/episan/RP07/region_meth/Features_meth_v5/feature_unionbeds/CpG_CDSs_unionbed_v5_25NAs.bed"
### Dependencies: either unioncount_reg_avgmet.py / unioncount_reg_weighmet.py / unioncount_reg_fracMetCs.py

### Steps:
### 3) Create index file: a txt file with one column, containing all possible scaffolds (sorted)
### 4) Intersect regions.bed and unionbed and discard regions with less than x=5 (can be defined) Cs covered (unreliable data)
### 5) Extract methylation of individual regions (genes/TEs/promoters) -> Adam's average_over_bed.py
### 6) Add strand, ID and coverage info to the output file ($fout)
### 7) merge Fw and Rv regions methylation results
### 8) Record average methylation distribution of analysed regions
### 9) Find differentially methylated regions (genes or promoters). This can be used for further analysis such as enrichment analysis
### 10)Extract average methylation of regions with average methylation higher than x%. Genes and promoters with really low average methylation are a source of noise (Kawakatsu et al. 2016 removed genes with less than 5% meth)

### IMPORTANT TO KEEP IN MIND !!!
## There are overlapping genes/promoters on the two different strands but not on the same strand!

## DEFINE CODE, INPUT AND PARAMETERS
unionbed=$1													# NB Unionbed filename has to be ${context}_${feature}_unionbed.bed !!!!!!!!!!
wDir=/scr/episan/RP07/region_meth/Features_meth_v5
Fw_regions=${wDir}/feature_beds/Ta_v5_Fw_promoters.bed
Rv_regions=${wDir}/feature_beds/Ta_v5_Rv_promoters.bed
avg_over_bed=${wDir}/unioncount_reg_avgmet.py				# DEFINE HERE: script to extract regions methylation stats from unioncount and regions files
cont=$(echo $unionbed | grep -o 'C[pH][HG]')				# Extract context
feature=$(basename $unionbed .bed | cut -d_ -f2)			# Extract feature from unionbed
x=5				# NB: Number of Cs that should be covered within each region feature in order to analyse it (NB: Context specific, so don't be too strict!!). This removes unreliable averages
minMet=5		# Minimum methylation % for regions to be considered for calculating the average. Genes/promoters with lower methylaiton will be excluded from the accessions average calculation

## DEFINE OUTPUT, CONTEXT AND TEMPORARY FILES
outdir=${wDir}/${feature}_meth
mkdir -p $outdir
index=${outdir}/${cont}_${feature}_index.txt							# Temporary file, will be created and deleted
regions_unionbed=${outdir}/union_m${cont}_indiv_${feature}.bed			# Output unionbed file with methylation of individual regions (Fw and Rv so some overlap)
reg_meth_distr=${outdir}/Distribution_m${cont}_${feature}.txt			# Distribution summary, will be created
feature_avmeth=${outdir}/${cont}_avg_${minMet}minmeth_${feature}.txt	# Will be created


## 3) Create index: a txt file with one column, containing all possible scaffolds (sorted)
## It is generated from the unionbed because after intersecting, it will include all region file scaffolds
tail -n+2 $unionbed | cut -f1 | uniq > $index

## Loop through Fw and Rv regions
for i in "Fw" "Rv";
do
	# DEFINE FILE NAMES
	if [ $i == "Fw" ];then reg_bed=$Fw_regions; else reg_bed=$Rv_regions;fi
	fin=${outdir}/covered_${cont}_${i}_${feature}.bed
	fout=${outdir}/m${cont}_${i}_${feature}.bed
	echo $fout report
	
## 4) Intersect region_bed and unionbed to discard uncovered regions (less than "x" Cs) and report coverage of covered ones for further filtering
	awk 'NR!=1{OFS="\t";print $1,$2,$3}' $unionbed | bedtools intersect -a $reg_bed -b stdin -c | awk -v x=$x '{if($NF>x){OFS="\t";print $1,$2,$3,$4,$5,$NF}}' | sort -k1,1V -k2,2n > $fin
	
## 5) Extract methylation of individual regions -> Adam's average_over_bed.py
	echo Running average_over_bed.py for $i strand ...
	python3 $avg_over_bed $fin $unionbed $index $fout

	# Add header to fin to have regions on same line as fout
	echo -e chrom"\t"start"\t"end"\t"Strand"\t"Info"\t"Covered_Cs > ${outdir}/tmp1_$(basename $fin)
	cat $fin >> ${outdir}/tmp1_$(basename $fin)
	mv ${outdir}/tmp1_$(basename $fin) $fin
	
## 6) Add strand, ID and coverage info to the methylation file ($fout)
	cut -f-6 $fin > ${outdir}/tmp2_$(basename $fin)
	cut -f4- $fout > ${outdir}/tmp3_$(basename $fout)
	paste ${outdir}/tmp2_$(basename $fin) ${outdir}/tmp3_$(basename $fout) > $fout
	rm ${outdir}/tmp*${cont}_${i}_${feature}.bed
done

## 7) merge Fw and Rv regions methylation results
head -1 ${outdir}/m${cont}_Fw_${feature}.bed > $regions_unionbed
awk 'FNR!=1' ${outdir}/m${cont}_*_${feature}.bed | sort -k1,1V -k2,2n >> $regions_unionbed

## Cleanup
rm ${outdir}/covered_${cont}_*_${feature}.bed	$index # rm fins and index

## 8) Record average methylation distribution of analysed regions
echo -e Meth"\t"Num_of_${feature} > $reg_meth_distr
for min in {0..98..2};
do
	echo -e -n ${min}-$(( $min + 2 ))%"\t" >> $reg_meth_distr
	awk -v min=$min 'NR!=1{c=0;s=0;max=(min+10);for(i=7;i<=NF;i++){if($i!="NA"){s+=$i;c++}};av=(s/c); if(av>=min && av<max){print}}' $regions_unionbed | wc -l >> $reg_meth_distr
done

## 9) OPTIONAL 1: Find differentially methylated features
MEF=0.05	# Define Minor Epiallele Frequency (proportion of samples which need to have differential methylation from the others)
# Define Epiallele Difference (minimum methylation difference to define different epialleles) NB: Use 15 or 20 for CpG and 10 or 15 for CHG and CHH
if [ ${cont} == "CpG" ];then ED=20;else ED=15;fi
awk -v MEF=$MEF -v ED=$ED 'OFS="\t"{if(NR==1){print}else{for(i=7;i<=NF;i++){if($i!="NA"){list[i]=$i}};asort(list);
spls=length(list); MEC=int(MEF*spls+0.999); MECdiff=(list[(spls-MEC+1)]-list[MEC]);
if(MECdiff>ED){print}; delete list}}' $regions_unionbed > ${outdir}/${cont}_${ED}ED_5MEF_${feature}.bed

## 10) OPTIONAL 2: Extract average methylation of regions with average methylation higher than x%
## NB: Not to count Cs in overlapping genes twice we go back to the unionbed files for calculating the averages
awk -v x=$minMet 'OFS="\t"{if(NR>1){s=0;c=0;for(i=7;i<=NF;i++){if($i!="NA"){s+=$i;c+=1}};av=(s/c);if(av>x){print $1,$2,$3}}}' $regions_unionbed | bedtools merge -i stdin > ${outdir}/${cont}_${minMet}minmeth_${feature}_merged.bed	# Filter regions and merge them
head -1 $unionbed | cut -f4- | tr "\t" "\n" > ${outdir}/${cont}_${feature}_tmp4.txt									#Print sample names in a column
## Calculate avg and weighted methylation in methylated regions
tail -n+2 $unionbed | bedtools intersect -a stdin -b ${outdir}/${cont}_${minMet}minmeth_${feature}_merged.bed -u | \
awk '{for(i=4;i<=NF;i++){if($i!="NA"){split($i,a,"/");metsum[i]+=a[1];totsum[i]+=a[2];AVGsum[i]+=(a[1]/a[2]);c[i]++}}} \
END{for(i=4;i<=NF;i++){printf "%.2f\t%.2f\n",(100*AVGsum[i]/c[i]),(100*metsum[i]/totsum[i])}}' > ${outdir}/${cont}_${feature}_tmp5.txt

echo -e Sample"\t"Avg_m${cont}_${feature}_${minMet}minmeth"\t"Weighted_m${cont}_${feature}_${minMet}minmeth > $feature_avmeth
paste ${outdir}/${cont}_${feature}_tmp4.txt ${outdir}/${cont}_${feature}_tmp5.txt >> $feature_avmeth
rm ${outdir}/${cont}_${feature}_tmp*.txt ${outdir}/${cont}_${minMet}minmeth_${feature}_merged.bed


