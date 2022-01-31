#!/bin/bash

### Author: Dario Galanti Feb 2021
### Aim: Extract methylation of individual regions (genes or promoters) from unionbed file, using Adam Nunn's average_over_bed.py (https://github.com/Dario-Galanti/WGBS_downstream/blob/main/average_over_bed.py)
### Explanation: genes/promoters on different strands often overlap and average_over_bed.py needs discrete non-overlapping regions, so we run it two times, once for each strand.
### Input 1): position-sorted unionbed file (properly sorted, so "sort -k1,1V -k2,2n")
### Input 2): Two bed files for Forward and Reverse regions (genes, TEs or promoters). Columns 4 and 5 should contain "strand" and "geneID" info (extra columns might be present but will be ignored)
### Run: bash Regions_methylation.sh ${context}_unionbedg.bed regions.bed feature_name
### Run: sbatch --partition test --cpus-per-task 4 --mem 40G --time 18:00:00 --wrap "bash Regions_methylation.sh unionbed_v3_3cov_0.25NAs/CpG_unionbed_v3_0.25NAs.bed promoter"
### Dependencies: average_over_bed.py (https://github.com/Dario-Galanti/WGBS_downstream/blob/main/average_over_bed.py)

### Steps:
### 1) Subset unionbed to regions of interest (reduce RAM usage and running time)
### 2) Extract average methylation of all regions
### 3) Create index file: a txt file with one column, containing all possible scaffolds (sorted)
### 4) Intersect regions.bed and unionbed and discard promoters with less than x=5 (can be defined) Cs covered (unreliable data)
### 5) Extract methylation of individual regions (genes/promoters) -> Adam's average_over_bed.py
### 6) Add strand, ID and coverage info to the output file ($fout)
### 7) merge Fw and Rv regions methylation results
### 8) Record average methylation distribution of analysed regions
### 9) Find differentially methylated regions (genes or promoters). This can be used for further analysis such as enrichment analysis
### 10)Extract average methylation of regions with average methylation higher than x%. Genes and promoters with really low average methylation are a source of noise (Kawakatsu et al. 2016 removed genes with less than 5% meth)

### IMPORTANT TO KEEP IN MIND !!!
## Unionbed file has to be sorted correctly (Scaff_1 > Scaff_2 > ... > Scaff_10 > Scaff_11 ...), this can be obtained with "sort -k1,1V -k2,2n"
## There are overlapping genes/promoters on the two different strands but not on the same strand!


## DEFINE INPUT AND PARAMETERS 
wDir=/scr/episan/RP07/region_meth/promoter_meth/prom_meth
Fw_regions=${wDir}/Fw_promoter.bed
Rv_regions=${wDir}/Rv_promoter.bed
unionbed=$1
feature=$2
x=5				# NB: Number of Cs that should be covered within each region feature in order to analyse it (NB: Context specific, so don't be too restrictive!!). This removes unreliable averages
minMet=5		# Minimum methylation % for regions to be considered for calculating the average. Genes/promoters with lower methylaiton will be excluded from the accessions average calculation

## DEFINE OUTPUT, CONTEXT AND TEMPORARY FILES
average_over_bed=/scr/episan/RP07/region_meth/gene_meth/gene_methIMP/average_over_bed.py	# Adam Nunn script to extract regions methylation from unionbed
cont=$(echo $unionbed | grep -o 'C[pH][HG]')							# Extract context
subset_unionbed=${wDir}/${cont}_${feature}_unionbed_v3_0.25NAs.bed		# Will be created
index=${wDir}/${cont}_${feature}_index.txt								# Temporary file, will be created and deleted
merged_regions=${wDir}/Ta_${feature}_regions_merged.bed					# Will be created if not already existing
regions_AvgMeth=${wDir}/${cont}_AvgMeth_${feature}s.bed					# Output unionbed file with methylation of all regions (Fw and Rv so some will be overlapping)
regions_meth=${wDir}/${cont}_meth_individual_${feature}s.bed			# Output unionbed file with methylation of all regions (Fw and Rv so some will be overlapping)
reg_meth_distr=${wDir}/${cont}_meth_${feature}s_distribution.txt		# Will be created
feature_avmeth=${wDir}/${cont}_avg_${minMet}minmeth_${feature}s.txt	# Will be created

## 1) Subset unionbed to regions of interest (genes or promoters). This decreases RAM usage and running time.
head -1 $unionbed > $subset_unionbed
if [ ! -f $merged_regions ];then cat $Fw_regions $Rv_regions | sort -k1,1V -k2,2n | bedtools merge -i stdin -c 4,5 -o collapse,collapse > $merged_regions;fi	#NB: Cols 4 and 5 should have "strand" and "geneID" info
tail -n+2 $unionbed | bedtools intersect -a stdin -b $merged_regions >> $subset_unionbed

## 2) Extract average methylation of all regions (eg. average of Cs in all genes/promoters)
head -1 $subset_unionbed | cut -f4- | tr "\t" "\n" > ${wDir}/${cont}_tmp1.txt									#Print sample names in a column
awk '{for(i=4;i<=NF;i++){if($i!="NA"){meth[i]+=$i;c[i]++}}}END{for(i=4;i<=NF;i++){printf "%.2f\n",(meth[i]/c[i])}}' $subset_unionbed > ${wDir}/${cont}_tmp2.txt
echo -e Sample"\t" ${cont}_avg_meth_${feature}s > $regions_AvgMeth
paste ${wDir}/${cont}_tmp1.txt ${wDir}/${cont}_tmp2.txt >> $regions_AvgMeth
rm ${wDir}/${cont}_tmp*.txt

## 3) Create index: a txt file with one column, containing all possible scaffolds (sorted)
## It is generated from the unionbedg because after intersecting, it will include all region file scaffolds
tail -n+2 $subset_unionbed | cut -f1 | uniq > $index

## Loop through Fw and Rv genes
for i in "Fw" "Rv";
do
	# DEFINE FILE NAMES
	if [ $i == "Fw" ];then reg_bed=$Fw_regions; else reg_bed=$Rv_regions;fi
	fin=${wDir}/covered_${cont}_${i}_${feature}s.bed
	fout=${wDir}/meth_${cont}_${i}_${feature}s.bed
	echo $fout report
	
## 4) Intersect region_bed and unionbed to discard uncovered regions (less than "x" Cs) and report coverage of covered ones for further filtering
	awk 'NR!=1{OFS="\t";print $1,$2,$3}' $subset_unionbed | bedtools intersect -a $reg_bed -b stdin -c | awk -v x=$x '{if($NF>x){OFS="\t";print $1,$2,$3,$4,$5,$NF}}' | sort -k1,1V -k2,2n > $fin
	
## 5) Extract methylation of individual regions -> Adam's average_over_bed.py
	echo Running average_over_bed.py for $i strand ...
	python3 $average_over_bed $fin $subset_unionbed $index $fout

	# Add header to fin to have regions on same line as fout
	echo -e chrom"\t"start"\t"end"\t"Strand"\t"Info"\t"Covered_Cs > ${wDir}/tmp1_$(basename $fin)
	cat $fin >> ${wDir}/tmp1_$(basename $fin)
	mv ${wDir}/tmp1_$(basename $fin) $fin
	
## 6) Add strand, ID and coverage info to the methylation file ($fout)
	cut -f-6 $fin > ${wDir}/tmp2_$(basename $fin)
	cut -f4- $fout > ${wDir}/tmp3_$(basename $fout)
	paste ${wDir}/tmp2_$(basename $fin) ${wDir}/tmp3_$(basename $fout) > $fout
	rm ${wDir}/tmp*${cont}_${i}_${feature}s.bed
done

## 7) merge Fw and Rv regions methylation results
head -1 ${wDir}/meth_${cont}_${i}_${feature}s.bed > $regions_meth
awk 'FNR!=1' ${wDir}/meth_${cont}_*_${feature}s.bed | sort -k1,1V -k2,2n >> $regions_meth

## Cleanup
rm ${wDir}/covered_${cont}_*_${feature}s.bed	# rm fins
rm $index

## 8) Record average methylation distribution of analysed regions
echo -e Meth"\t"Num_of_${feature}s > $reg_meth_distr
for min in {0..98..2};
do
	echo -e -n ${min}-$(( $min + 2 ))%"\t" >> $reg_meth_distr
	awk -v min=$min 'NR!=1{c=0;s=0;max=(min+10);for(i=7;i<=NF;i++){if($i!="NA"){s+=$i;c++}};av=(s/c); if(av>=min && av<max){print}}' $regions_meth | wc -l >> $reg_meth_distr
done

## 9) Find differentially methylated regions
MEF=0.05	# Define Minor Epiallele Frequency (proportion of samples which need to have differential methylation from the others)
# Define Epiallele Difference (minimum methylation difference to define different epialleles) NB: Use 15 or 20 for CpG and 10 or 15 for CHG and CHH
if [ ${cont} == "CpG" ];then ED=20;else ED=15;fi
awk -v MEF=$MEF -v ED=$ED 'OFS="\t"{if(NR==1){print}else{for(i=7;i<=NF;i++){if($i!="NA"){list[i]=$i}};asort(list);
spls=length(list); MEC=int(MEF*spls+0.999); MECdiff=(list[(spls-MEC+1)]-list[MEC]);
if(MECdiff>ED){print}; delete list}}' $regions_meth > ${cont}_${ED}ED_${feature}s.bed

## 10) Extract average methylation of regions with average methylation higher than x%
## NB: Not to count Cs in overlapping genes twice we go back to the unionbed files for calculating the averages
awk -v x=$minMet 'OFS="\t"{if(NR>1){s=0;c=0;for(i=7;i<=NF;i++){if($i!="NA"){s+=$i;c+=1}};av=(s/c);if(av>x){print $1,$2,$3}}}' $regions_meth | bedtools merge -i stdin > ${cont}_${minMet}minmeth_${feature}s_merged.bed	# Filter regions and merge them
head -1 $subset_unionbed | cut -f4- | tr "\t" "\n" > ${cont}_${feature}_tmp1.txt									#Print sample names in a column
tail -n+2 $subset_unionbed | bedtools intersect -a stdin -b ${cont}_${minMet}minmeth_${feature}s_merged.bed | awk '{for(i=4;i<=NF;i++){if($i!="NA"){meth[i]+=$i;c[i]++}}}END{for(i=4;i<=NF;i++){printf "%.2f\n",(meth[i]/c[i])}}' > ${cont}_${feature}_tmp2.txt 
echo -e Sample"\t"${cont}_${feature}s_avmeth_${minMet}minmeth > $feature_avmeth
paste ${cont}_${feature}_tmp1.txt ${cont}_${feature}_tmp2.txt >> $feature_avmeth
rm ${cont}_${feature}_tmp*.txt ${cont}_${minMet}minmeth_${feature}s_merged.bed

