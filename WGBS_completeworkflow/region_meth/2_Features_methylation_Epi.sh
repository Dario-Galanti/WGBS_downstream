### Author: Dario Galanti May 2022
### Aim: Extract avg and weighted methylation of whole genomic features (all genes, CDSs, promoters, TEs ...) from unioncount file
### Input 1): position-sorted unioncount file (properly sorted, so "sort -k1,1V -k2,2n") with "methylated/total" read counts per position
### Input 2): Feature bed files from 1_Features_to_bed.sh
### Run: bash 2_Features_methylation_Epi.sh ${context}_unionbed.bed
### Run: bash 2_Features_methylation_Epi.sh /scr/episan/RP07/bam_aligned/WGBS/unioncounts_v5_3cov_25NAs/CpG_unioncount_v5_25NAs.bed

### Steps:
### 1) Loop through feature files
### 2) Subset unioncount to each feature of interest
### 2) Extract average and weighted methylation of each feature (see Schultz et al. 2013 for avg vs weighted calculation)

### IMPORTANT TO KEEP IN MIND !!!
## NB: Feature filenames must contain feature name just before file extention: Ta_v5_${feature}s.bed. Eg. Ta_v5_genes.bed

## DEFINE UNIONBED AND FEATURE FILES
unionbed=$1						# /scr/episan/RP07/bam_aligned/WGBS/unioncounts_binom_v5_3cov_25NAs/CpG_unionbinom_v5_25NAs.bed
inDir=/scr/episan/RP07/region_meth/Features_meth_v5/feature_beds
Genes=${inDir}/Ta_v5_genes.bed					# NB: Feature filenames must contain feature name just before file extention!
CDSs=${inDir}/Ta_v5_CDSs.bed					# NB: Feature filenames must contain feature name just before file extention!
introns=${inDir}/Ta_v5_introns.bed				# NB: Feature filenames must contain feature name just before file extention!
TEs=${inDir}/Ta_v5_TEs.bed						# NB: Feature filenames must contain feature name just before file extention!
promoters=${inDir}/Ta_v5_promoters.bed			# NB: Feature filenames must contain feature name just before file extention!
intergenic=${inDir}/Ta_v5_intergenic.bed		# NB: Feature filenames must contain feature name just before file extention!

## DEFINE OUTPUT, CONTEXT AND WORK FILES
cont=$(echo $unionbed | grep -o 'C[pH][HG]' | cut -f1)					# Extract context
outdir=/scr/episan/RP07/region_meth/Features_meth_v5/featureMet_results
sub_union_dir=/scr/episan/RP07/region_meth/Features_meth_v5/feature_unionbeds
mkdir -p $outdir
mkdir -p $sub_union_dir
mkdir -p /scr/episan/RP07/work
mkdir -p /scr/episan/RP07/logs

## FEATURE FILES ARRAY
feature_files=($Genes $CDSs $introns $TEs $promoters $intergenic)

### LOOP THROUGH FEATURE FILES
for fin in ${feature_files[@]};
do
	feature=$(basename $fin .bed | rev | cut -d_ -f1 | rev)
	jobName=/scr/episan/RP07/work/m${cont}_calc_"$(basename $fin .bed)".sh
	(
	echo "#!/bin/bash"
	echo "#SBATCH --partition crunch"		#crunch -> epi or test -> diverse
	echo "#SBATCH --cpus-per-task 8"
	echo "#SBATCH --mem 32G"
	echo "#SBATCH --time 20:00:00"
	echo "#SBATCH --output /scr/episan/RP07/logs/m${cont}_${feature}_calc.%A.out"
	echo "#SBATCH --error /scr/episan/RP07/logs/m${cont}_${feature}_calc.%A.err"
	echo "#SBATCH --job-name m${cont}_calculation"
	echo ""
	## Define stuff
	echo "subset_unionbed=${sub_union_dir}/${cont}_${feature}_unionbed_v5_25NAs.bed"
	echo "met_fout=${outdir}/m${cont}_${feature}_v5_25NAs.txt"
	echo ""
	## Intersect unionbed
	echo "head -1 $unionbed > \$subset_unionbed"
	echo "tail -n+2 $unionbed | bedtools intersect -a stdin -b $fin -u >> \$subset_unionbed"
	echo ""
	## Calculate average and weighted methylation
	echo "head -1 \$subset_unionbed | cut -f4- | tr '\t' '\n' > ${outdir}/${cont}_${feature}_tmp1.txt"				#Print sample names in a column
	echo "awk 'NR>1{for(i=4;i<=NF;i++){if(\$i!=\"NA\"){split(\$i,a,\"/\");metsum[i]+=a[1];totsum[i]+=a[2];AVGsum[i]+=(a[1]/a[2]);c[i]++}}} END{for(i=4;i<=NF;i++){printf \"%.2f\t%.2f\n\",(100*AVGsum[i]/c[i]),(100*metsum[i]/totsum[i])}}' \$subset_unionbed > ${outdir}/${cont}_${feature}_tmp2.txt"
	echo "echo -e Sample'\t'Avg_m${cont}_${feature}_25NAs'\t'Weighted_m${cont}_${feature}_25NAs > \$met_fout"
	## or Fraction of methylated Cs (only after performing binomial test and converting to 0 met reads of non significant positions)
	#echo "awk 'NR>1{for(i=4;i<=NF;i++){if(\$i!=\"NA\"){split(\$i,a,\"/\");if(a[1]>0){metsum[i]++};totsum[i]++}}} END{for(i=4;i<=NF;i++){printf \"%.2f\n\",(100*metsum[i]/totsum[i])}}' \$subset_unionbed > ${outdir}/${cont}_${feature}_tmp2.txt"
	#echo "echo -e Sample'\t'Fraction_m${cont}_${feature}_25NAs > \$met_fout"
	echo "paste ${outdir}/${cont}_${feature}_tmp1.txt ${outdir}/${cont}_${feature}_tmp2.txt >> \$met_fout"
	echo "rm ${outdir}/${cont}_${feature}_tmp*.txt"
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		sbatch ${jobName}
done
exit

## COMBINE ALL RESULTS
#cat ${outdir}/mCHG_CDSs_v5_25NAs.txt | cut -f1 > ${outdir}/Avg_feature_meth.txt
#for f in ${outdir}/mC*_25NAs.txt;
#do
# cat $f | cut -f2 > ${outdir}/Avg_tmpfile.txt
# cat $f | cut -f3 > ${outdir}/Weighted_tmpfile.txt
# paste ${outdir}/Avg_feature_meth.txt ${outdir}/Avg_tmpfile.txt > ${outdir}/Avg_feature_meth_tmp.txt
# paste ${outdir}/Weighted_feature_meth.txt ${outdir}/Weighted_tmpfile.txt > ${outdir}/Weighted_feature_meth_tmp.txt
# mv ${outdir}/Avg_feature_meth_tmp.txt ${outdir}/AvgMet_features.txt
# mv ${outdir}/Weighted_feature_meth_tmp.txt ${outdir}/WeighMet_features.txt
#done

