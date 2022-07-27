#!/bin/bash
### Author: Dario Galanti May 2022
### Aim: Correct for non-conversion rates as described in Schultz et al 2013. Perform binomial test on each C to test if methylation is higher than background frequency (non conv rate) and convert meth counts to zero if binomial test fails.
### Input 1): position-sorted unioncount file (sort -k1,1V -k2,2n), with "methylated/total" read counts for all individuals.
### Input 2): samples.txt file with two columns: sample non_conv_rate
### Output: individual scaffold unioncount files which later have to be combined with 3_comb_unioncount.sh
### Run: bash 2a_unioncount_binom_byscaff.sh /scr/episan/RP07/bam_aligned/WGBS/unioncounts_v5_3cov_25NAs/CpG_unionbed_v5_25NAs.bed /scr/episan/RP07/bam_aligned/WGBS/conversion/final/Sc364_unmetRegs_10cov_NonConvRates_final.txt

## PROCESS:
# To speed up the process we parallelize by scaffold and output individual scaffold unioncount files which later have to be combined (3_comb_unioncount.sh)
# Each scaffold is processed by an R script that performs the binomial test
# NB: Give enough time for the scaffold with more positions to finish

## Define tools and scripts in Epi
Rscript=~/conda/R4/bin/Rscript		# Use conda R4 to access libraries
#binom_rscript=/scr/episan/RP07/region_meth/Features_meth_v5/union_binom.R
binom_rscript=/scr/episan/RP07/bam_aligned/WGBS/unioncount_binom.R

## Define input and outputs
unionbed=$1					# /scr/episan/RP07/region_meth/Features_meth_v5/feature_unionbeds/CpG_CDSs_unionbed_v5_25NAs.bed
spls_NonConv=$2				# /scr/episan/RP07/bam_aligned/WGBS/conversion/final/Sc364_unmetRegs_10cov_NonConvRates_final.txt
cont=$(basename $unionbed | grep -o 'C[pH][HG]')	# Extract context
outdir=/scr/episan/RP07/bam_aligned/WGBS/unioncounts_binom_v5_3cov_25NAs
scaffdir=${outdir}/Scaff_res_${cont}
mkdir -p $scaffdir
chromosomes=7				# Number of chromosomes or big scaffolds. Only used to assign more resources to these

### Test if IDs are the same
union_spls=($(head -1 $unionbed | cut -f4- | tr "\t" "\n"))
NonConv_spls=($(cat $spls_NonConv | cut -f1))
spls=$(echo "${#NonConv_spls[@]}")
for i in $(seq 1 $spls);
do
 if [ "${union_spls[$i]}" != "${NonConv_spls[$i]}" ];then echo ERROR: check sample names and their order; exit 1;fi
done
### Save header
head -1 $unionbed > ${scaffdir}/Scaffold_0_${cont}_header.bed

### SPLIT UNIONBED BY SCAFFOLD
scaff_arr=($(tail -n+2 $unionbed | cut -f1 | uniq))

for scaf in ${scaff_arr[@]};
do
	scaf_union=${scaffdir}/${scaf}_$(basename $unionbed)
	fout=${scaffdir}/${scaf}_binom_$(basename $unionbed)
	## Define resources for big and small scaffolds
	scaf_num=$(echo $scaf | grep -o -E '[0-9]+')
	if [ $scaf_num -le $chromosomes ];then G=20; h=120;else G=8; h=20;fi
	jobName=/scr/episan/RP07/work/${cont}_binom_${scaf}.sh
	(
	echo "#!/bin/bash"
	echo "#SBATCH --partition test"		#crunch -> epi or test -> diverse
	echo "#SBATCH --cpus-per-task 2"
	echo "#SBATCH --mem ${G}G"
	echo "#SBATCH --time ${h}:00:00"		# 1M lines in about 1:40min
	echo "#SBATCH --output /scr/episan/RP07/logs/${cont}_binom_${scaf}.%A.out"
	echo "#SBATCH --error /scr/episan/RP07/logs/${cont}_binom_${scaf}.%A.err"
	echo "#SBATCH --job-name ${cont}_binom_${scaf}"
	echo ""
	echo "tail -n+2 $unionbed | grep -P ${scaf}'\t' > $scaf_union"
	echo "$Rscript $binom_rscript ${scaf_union} ${spls_NonConv} $fout"
	echo "rm ${scaf_union}"
	) > $jobName
		chmod +x $jobName
		echo "Bash file $jobName created"
		sbatch ${jobName}
done
exit

### COMBINE SCAFFOLD RESULTS
#unionbed=feature_unionbeds/CpG_CDSs_unionbed_v5_25NAs.bed
#cont=$(basename $unionbed | grep -o 'C[pH][HG]')
#fout=feature_union_binom/${cont}_CDSs_unionbinom_v5_25NAs.bed
#scaff_files=($(ls feature_union_binom/Scaff_results/Scaffold_*_binom_$(basename $unionbed) | sort -V))
#head -1 ${unionbed} > $fout
#for f in ${scaff_files[@]};do cat $f >> $fout;done

### CALCULATE Methylated and Total Cs in CDSs after the binomial test
#echo -e Sample"\t"m${cont}s"\t"tot${cont}s > feature_union_binom/FreqMet${cont}s_CDSs_binom_v5_25NAs.txt
#head -1 $fout | cut -f4- | tr "\t" "\n" > ${cont}_tmp1.txt
#tail -n+2 $fout | awk '{for(i=4;i<=NF;i++){if($i!="NA"){split($i,a,"/");if(a[1]!=0){mCs[i]++};totCs[i]++}}} END{for(i=4;i<=NF;i++){printf "%s\t%s\n",mCs[i],totCs[i]}}' > ${cont}_tmp2.txt
#paste ${cont}_tmp1.txt ${cont}_tmp2.txt >> feature_union_binom/FreqMet${cont}s_CDSs_binom_v5_25NAs.txt
#rm ${cont}_tmp1.txt ${cont}_tmp2.txt



