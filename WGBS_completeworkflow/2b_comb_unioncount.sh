#!/bin/bash
## Author: Dario Galanti June 2022
## Aim: Combine scaffold unioncount files from 2a_unioncount_binom_byscaff.sh and extract average and weighted methylation (Schultz et al. 2012)
## Documentation (EpiDiverse WGBS pipeline): https://github.com/EpiDiverse/wgbs
## Input: individual scaffold unioncount files with "methylated/total" read counts (possibly with NAs depending on filtNAs parameter)
## Output 1): unique unioncount file with "methylated/total" read counts (possibly with NAs depending on filtNAs parameter)
## Output 2): txt file with average and weighted methylation (Schultz et al. 2012)
## Run: bash 2b_comb_unioncount.sh indir outdir
## Run: sbatch --partition test --cpus-per-task 4 --mem 40G --time 24:00:00 --wrap "bash 2b_comb_unioncount.sh /scr/episan/RP07/bam_aligned/WGBS/unioncounts_binom_v5_3cov_25NAs/Scaff_res_CpG /scr/episan/RP07/bam_aligned/WGBS/unioncounts_binom_v5_3cov_25NAs"

## FILE DIRECTORIES AND NA FILTERING:
inDir=$1
outDir=$2
#filtNAs=25
cont=$(basename $inDir | grep -o 'C[pH][HG]')
fout=${outDir}/${cont}_unionbinom_v5_25NAs.bed
av_meth_fout=${outDir}/${cont}_binom_avgMet_25NAs.txt

scaff_files=($(ls -d ${inDir}/* | sort -V))
rm $fout			#Fresh start
for f in ${scaff_files[@]};do cat $f >> $fout;done
#rm -r $inDir

## CALCULATE AVG AND WEIGHTED METHYLATION AFTER NA FILTERING
# We output also met and total counts so that weighted methylation can be calculated also accross contexts
head -1 $fout | cut -f4- | tr "\t" "\n" > $outDir/${cont}_tmp1.txt										#Print sample names in a column
#tail -n+2 $fout | awk '{for(i=4;i<=NF;i++){if($i!="NA"){meth[i]+=$i;c[i]++}}}END{for(i=4;i<=NF;i++){printf "%.2f\n",(meth[i]/c[i])}}' > $outDir/${cont}_tmp2.txt
awk 'NR>1{for(i=4;i<=NF;i++){if($i!="NA"){split($i,a,"/");if(a[1]>0){mCs[i]++};metsum[i]+=a[1];totsum[i]+=a[2];AVGsum[i]+=(a[1]/a[2]);c[i]++}}} END{for(i=4;i<=NF;i++){printf "%.2f\t%s\t%s\t%s\t%s\t%.2f\n",(100*AVGsum[i]/c[i]),mCs[i],c[i],metsum[i],totsum[i],(100*metsum[i]/totsum[i])}}' $fout > $outDir/${cont}_tmp2.txt
echo -e Sample"\t"Avg_m${cont}_${filtNAs}NAs"\t"m${cont}s"\t"tot${cont}s"\t"Met_${cont}counts"\t"Tot_${cont}counts"\t"Weighted_m${cont}_${filtNAs}NAs > $av_meth_fout
paste $outDir/${cont}_tmp1.txt $outDir/${cont}_tmp2.txt >> $av_meth_fout
rm $outDir/${cont}_tmp*.txt

