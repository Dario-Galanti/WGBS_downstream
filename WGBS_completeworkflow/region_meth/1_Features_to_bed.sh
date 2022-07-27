#!/bin/bash
### Aim: Using Thlaspi .gff annotations (Nunn et al 2022), extract whole genes, CDSs, introns, promoters and TEs to bed files (split to Fw and Rv, and all together)
### NB: This script is specific for the Thlaspi annotations (Nunn et al 2022), it might need to be adapted for other annotations
### Author: Dario Galanti June 2022
### Input: gene and TE annotation files (/scr/episan/RP03/thlaspi/REVISION/FINAL/final.gff) v5, genome index
### Output: 2 bed files (Fw and Rv genes) with no headers. Columns are "chr" "start" "end" "strand" "geneID"
### Run: bash 1_Features_to_bed.sh outdir
### Run: bash 1_Features_to_bed.sh /scr/episan/RP07/region_meth/Features_meth_v5/feature_beds

## Dependencies: python3 and get_promoters.py script

### IMPORTANT TO KEEP IN MIND !!!
## Start position in gff is 1 based, in bed is 0 based!!!
## There are overlapping genes on the two different strands, but not on the same strand!

## DEFINE INPUT FILE AND OUTDIR
python3=~/conda/Python/bin/python3
get_promoters=/scr/episan/RP07/region_meth/Features_meth_v5/get_promoters.py
genome_index=/scr/episan/RP03/thlaspi/REVISION/FINAL/final.fasta.fai
annotation=/scr/episan/RP03/thlaspi/REVISION/FINAL/final.gff				# Nunn et al 2022
TE_anno=/scr/episan/RP03/thlaspi/REVISION/FINAL/final.TEs.gff				# TE annotation Nunn et al. 2022
outDir=$1
mkdir -p $outDir

## DEFINE OUTPUT FILES
all_CDSs=${outDir}/Ta_v5_CDSs.bed
Fw_CDSs=${outDir}/Ta_v5_Fw_CDSs.bed
Rv_CDSs=${outDir}/Ta_v5_Rv_CDSs.bed

all_genes=${outDir}/Ta_v5_genes.bed
Fw_genes=${outDir}/Ta_v5_Fw_genes.bed
Rv_genes=${outDir}/Ta_v5_Rv_genes.bed

introns=${outDir}/Ta_v5_introns.bed		# Regions which are always introns (eg. if there is an exon or UTR on the opposite strand, the region is excluded)
UTRs=${outDir}/Ta_v5_UTRs.bed			# All UTR regions, even if exon on the other strand

all_TEs=${outDir}/Ta_v5_TEs.bed
Fw_TEs=${outDir}/Ta_v5_Fw_TEs.bed
Rv_TEs=${outDir}/Ta_v5_Rv_TEs.bed

all_promoters=${outDir}/Ta_v5_promoters.bed
Fw_promoters=${outDir}/Ta_v5_Fw_promoters.bed
Rv_promoters=${outDir}/Ta_v5_Rv_promoters.bed

intergenic=${outDir}/Ta_v5_intergenic.bed

stats=${outDir}/stats.txt

## EXTRACT CDSs
awk -F"\t|=" 'OFS="\t"{if($3=="CDS" && $2=="T_arvense_v2"){split($10,id,"-");print $1,($4-1),$5,$7,id[1]}}' $annotation | sort -k1,1V -k2,2n > $all_CDSs
cat $all_CDSs | grep + > $Fw_CDSs
cat $all_CDSs | grep - > $Rv_CDSs

## EXTRACT GENES
awk -F"\t|;" 'OFS="\t"{if($3=="gene" && $2=="T_arvense_v2"){split($9,id,"=");print $1,($4-1),$5,$7,id[2]}}' $annotation | sort -k1,1V -k2,2n > $all_genes
cat $all_genes | grep + > $Fw_genes
cat $all_genes | grep - > $Rv_genes

## OBTAIN INTRONS
grep UTR $annotation | awk -F"\t|=" 'OFS="\t"{if($2=="T_arvense_v2"){split($10,id,"-");print $1,($4-1),$5,$7,id[1]}}' | sort -k1,1V -k2,2n > $UTRs
bedtools subtract -a $all_genes -b $all_CDSs | bedtools subtract -a stdin -b $UTRs | sort -k1,1V -k2,2n > $introns

## EXTRACT TEs
## NB: TEs are overlapping and nested even on the same strand, so we split the strands first, and merge all strand specific TEs for further analysis
awk '{if($7 == "+"){OFS="\t";split($9,id,"=");sf=substr(id[3],1,3);print $1,($4-1),$5,$7,sf,id[3]}}' $TE_anno | bedtools merge -i stdin -d -1 -c 4,5,6 -o distinct,distinct,distinct | sort -k1,1V -k2,2n > $Fw_TEs
awk '{if($7 == "-"){OFS="\t";split($9,id,"=");sf=substr(id[3],1,3);print $1,($4-1),$5,$7,sf,id[3]}}' $TE_anno | bedtools merge -i stdin -d -1 -c 4,5,6 -o distinct,distinct,distinct | sort -k1,1V -k2,2n > $Rv_TEs
cat $Fw_TEs $Rv_TEs | sort -k1,1V -k2,2n > $all_TEs

## EXTRACT PROMOTERS
$python3 $get_promoters $Fw_genes $Rv_genes $genome_index $outDir
cat $Fw_promoters $Rv_promoters | sort -k1,1V -k2,2n > $all_promoters

## EXTRACT INTERGENIC REGIONS
## First we combine everything and then subtract it from the whole genome_index
cat $all_genes $all_TEs $all_promoters | cut -f1-3 | sort -k1,1V -k2,2n | bedtools merge -i stdin | bedtools complement -i stdin -g $genome_index > $intergenic

## STATS
echo -e file"\t"seq_length > $stats
for f in ${outDir}/*.bed;
do
 file=$(basename $f)
 bedtools merge -i $f | awk -v f=$file 'OFS="\t"{len+=($3-$2)} END {print f,len}' >> $stats
done


