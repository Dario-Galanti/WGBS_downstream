#!/usr/bin/python3
## Author: Dario Galanti, Nov 2020
## Input 1: Two bed files of forward and reverse genes
## Input 2: Ref genome index file to retrieve scaffold lengths
## Aim: Extract promoter regions in bed format. 2kb are used by default but promoters are cut shorter if overlapping with upstream genes or if exiting scaffold limits
## Aim: In addition draw a distribution of promoter lengths
## Output: 2 bed files (one per strand) with promoter regions in format: chromosome"\t"start"\t"end"\t"strand"\t"geneID"\t"length
## Run: python3 get_promoters.py Fw_genes.bed Rv_genes.bed /scr/epi/genomes/thlaspi_arvense/thlaspi.fa.fai

## Dependencies: conda activate Python

## Import modules
import sys
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

## Define input
fin_Fw = str(sys.argv[1])
fin_Rv = str(sys.argv[2])
index = str(sys.argv[3])

## Populate index dictionary (Scaffold:length)
genome_dict = {}
with open(index,"r") as ind:
	for line in ind:
		line = line.strip()
		line = line.split("\t")
		genome_dict[line[0]] = line[1]

## Define output files and list to store promoter lengths
Fw_prom = "Fw_promoters.bed"
Rv_prom = "Rv_promoters.bed"
prom_len_hist = "Prom_length_dist.png"
prom_len = []

prev_gene = ["start"]
with open(fin_Fw, "r") as fin, open(Fw_prom, "w") as fout:	# Fw strand (5'->3'), so extract 2kb before gene (upstream)
	for line in fin:
		line = line.strip()
		line = line.split("\t")
	## Retrieve previous gene end if was on same scaffold
		if prev_gene[0] == line[0]:
			prev_gene_end = int(prev_gene[2])
		else:
			prev_gene_end = 0
	## Define promoter region
		st = max(prev_gene_end,(int(line[1])-2000))	#Input is in bed format, so start is 0 based
		end = int(line[1])							#Input is in bed format, so start is 0 based -> no need to subtract one extra base
		id = str(line[4])
		length = end-st
		if length > 0:								# Just in case there are overlapping genes on the same strand
			prom_len.append(length)
			print(line[0],st,end,line[3],id,length, sep="\t", file=fout)
	## Store previous gene info
		prev_gene = line[0:3]
		
		
scaff = "start"
line_num = 0
with open(fin_Rv, "r") as fin, open(Rv_prom, "w") as fout:	# Rv strand (3'->5'), so extract 2kb after gene (upstream)
	for line in fin:
		line = line.strip()
		line = line.split("\t")
		line_num += 1
	## Skip first line
		if line_num > 1:
	## If current line has same scaffold as previous, cut prev line promoter if overlapping with newline gene.
	## If current line has new scaffold, check whether prev line promoter is not exceeding scaffold length.
			if line[0] == scaff:
				end = min(end, int(line[1]))
			else:
				end = min(end, scaff_end)
			
			length = end-st
			if length > 0:								# Just in case there are overlapping genes on the same strand
				prom_len.append(length)
				print(scaff,st,end,strand,id,length, sep="\t", file=fout)
	## Store line info to print in next iteration
		scaff = line[0]
		scaff_end = int(genome_dict[scaff])					#Retrieve scaffold end
		st = int(line[2])									#Input is in bed format, so end is 1 based but start should be 0 based -> no need to add one base
		end = st+2000
		id = str(line[4])
		strand = line[3]
	
	## Print last iteration
	end = min(end, scaff_end)
	length = end-st
	print(scaff,st,end,strand,id,length, sep="\t", file=fout)


print("Min promoter lenght: " + str(min(prom_len)) + "\t" + "Max promoter lenght: " + str(max(prom_len)))
## Draw histogram of promoter lengths
plt.hist(prom_len, 20, color='#0504aa')
#plt.grid(axis='y', alpha=0.75)
#plt.xlim(0,2000)
plt.xlabel('length(bp)')
plt.ylabel('Count')
plt.title('Promoter length distribution')
#maxfreq = n.max()
#plt.ylim(ymax=np.ceil(maxfreq / 10) * 10 if maxfreq % 10 else maxfreq + 10)

#sns.histplot(data=prom_len, x="Promoter_length(bp)")
plt.savefig(prom_len_hist)
