#!/usr/bin/python3
## Author: Dario Galanti, September 2020
## Aim: Plot the average (and optionally st.dev) distribution of regions in each context at the same time.
## Input: Region unionbed files with methylation values for CpG, CHG and CHH with common name part. The first 6 columns should contain region info ("chr" "start" "end" "strand" "info" "Covered_Cs")
## NB: This script is only good for fairly short unionbed files (100 000 rows max). For long unionbed files see /Bioinformatics/WGBS_DMR_analysis/METHimpute/unionbed_AV_SD_allcont.py
## Run: python3 unionbed_reg_AVG_allcont.py input_file_commonpart feature_name
## Run: sbatch --partition test --cpus-per-task 2 --mem 2G --time 10:00:00 --wrap "python3 unionbed_reg_AVG_allcont.py meth_genes.bed gene"
## Run: NB: Make sure the same files for CHG and CHH are available!!

## Matplotlib code
## https://www.machinelearningplus.com/plots/matplotlib-histogram-python-examples/
## https://matplotlib.org/3.1.1/gallery/statistics/hist.html

## Dependencies: conda activate Python
## NB: This script handles datasets with some NAs! But not huge files, read above.

## Import modules (this modules should be installed!!)
import re
import sys
from statistics import mean
from statistics import pstdev
import matplotlib.pyplot as plt				# To draw SD histogram
from matplotlib import colors
import seaborn as sns
import math									# To use math.ceil()

## Define input files
fin_CpG = str(sys.argv[1])					#CpG unionbed. Make sure corresponding CHG and CHH also exist
fin_CHG = re.sub("CpG", "CHG", fin_CpG)
fin_CHH = re.sub("CpG", "CHH", fin_CpG)
feature = str(sys.argv[2])

## Extract CpG avgs and st.dev lists
line_num = 0
CpG_AV = []		# List of averages to draw histogram
#CpG_SD = []		# List of standard deviations to draw histogram
with open(fin_CpG, "r") as fin:
	for line in fin:
		line_num += 1
		if line_num > 1:
			line = line.strip()
			splitline = line.split("\t")
			val = splitline[6:]
			val = [float(i) for i in val if str(i) != 'NA']
			#sd = round(pstdev(val), 4)
			#CpG_SD.append(sd)
			av = round(mean(val), 4)
			CpG_AV.append(av)

## Extract CHG avgs and st.dev lists
line_num = 0
CHG_AV = []		# List of averages to draw histogram
#CHG_SD = []		# List of standard deviations to draw histogram
with open(fin_CHG, "r") as fin:
	for line in fin:
		line_num += 1
		if line_num > 1:
			line = line.strip()
			splitline = line.split("\t")
			val = splitline[6:]
			val = [float(i) for i in val if str(i) != 'NA']
			#sd = round(pstdev(val), 4)
			#CHG_SD.append(sd)
			av = round(mean(val), 4)
			CHG_AV.append(av)

## Extract CHH avgs and st.dev lists
line_num = 0
CHH_AV = []		# List of averages to draw histogram
#CHH_SD = []		# List of standard deviations to draw histogram
with open(fin_CHH, "r") as fin:
	for line in fin:
		line_num += 1
		if line_num > 1:
			line = line.strip()
			splitline = line.split("\t")
			val = splitline[6:]
			val = [float(i) for i in val if str(i) != 'NA']
			#sd = round(pstdev(val), 4)
			#CHH_SD.append(sd)
			av = round(mean(val), 4)
			CHH_AV.append(av)

print("files iteration finished")

## Store lists into file to download and make graphs locally. This output can be deleted!!!
#with open("Positions_AV_SD.txt", "w") as fout:
#	print("CpG_AV\t" + '\t'.join(map(str,CpG_AV)), file=fout)
#	print("CpG_SD\t" + '\t'.join(map(str,CpG_SD)), file=fout)
#	print("CHG_AV\t" + '\t'.join(map(str,CHG_AV)), file=fout)
#	print("CHG_SD\t" + '\t'.join(map(str,CHG_SD)), file=fout)
#	print("CHH_AV\t" + '\t'.join(map(str,CHH_AV)), file=fout)
#	print("CHH_SD\t" + '\t'.join(map(str,CHH_SD)), file=fout)

## Add a extreme numbers to all lists, to obtain exact bins. Unfortunately bin size cannot be defined
CpG_AV.extend([0,94])
CHG_AV.extend([0,94])
CHH_AV.extend([0,94])

## Print histogram of averages
sns.set_style("white")
kwargs = dict(hist_kws={'alpha':.5}, kde_kws={'linewidth':2})	#The lower alpha the more transparent and blury
plt.figure(figsize=(10,7), dpi= 100)
bins = math.ceil(max(CpG_AV)/2)
sns.distplot(CpG_AV, color="#CD5B45", label="CG", **kwargs, bins=bins, kde=False)	# Old col: "dodgerblue"
bins = math.ceil(max(CHG_AV)/2)
sns.distplot(CHG_AV, color="#EEC900", label="CHG", **kwargs, bins=bins, kde=False)	# Old col: "orange"
bins = math.ceil(max(CHH_AV)/2)
sns.distplot(CHH_AV, color="#1E90FF", label="CHH", **kwargs, bins=bins, kde=False)	# Old col: "deeppink"
plt.tick_params(axis='both', which='major', labelsize=16)
plt.xlabel('Avg methylation %', fontsize=18)
plt.ylabel("".join(['N° ',feature,'s']), fontsize=18)
sns.set_context("paper", rc={"axes.linewidth": 1.2})
#plt.set_yscale('log')
#plt.yscale('log')			# Log the y scale
#plt.xlim(0,100)
plt.legend(fontsize=18)
plt.savefig("".join(['Av_Meth_hist_',feature,'.pdf']), format="pdf")

## OPTIONAL: Print histogram of St deviations (Requires to make SD lists as well while iterating through the files)
#plt.clf()
#sns.set_style("white")
#kwargs = dict(hist_kws={'alpha':.5}, kde_kws={'linewidth':2})	#The lower alpha the more transparent and blury
#plt.figure(figsize=(10,7), dpi= 100)
#bins = math.ceil(max(CpG_AV)/2)
#sns.distplot(CpG_SD, color="#CD5B45", label="CpG", **kwargs, bins=bins, kde=False)	# Old col: "dodgerblue"
#bins = math.ceil(max(CHG_AV)/2)
#sns.distplot(CHG_SD, color="#EEC900", label="CHG", **kwargs, bins=bins, kde=False)	# Old col: "orange"
#bins = math.ceil(max(CHH_AV)/2)
#sns.distplot(CHH_SD, color="#1E90FF", label="CHH", **kwargs, bins=bins, kde=False)	# Old col: "deeppink"
#plt.tick_params(axis='both', which='major', labelsize=16)
#plt.xlabel('St.dev of methylation %', fontsize=18)
#plt.ylabel("".join(['N° ',feature,'s']), fontsize=18)
#plt.xlim(0,0.45)
#plt.legend(fontsize=18)
#plt.savefig('St.dev_meth_hist.pdf', format="pdf")



