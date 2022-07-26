#!/usr/bin/env python

'''
Title: average_over_bed.py
Date: 20191012
Author: Adam Nunn (modified by Dario Galanti)
Original version: https://github.com/EpiDiverse/ewas/blob/master/bin/average_over_bed.py
Description:
	This script iterates over a positions file generated by bedtools unionbedg
	and a corresponding regions file, to give the average values of all positions
	contained within each region. Only works with position sorted bed files with
	distinct, non-overlapping regions, such as the DMRs from Metilene.
List of functions:
	main()
	buildIndex()
	moveToNextRegion()
Procedure:
	1. Open the regions and positions files for reading
	2. Iterate through the positions file
	3. On each position, test if we need to move to the next region
	4. On each position, evaluate or skip to next
Input:
	1) Bedfile harbouring location-ordered, discrete non-overlapping regions
	2) unionbedg file with headers (can contain NAs)
	3) Index file: a txt file with one column, containing all possible scaffolds (sorted) from both files.
Usage:
	average_over_bed.py [REGIONS] [POSITIONS] [INDEX] [FOUT]
eg. average_over_bed.py regions.bed positions.bed index.txt reg_meth.bed
'''

###################
## INIT ENVIRONMENT

import sys


#################
## BEGIN __MAIN__
def main(REGIONS,POSITIONS,INDEX,FOUT):										# Some code explanations

	# 0) Build scaffold index
	index = buildIndex(INDEX)

	# 1) Open the regions and positions files for reading
	with open(REGIONS,'r') as regions, open(POSITIONS,'r') as positions, open(FOUT, 'w') as fout:

		# setup vars for initial region
		try: region = next(regions)											# This try block checks whether there are any lines in the regions file
		except StopIteration: 
			print("No regions detected in bedfile")
			raise SystemExit(0)
		
		region = region.rstrip()
		region = region.split("\t")

		Chr = region[0]
		iChr = index.index(Chr) 											# Retrieve the Chr/Scaff index in the index list
		Start = int(region[1])
		End = int(region[2])

		# define booleans for header and region detection
		Head = True
		Found = False

		# 2) Iterate through the positions file
		for position in positions:

			# split the position line
			position = position.rstrip()
			position = position.split("\t")

			# stage variables and print the header
			if Head:

				# establish initial variables
				length = len(position) - 3
				totals = [(0, 0) for i in range(0,length)]

				# print header and skip to next position
				print("\t".join(position), file=fout)
				Head = False
				continue

			#Chr = region[0]
			#iChr = index.index(Chr)                                                                                         # Retrieve the Chr/Scaff index in the index list
			#Start = int(region[1])
			#End = int(region[2])

			print(totals[0])

			##############################################
			# establish current position
			current_Chr = position[0]
			current_iChr = index.index(current_Chr)
			current_Pos = int(position[1]) # should be START pos

			# 3) On each new position, test if we have reached the end of a Found region
			if Found and ((current_Chr != Chr) or ((current_Chr == Chr) and (current_Pos >= End))):
				#print("{}\t{}\t{}\t{}\t{}\t{}\tmoving to next region".format(Chr,Start,End,current_Chr,current_Pos,position[2]),file=sys.stderr())

				# print averages for current region
				averages = [str(format(n[0]/n[1], '.5f')) if n[1] != 0 else "NA" for n in totals]
				#print(totals)
				#print(averages)
				print("{}\t{}\t{}\t{}".format(Chr, Start, End, "\t".join(averages)), file=fout)

				# reset counts
				totals = [(0, 0) for i in range(0,length)]
				Found = False

				# move the region until we are ahead of current position
				region = moveToNextRegion(index, regions, current_Chr, current_iChr, current_Pos)
				if region is None: break
				else:
					Chr = region[0]
					Start = int(region[1])
					End = int(region[2]) 

			# 4) On each position, evaluate current position or skip to next
			if (current_Chr == Chr) and (current_Pos in range(Start, End)):
				#print("{}\t{}\t{}\t{}\t{}\t{}\tposition is inside region".format(Chr,Start,End,current_Chr,current_Pos,position[2]),file=sys.stderr())
				totals = [(totals[i][0],totals[i][1]) if position[i+3] == "NA" else (totals[i][0]+float(position[i+3]),totals[i][1]+1) for i in range(0, length)]
				Found = True

			# we already moved the region, so now move position
			else:
				#print("{}\t{}\t{}\t{}\t{}\t{}\tsearching for first position in region".format(Chr,Start+1,End+1,current_Chr,current_Pos,position[2]),file=sys.stderr())
				continue

		# unclosed For loop
		if Found:
			averages = [str(format(n[0]/n[1], '.5f')) if n[1] != 0 else "NA" for n in totals]
			print("{}\t{}\t{}\t{}".format(Chr, Start, End, "\t".join(averages)), file=fout)


## END OF __MAIN__
##################


###################
## DEFINE FUNCTIONS

# build index of scaffold IDs
def buildIndex(INDEX):			# This simply returns a list of all lines in index.txt (so all Chrom and Scaffolds)

	with open(INDEX, 'r') as f: lines = f.read().splitlines()
	return lines

# move to the next region until we are ahead of position
def moveToNextRegion(index, regions, Chr, iChr, Pos):

	# move the region
	try: region = next(regions)
	except StopIteration: return None

	region = region.rstrip()
	region = region.split("\t")

	# if the new region is behind the current position, we need to run recursively
	if (index.index(region[0]) < iChr) or ((region[0] == Chr) and (Pos+1 > int(region[2]))):
		#print("{}\t{}\t{}\tregion removed due to coverage filter".format(region[0],region[1],region[2]),file=sys.stderr)
		region = moveToNextRegion(index, regions, Chr, iChr, Pos)
	
	return region


## END OF FUNCTIONS
###################

#############
## RUN SCRIPT

# run main()
if __name__ == '__main__':
	main(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4])

## END OF SCRIPT
################
