#!/usr/bin/env R
## Author: Dario Galanti May 2022
## Aim: Perform binomial test on each C to test if methylation is higher than error rate (NonConversionRate) and convert meth counts to zero if binomial test fails.
## Input 1: simil-unionbed file with methylated/total counts
## Input 2: samples.txt: txt file with 2 columns: Samples and NonConvRate
## Run: Rscript --vanilla unioncount_binom.R Scaffold_1_CpG_CDSs_unionbed_v5_5NAs.bed NonConvRates_final.txt Scaffold_1_CpG_CDSs_union_binom.bed 
## Run: sbatch --partition test --cpus-per-task 4 --mem 40G --time 30:00:00 --wrap "Rscript --vanilla unioncount_binom.R /scr/episan/RP07/region_meth/Features_meth_v5/feature_unionbeds/CpG_CDSs_unionbed_v5_25NAs.bed Sc364_unmetRegs_10cov_NonConvRates_final.txt CpG_CDSs_union_binom.bed"

## References for the method: Takuno and Gaut 2013. Niederhuth et al. 2016.

library(data.table)
library(stringr)
library(stats)

## READ INPUT
args <- commandArgs(trailingOnly=T)
#df <- fread(args[1], header=T, data.table=FALSE)  # Whole unionbed with headers
df <- fread(args[1], header=F, data.table=FALSE)  #Scaffold unionbed without headers
NonConvRates <- read.table(args[2], header=F)
fout <- args[3]

## FORMAT
positions <- df[,1:3]
df <- df[,4:ncol(df)]
NonConvRates[2] <- NonConvRates[2]/100

## DEFINE FUNCTION
binom_cleanup <- function(m, NCR) {
  if (is.na(m)==F){
    met <- as.integer(strsplit(m, split = "/")[[1]][1])
    #met <- as.integer(str_extract(m, '^[0-9]+'))  #slower
    tot <- as.integer(strsplit(m, split = "/")[[1]][2])
    #tot <- as.integer(str_extract(m, '[0-9]+$'))  #slower
    if(met > 0){
      b.test <- binom.test(met,tot,p=NCR, alternative="greater")
      if(b.test$p.value > 0.01){m <- paste(0,tot,sep="/")}
    }
  }
  return(m)
}


## FOR LOOP THROUGH COLUMNS AND LAPPLY
print(paste("Beginning binomial tests at", Sys.time()))
for (c in 1:ncol(df)){
  NCRate <- NonConvRates[c,2]
  newcol <- unlist(lapply(df[,c], binom_cleanup, NCR=NCRate))
  df[c] <- newcol
}
print(paste("Binomial tests finished at", Sys.time()))

df <- cbind(positions,df)

## PRINT RESULTS
fwrite(df, fout, sep="\t",row.names=FALSE, col.names=FALSE, quote=FALSE, na="NA")



## OPTIONAL: calculate weighted methylation
#met_calls <- as.integer(strsplit(m, split = "/")[[1]][1])
#substr_met  <- function(x){
#  if (is.na(x)==F){x <- as.integer(strsplit(x, split = "/")[[1]][1])}
#  return(x)
#}
#substr_tot  <- function(x){
#  if (is.na(x)==F){x <- as.integer(strsplit(x, split = "/")[[1]][2])}
#  return(x)
#}
#df_met <- data.frame(lapply(df,substr_met))



