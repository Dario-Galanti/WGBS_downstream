#!/usr/bin/env R
## Author: Dario Galanti May 2022
## Aim: Perform binomial test on each CDS, to test if methylation is more dense than average of all CDS.
## Input: unioncount file with methylated/total Cs in each CDS. Whether a C is methylated or not is determinad by a binomial test applied in previous steps.
## Run locally

## References for the method: Takuno and Gaut 2013. Niederhuth et al. 2016.

library(data.table)
library(stringr)
library(stats)

## READ INPUT
#Run one context at a time, extract genes that can be analysed in all contexts and then only carry on analysis with those genes
fin <- "FracMetCs_CDSs/fracMetCpGs_indiv_CDSs.bed"
fin <- "FracMetCs_CDSs/fracMetCHGs_indiv_CDSs.bed"
fin <- "FracMetCs_CDSs/fracMetCHHs_indiv_CDSs.bed"
df <- fread(fin, header=T, data.table=FALSE)  # Whole unionbed with headers
cont <- str_extract(fin, "C[p,H][G,H]")
rownames(df) <- df$Info

# DEFINE Background frequencies (average fraction of met Cs in CDS)
CDS_freqMetCs <- c(0.034755, 0.022100, 0.016805)  #For CpG, CHG and CHH


## FORMAT
#positions <- df[,1:6]
#df <- df[,7:ncol(df)]
spls <- ncol(df)-6
min_cov_spls <- ceiling(spls*0.90)

## FILTERING
## We only analysed genes with more than 10Cs covered in at least 90% of the samples
## We only calculated fraction of methylated Cs in samples with at least 6 citosines covered
## For loop is probably a bit slow, but SO MUCH EASIER SYNTAX then fucking horrible apply and shit R functions
flt_df <- df[1,]
flt_df <- flt_df[-c(1),]
for (r in 1:nrow(df)){
  covered_spls <- 0
  #methylated_spls <- 0
  for(c in 7:ncol(df)){
    if (is.na(df[r,c])==F){
      #met <- as.integer(strsplit(df[r,c], split = "/")[[1]][1])
      tot <- as.integer(strsplit(df[r,c], split = "/")[[1]][2])
      if(tot < 6){
        df[r,c] <- NA
      } else if (tot >= 10) {covered_spls <- (covered_spls + 1)}
    }
  }
  if(covered_spls >= min_cov_spls){
    flt_df <- rbind(flt_df,df[r,])
  }
}
## Or just
flt_df <- df
for (r in 1:nrow(flt_df)){
  for(c in 7:ncol(flt_df)){
    if (is.na(flt_df[r,c])==F){
      tot <- as.integer(strsplit(flt_df[r,c], split = "/")[[1]][2])
      if(tot < 6){
        flt_df[r,c] <- NA
      }
    }
  }
}
## Check that rows with too many NAs were deleted
max(rowSums(is.na(flt_df)))

CpG_analysed_genes <- flt_df$Info
CHG_analysed_genes <- flt_df$Info
CHH_analysed_genes <- flt_df$Info
analysed_genes <- intersect(intersect(CpG_analysed_genes, CHG_analysed_genes), CHH_analysed_genes)
write.table(CHH_analysed_genes,"CHH_analysed_genes.txt", quote = F, row.names = F, col.names = F, sep = "\t")
#NB: Before filtering we have 23187 CDSs in CpG and 23335 in CHG
# After filtering cov10 in at least 90% samples, we have 22960 CDSs in CpG, 23103 in CHG and 23720 in CHH.
# CpG-CHG intersection is 22705 CDSs. CpG-CHG-CHH intersection is 22703 CDSs.

## INTERSECT GENES WITH GENES TO ANALYSE
vec <- rownames(flt_df) %in% analysed_genes
flt_df <- flt_df[vec,]


## DEFINE FUNCTION: Binom test
binom_test <- function(m, prob) {
  if (is.na(m)==F){
    met <- as.integer(strsplit(m, split = "/")[[1]][1])
    tot <- as.integer(strsplit(m, split = "/")[[1]][2])
    if(met > 0){
      b.test <- binom.test(met,tot,p=prob, alternative="greater")
      m <- b.test$p.value
    } else {m <- 1}
  }
  return(m)
}

## FOR LOOP THROUGH COLUMNS AND LAPPLY
gbM_df <- flt_df
print(paste("Beginning binomial tests at", Sys.time()))
for (c in 7:ncol(df)){
  newcol.p <- unlist(lapply(flt_df[,c], binom_test, prob=freq))
  newcol.q <- p.adjust(newcol.p, "BH")
  newcol <- ifelse(newcol.q < 0.05, 1, 0)
  gbM_df[c] <- newcol
}
print(paste("Binomial tests finished at", Sys.time()))



## PRINT RESULTS
fwrite(gbM_df, paste(cont,"_gbM_binom_fdr.txt",sep=""), sep="\t",row.names=F, col.names=T, quote=F, na="NA")

## EXTRACT GENES METHYLATED IN X PORTION OF SAMPLES
gbM_df <- fread(paste(cont,"_gbM_binom_fdr.bed",sep=""), header=T, data.table=F)
rownames(gbM_df) <- gbM_df$Info
gbM_df <- gbM_df[7:ncol(gbM_df)]
gbM_df[] <- lapply(gbM_df, as.integer)


metgenes <- gbM_df[(rowSums(gbM_df, na.rm=T)/spls-rowSums(is.na(gbM_df))) > 0.7,]
write.table(rownames(metgenes),paste(cont,"_met_genes_70perc.txt",sep=""), quote = F, row.names = F, col.names = F, sep = "\t")

varmetgenes <- gbM_df[(rowSums(gbM_df, na.rm=T)/spls-rowSums(is.na(gbM_df))) > 0.3,]
varmetgenes <- varmetgenes[(rowSums(varmetgenes, na.rm=T)/spls-rowSums(is.na(varmetgenes))) < 0.7,]
write.table(rownames(varmetgenes),paste(cont,"_var_met_genes_30-70perc.txt",sep=""), quote = F, row.names = F, col.names = F, sep = "\t")

### Calculate number of methylated genes per sample
spls_metgenes <- as.data.frame(colSums(gbM_df[7:ncol(gbM_df)], na.rm=T))
spls_metgenes$Sample <- rownames(spls_metgenes)
colnames(spls_metgenes) <- c(paste(cont,"_metgenes",sep=""),"Sample")
spls_metgenes <- spls_metgenes[,c(2,1)]
write.table(spls_metgenes,paste(cont,"_metgenes_persample.txt",sep=""), quote = F, row.names = F, col.names = T, sep = "\t")

### Calculate gbM genes per sample
CpG_metCDSs <- fread("CpG_gbM_binom_fdr.bed", header=T, data.table=FALSE)
CHG_metCDSs <- fread("CHG_gbM_binom_fdr.bed", header=T, data.table=FALSE)
CHH_metCDSs <- fread("CHH_gbM_binom_fdr.bed", header=T, data.table=FALSE)

spls <- colnames(CpG_metCDSs)[7:ncol(CpG_metCDSs)]
gbM_per_sample <- rep(0,length(spls))

for (r in 1:nrow(CpG_metCDSs)){
  for(c in 7:ncol(CpG_metCDSs)){
    try(if (CpG_metCDSs[r,c] == 1 && (CHG_metCDSs[r,c] == 0 || is.na(CHG_metCDSs[r,c])) && (CHH_metCDSs[r,c] == 0 || is.na(CHH_metCDSs[r,c]))){
      if(CHG_metCDSs[r,c] == 0 || CHH_metCDSs[r,c] == 0){
        spl_n <- (c-6)
        gbM_per_sample[spl_n] <- gbM_per_sample[spl_n]+1 
      }
    })
  }
}

gbM_persample <- cbind(spls,gbM_per_sample)
colnames(gbM_persample) <- c("id","gbM_genes")
write.table(gbM_persample,"gbM_genes_CGonly_per_sample.txt", quote = F, row.names = F, col.names = T, sep = "\t")

### Calculate TEm genes per sample
spls <- colnames(CpG_metCDSs)[7:ncol(CpG_metCDSs)]
TEm_per_sample <- rep(0,length(spls))

for (r in 1:nrow(CpG_metCDSs)){
  for(c in 7:ncol(CpG_metCDSs)){
    try(if (CHG_metCDSs[r,c] == 1 || CHH_metCDSs[r,c] == 1 ){
        spl_n <- (c-6)
        TEm_per_sample[spl_n] <- TEm_per_sample[spl_n]+1 
    })
  }
}

TEm_persample <- cbind(spls,TEm_per_sample)
colnames(TEm_persample) <- c("id","TEm_genes")
write.table(TEm_persample,"TEm_genes_per_sample.txt", quote = F, row.names = F, col.names = T, sep = "\t")


