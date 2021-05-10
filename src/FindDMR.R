#!/usr/bin/env Rscript

#Script to take in directory of txt files of CpG WGBS methylation and output DMR, ie. reads in data using Methylkit and calculates DMR using DSS. 
#Txt file format: chr, position, strand (arbitrary), frequency of C, coverage
#Usage: Rscript FindDMR.R <directory name with all cytosine files i.e chromosome>

#Get command line arguments
args = commandArgs(trailingOnly=TRUE)

#Load libraries
library(tidyverse)
library(methylKit)
library(DSS)
require(bsseq)

#Working dir
bs_path <- "/storage1/fs1/ha/Active/maherlab/saha.d/projects/BS_seq/mCRPC/Methylkit_files/"

#List full path for files in directory
bsfiles <- as.list(list.files(path = paste0(bs_path, args[1]), full.names = T))

#Get sample ids
sampleid <- as.list(list.files(path = paste0(bs_path, args[1])))
sampleid <- lapply(str_split(sampleid, ".merged."), FUN = function(i){i[1]})

#Read in clinical annotations df for WGS/WGBS samples and subset sampleids/files/clin list to those in common
clin <- readRDS("/storage1/fs1/ha/Active/maherlab/saha.d/projects/Signature_analysis/mCRPC_WGBS_clin_anno.rds")
bsfiles <- bsfiles[sampleid %in% row.names(clin)]
sampleid <- sampleid[sampleid %in% row.names(clin)]
clin <- clin[row.names(clin) %in% sampleid,]


#Treatment vector
tx <- rep(1, length(sampleid))


#Read in files in methylkit format -keep min coverage at 0 to get more CpGs
obj <- methRead(bsfiles, sample.id=sampleid, treatment = tx,
                assembly="hg38", context="CpG",
                pipeline =list(fraction=FALSE,chr.col=1,start.col=2,end.col=2, coverage.col=5,strand.col=3,freqC.col=4), 
                mincov = 0)

#Format each df per sample for DSS - <chr> <position> <Coverage> <# Methylation>
meth <- lapply(1:length(sampleid), FUN = function(i){
  
  df <- obj[[i]] %>% GRanges() %>% as.data.frame() %>% dplyr::select(c(seqnames, start, coverage, numCs))
  colnames(df) <- c("chr", "pos", "N", "X")
  return(df)  
  
})
names(meth) <- sampleid


#Create BS object
bsobj <- makeBSseqData(dat = meth,sampleNames = names(meth))

#Use clinical anno df as design
DMLfit = DMLfit.multiFactor(bsobj, design=clin, formula=~AR.loci.amplification
+ Type + Site, smoothing = T)

#Test if coef for AR is different from 0 (first term is intercept) for each DML  - return significant DML
DMLtest.cell = DMLtest.multiFactor(DMLfit, coef=2)

#Aggregate significant DML into DMR
dmr <- callDMR(DMLresult = DMLtest.cell, p.threshold = 1e-5)

#Write dmr bed to disk - need to use Aggregate_cpg.R to find methylation values
write.table(x = dmr, file = paste0("/storage1/fs1/ha/Active/maherlab/saha.d/projects/Signature_analysis/", args[1], ".dmr.AR.bed"), sep = "\t", quote = F, row.names = F)

