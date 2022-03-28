#!/usr/bin/env Rscript

#Script to take in rds dataframe of samples with features and label column - run RF classifier and output results
#Usage: Rscript Run_RF.R <rds df> <column of labels e.g CD8>
#Docker 


#Get command line arguments
args = commandArgs(trailingOnly=TRUE)



#Load libraries
library(tidyverse)
library(ranger)




print("Reading in data")
#Working dir
wd <- "/storage1/fs1/ha/Active/maherlab/saha.d/projects/mCRPC_classifier/"
#Read in dataframe- samples on rows - features on columns with one column for result
data <- readRDS(paste0(wd, "data/", args[1]))
#Set names of features and label columns in df
response <- args[2]
predictors <- data %>% dplyr::select( -c(OS.mCRPC, Event) ) %>% colnames() #Get cell fractions as features


n_features <- length(predictors) 
  


print("Fitting Random survival Forest classifier with default hparams")
rf_surv <- ranger(dependent.variable.name = "OS.mCRPC", 
                  status.variable.name = "Event",data =data, 
                  importance = "permutation", seed = 123 )



print("Writing results")
#Write to disk
saveRDS(rf_surv, file = paste0(wd, "results/rsf_results.rds"))



