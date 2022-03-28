#!/usr/bin/env Rscript

#Script to take in training data of DNA methylation as features and binary label, run glmnet and output results
#Usage: Rscript Run_glmnet.R <rds dataframe>
#Docker: docker(dsaha0295/glmnet)

#Get command line arguments
args = commandArgs(trailingOnly=TRUE)



#Load libraries
library(dplyr)#Data wrangling
library(glmnet)#Regression modeling
library(caret)#Model training
library(doParallel)#Parallelization
library(survival)


cl <- makePSOCKcluster(5)#Specify 5 cores
registerDoParallel(cl)

set.seed(123)


print("Reading in data")
#Working dir
wd <- "/storage1/fs1/ha/Active/maherlab/saha.d/projects/mCRPC_classifier/"

#Read in dataframe- samples on rows - features on columns with one column for result
data <- readRDS(paste0(wd, "data/", args[1]))


#Set names of features and label columns in df
x <- data[,c(celltypes, "Bone")] %>% as.matrix()


#Survival object as response
y  <- with(data, Surv(time = OS.mCRPC, event = Event))



print("Fitting Model")

#Run elastic net logistic regression with 10 fold CV for grid search over lambda and alpha parameters

cvfit <- cv.glmnet(x, y, family = "cox", type.measure = "C")



stopCluster(cl)

print("Writing results")
#Write results to disk
saveRDS(object = cvfit, file = paste0(wd, "results/coxph.res.rds" ))

