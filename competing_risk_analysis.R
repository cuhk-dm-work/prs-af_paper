#!/bin/bash


#######################################################################
### Consider competing risks of death using the Fine-Gray models
#######################################################################

#### load the required library

library(cmprsk)


#### read the data
#### "AF_DEATH_TIME" = time to event (year)
#### "AF_DEATH_END" = indictor of AF or mortality events (0: no; 1: AF event; 2: mortality event)
#### "COHORT" = 0 (HKDR)
#### "COHORT" = 1 (HKDB)
#### "YEAR_BASELINE" = enrollment year
#### "PC1","PC2","PC3","PC4" = the first four principal components
#### "SEX" = 1 (male)
#### "SEX" = 2 (female)
#### "AGE" = baseline age
#### "DMAGE" = baseline duration of diabetes
#### "PRS_AF" = polygenic risk score for AF

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")



#### Select subjects for whom complete data is available for all variables in the Cox regression model.

mydata <- mydata[complete.cases(mydata[,c("AF_DEATH_TIME","AF_DEATH_END","COHORT","YEAR_BASELINE","PC1","PC2","PC3","PC4","SEX","AGE","DMAGE","PRS_AF")]), ]



#### select the covariates

covariates = mydata[,c("COHORT" ,"YEAR_BASELINE", "PC1","PC2","PC3","PC4","SEX","AGE","DMAGE","PRS_AF")]



# Fit the Fine-Gray model

fg_model <- crr(ftime = mydata$AF_DEATH_TIME, fstatus = mydata$AF_DEATH_END, cov1 = covariates)

summary_fg = summary(fg_model)

