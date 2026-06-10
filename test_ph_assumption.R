#!/bin/bash


############################################################################################
#### Assess the proportional hazards assumption using the Schoenfeld residuals test 
############################################################################################

#### load the required library

library(survival)


#### read the data
#### "AF_TIME" = time to event (year)
#### "AF_END" = indictor of AF event (0: no; 1: yes)
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

mydata <- mydata[complete.cases(mydata[,c("AF_TIME","AF_END","COHORT","YEAR_BASELINE","PC1","PC2","PC3","PC4","SEX","AGE","DMAGE","PRS_AF")]), ]


#### Fit Cox regression model

cox_model <- coxph(Surv(AF_TIME, AF_END) ~ COHORT + YEAR_BASELINE + PC1 + PC2 + PC3 + PC4 + SEX + AGE + DMAGE + PRS_AF, data = mydata)


# Test proportional hazards assumption

ph_test <- cox.zph(cox_model)

summary(ph_test)

