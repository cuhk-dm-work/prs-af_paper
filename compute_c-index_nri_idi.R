#!/bin/bash


#### load the required library

library(survcomp)
library(survival)
library(dplyr)
library(boot)
library(nricens)
library(survIDINRI)


##############################################################################################################################
### Compute c-index in HKDB
##############################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### Fit Cox regression models and compute the predicted value
#### Model 1 (m1): PRS-AF only
#### Model 2 (m2): CHARGE-AF score only
#### Model 3 (m3): PRS-AF and CHARGE-AF score

cox_m1 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ PRS_AF, data=hkdb_data)
pred1 <- predict(cox_m1, type="lp")

cox_m2 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE, data=hkdb_data)
pred2 <- predict(cox_m2, type="lp")

cox_m3 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE + PRS_AF, data=hkdb_data)
pred3 <- predict(cox_m3, type="lp")


#### compute the c-index in HKDB

cindex1 <- concordance.index(x=pred1, surv.time=hkdb_data$AF_TIME_10YR, surv.event=hkdb_data$AF_END_10YR, method="noether")
cindex1$n
cindex1$c.index
cindex1$se
cindex1$lower
cindex1$upper
cindex1$p.value

cindex2 <- concordance.index(x=pred2, surv.time=hkdb_data$AF_TIME_10YR, surv.event=hkdb_data$AF_END_10YR, method="noether")
cindex2$n
cindex2$c.index
cindex2$se
cindex2$lower
cindex2$upper
cindex2$p.value

cindex3 <- concordance.index(x=pred3, surv.time=hkdb_data$AF_TIME_10YR, surv.event=hkdb_data$AF_END_10YR, method="noether")
cindex3$n
cindex3$c.index
cindex3$se
cindex3$lower
cindex3$upper
cindex3$p.value


#### count the sample sizes according to AF event status in HKDB

n_case <- nrow(hkdb_data %>% filter(AF_END_10YR == 1))
n_ctrl <- nrow(hkdb_data %>% filter(AF_END_10YR == 0))
n_case
n_ctrl


### compare the c-index between Model 3 and Model 2 in HKDB

cindex.comp(cindex3, cindex2)




####################################################################################################################################
### Compute 95% CI for c-index in HKDB using bootstrap resampling approach
####################################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### create the array for subsequent analysis

time <- hkdb_data$AF_TIME_10YR
status <- hkdb_data$AF_END_10YR
x1 <- hkdb_data$PRS_AF
x2 <- hkdb_data$CHARGE_AF_SCORE


#### Fit Cox regression models and compute the predicted value
#### Model 1 (m1): PRS-AF only
#### Model 2 (m2): CHARGE-AF score only
#### Model 3 (m3): PRS-AF and CHARGE-AF score

cox_m1 <- coxph(Surv(time, status) ~ x1)
pred1 <- predict(cox_m1, type="lp")

cox_m2 <- coxph(Surv(time, status) ~ x2)
pred2 <- predict(cox_m2, type="lp")

cox_m3 <- coxph(Surv(time, status) ~ x1 + x2)
pred3 <- predict(cox_m3, type="lp")


#### Calculate the original c-index for each model

c_orig1 <- concordance.index(x=pred1, surv.time=time, surv.event=status, method="noether")$c.index

c_orig2 <- concordance.index(x=pred2, surv.time=time, surv.event=status, method="noether")$c.index

c_orig3 <- concordance.index(x=pred3, surv.time=time, surv.event=status, method="noether")$c.index


# Define the bootstrap function

boot_c <- function(data, i) {
  
  #### Fit Cox regression models and compute the predicted value
  #### Model 1 (m1): PRS-AF only
  #### Model 2 (m2): CHARGE-AF score only
  #### Model 3 (m3): PRS-AF and CHARGE-AF score
  cox_m1_boot <- coxph(Surv(time[i], status[i]) ~ x1[i])
  pred1_boot <- predict(cox_m1_boot, type="lp")
  
  cox_m2_boot <- coxph(Surv(time[i], status[i]) ~ x2[i])
  pred2_boot <- predict(cox_m2_boot, type="lp")
  
  cox_m3_boot <- coxph(Surv(time[i], status[i]) ~ x1[i] + x2[i])
  pred3_boot <- predict(cox_m3_boot, type="lp")

  #### calculate the C-index
  c_index_1 <- concordance.index(x=pred1_boot, surv.time=time[i], surv.event=status[i], method="noether")$c.index
  c_index_2 <- concordance.index(x=pred2_boot, surv.time=time[i], surv.event=status[i], method="noether")$c.index
  c_index_3 <- concordance.index(x=pred3_boot, surv.time=time[i], surv.event=status[i], method="noether")$c.index
  
  return(c(c_index_1, c_index_2, c_index_3))
}


# Perform the bootstrap
set.seed(123)
boot_results <- boot(data.frame(time = time, status = status, x1 = x1, x2 = x2), boot_c, R = 10000)


# Calculate 95% CI for the C-index
c_index_ci_m1 <- quantile(boot_results$t[,1], probs = c(0.025, 0.975))
c_index_ci_m2 <- quantile(boot_results$t[,2], probs = c(0.025, 0.975))
c_index_ci_m3 <- quantile(boot_results$t[,3], probs = c(0.025, 0.975))


# View bootstrapped C-index and 95% CI
print(boot_results)
print(c_index_ci_m1)
print(c_index_ci_m2)
print(c_index_ci_m3)





########################################################################################################################################
#### 5-fold cross-validation in HKDB study
#### training model in 4-fold data
#### testing model in 1-fold data
########################################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### 5-fold cross-validation
set.seed(123)
folds <- cut(seq(1, nrow(hkdb_data)), breaks=5, labels=FALSE)

c_index_results <- numeric(5)

for(i in 1:5) {
  test_indices <- which(folds == i, arr.ind = TRUE)
  train_data <- hkdb_data[-test_indices, ]
  test_data <- hkdb_data[test_indices, ]
  
  cox_model <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE + PRS_AF, data = train_data)
  
  risk_scores <- predict(cox_model, newdata = test_data)
  concordance <- concordance.index(x=risk_scores, surv.time=test_data$AF_TIME_10YR, surv.event=test_data$AF_END_10YR, method="noether")$c.index
  
  c_index_results[i] <- concordance
}


# Mean C-index
mean_c_index <- mean(c_index_results)
cat("Mean C-index:", mean_c_index, "\n")





########################################################################################################################################################
### Validate Models 1-3 in HKDR for predicting 10-year AF risk
########################################################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### select the HKDR cohort

hkdr_data <- mydata %>% filter(COHORT == 0)


#### fit cox regression model in HKDB (training cohort)
#### Model 1 (m1): PRS-AF only
#### Model 2 (m2): CHARGE-AF score only
#### Model 3 (m3): PRS-AF and CHARGE-AF score

cox_m1 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ PRS_AF, data=hkdb_data, x=TRUE)
cox_m2 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE, data=hkdb_data, x=TRUE)
cox_m3 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE + PRS_AF, data=hkdb_data, x=TRUE)


#### in HKDR (validation cohort), compute the predicted value for each model estimated in HKDB

pred1 <- predict(cox_m1, hkdr_data, type="lp")
pred2 <- predict(cox_m2, hkdr_data, type="lp") 
pred3 <- predict(cox_m3, hkdr_data, type="lp") 


#### compute the c-index in HKDR

cindex1 <- concordance.index(x=pred1, surv.time=hkdr_data$AF_TIME_10YR, surv.event=hkdr_data$AF_END_10YR, method="noether")
cindex1$n
cindex1$c.index
cindex1$se
cindex1$lower
cindex1$upper
cindex1$p.value

cindex2 <- concordance.index(x=pred2, surv.time=hkdr_data$AF_TIME_10YR, surv.event=hkdr_data$AF_END_10YR, method="noether")
cindex2$n
cindex2$c.index
cindex2$se
cindex2$lower
cindex2$upper
cindex2$p.value

cindex3 <- concordance.index(x=pred3, surv.time=hkdr_data$AF_TIME_10YR, surv.event=hkdr_data$AF_END_10YR, method="noether")
cindex3$n
cindex3$c.index
cindex3$se
cindex3$lower
cindex3$upper
cindex3$p.value


#### count the sample sizes according to AF event status in HKDR

n_case <- nrow(hkdr_data %>% filter(AF_END_10YR == 1))
n_ctrl <- nrow(hkdr_data %>% filter(AF_END_10YR == 0))
n_case
n_ctrl


### compare the c-index between Model 3 and Model 2 in HKDR

cindex.comp(cindex3, cindex2)





####################################################################################################################################
#### Compute 95% CI for c-index in HKDR using bootstrap resampling approach
#### Validate Models 1-3 in HKDR for predicting 10-year AF risk
####################################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### select the HKDR cohort

hkdr_data <- mydata %>% filter(COHORT == 0)


#### fit cox regression model in HKDB (training cohort)
#### Model 1 (m1): PRS-AF only
#### Model 2 (m2): CHARGE-AF score only
#### Model 3 (m3): PRS-AF and CHARGE-AF score

cox_m1 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ PRS_AF, data=hkdb_data, x=TRUE)
cox_m2 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE, data=hkdb_data, x=TRUE)
cox_m3 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE + PRS_AF, data=hkdb_data, x=TRUE)


#### Define the bootstrapping function for the validation dataset

c_index_boot_valid <- function(data, indices) {
  # Create a bootstrapped dataset from the validation data
  boot_data <- data[indices, ]
  boot_data <- boot_data[complete.cases(boot_data[,c("AF_TIME_10YR","AF_END_10YR","CHARGE_AF_SCORE","PRS_AF")]), ]
  
  time <- boot_data$AF_TIME_10YR
  status <- boot_data$AF_END_10YR
  x1 <- boot_data$PRS_AF
  x2 <- boot_data$CHARGE_AF_SCORE
   
  # Predict risk scores using the Cox model fitted on the validation data
  pred1 <- predict(cox_m1, newdata = boot_data)
  pred2 <- predict(cox_m2, newdata = boot_data)
  pred3 <- predict(cox_m3, newdata = boot_data)
  
  # Calculate the C-index
  c_index_1 <- concordance.index(x=pred1, surv.time=time, surv.event=status, method="noether")$c.index
  c_index_2 <- concordance.index(x=pred2, surv.time=time, surv.event=status, method="noether")$c.index
  c_index_3 <- concordance.index(x=pred3, surv.time=time, surv.event=status, method="noether")$c.index
  
  return(c(c_index_1, c_index_2, c_index_3))
}


# Perform bootstrapping
set.seed(123)
boot_results_valid <- boot(data = hkdr_data, statistic = c_index_boot_valid, R = 10000)


# Calculate 95% CI for the C-index
c_index_ci_valid_m1 <- quantile(boot_results_valid$t[,1], probs = c(0.025, 0.975))
c_index_ci_valid_m2 <- quantile(boot_results_valid$t[,2], probs = c(0.025, 0.975))
c_index_ci_valid_m3 <- quantile(boot_results_valid$t[,3], probs = c(0.025, 0.975))


# View bootstrapped C-index and 95% CI
print(boot_results_valid)
print(c_index_ci_valid_m1)
print(c_index_ci_valid_m2)
print(c_index_ci_valid_m3)




####################################################################################################################################
#### compare two models in HKDR using NRI and IDI (for predicting 10-year AF risk) 
#### Model 1 (m1): CHARGE-AF score only
#### Model 2 (m2): PRS-AF and CHARGE-AF score
####################################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### select the HKDR cohort

hkdr_data <- mydata %>% filter(COHORT == 0)


#### fit cox regression model in HKDB (training cohort)
#### Model 1 (m1): CHARGE-AF score only
#### Model 2 (m2): PRS-AF and CHARGE-AF score

cox_m1 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE, data=hkdb_data, x=TRUE)
cox_m2 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE + PRS_AF, data=hkdb_data, x=TRUE)


#### in HKDR (validation cohort), compute the predicted value for each model estimated in HKDB

pred1 <- predict(cox_m1, hkdr_data, type="lp")
pred2 <- predict(cox_m2, hkdr_data, type="lp") 


#### create the array for subsequent analysis

hkdr_time <- hkdr_data$AF_TIME_10YR
hkdr_status <- hkdr_data$AF_END_10YR
pheno <- hkdr_data[,c("AF_TIME_10YR","AF_END_10YR")]

covs0 <- as.matrix(pred1)    # old model
covs1 <- as.matrix(pred2)  # new model


#### compute IDI and NRI
#### 10-year risk
#### no. of iteration

t10_n10000 <- IDI.INF(pheno, covs0, covs1, t0 = 10, npert=10000)


#### results

IDI.INF.OUT(t10_n10000)





########################################################################################################################################################
### Validate Models 1-3 in HKDR for predicting 29-year AF risk
########################################################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "AF_TIME_29YR" = time to event (censoring time is 29 years)
#### "AF_END_29YR" = indictor of AF event during the 29-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","AF_TIME_29YR","AF_END_29YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### select the HKDR cohort

hkdr_data <- mydata %>% filter(COHORT == 0)


#### fit cox regression model in HKDB (training cohort)
#### Model 1 (m1): PRS-AF only
#### Model 2 (m2): CHARGE-AF score only
#### Model 3 (m3): PRS-AF and CHARGE-AF score

cox_m1 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ PRS_AF, data=hkdb_data, x=TRUE)
cox_m2 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE, data=hkdb_data, x=TRUE)
cox_m3 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE + PRS_AF, data=hkdb_data, x=TRUE)


#### in HKDR (validation cohort), compute the predicted value for each model estimated in HKDB

pred1 <- predict(cox_m1, hkdr_data, type="lp")
pred2 <- predict(cox_m2, hkdr_data, type="lp") 
pred3 <- predict(cox_m3, hkdr_data, type="lp") 


#### compute the c-index in HKDR

cindex1 <- concordance.index(x=pred1, surv.time=hkdr_data$AF_TIME_29YR, surv.event=hkdr_data$AF_END_29YR, method="noether")
cindex1$n
cindex1$c.index
cindex1$se
cindex1$lower
cindex1$upper
cindex1$p.value

cindex2 <- concordance.index(x=pred2, surv.time=hkdr_data$AF_TIME_29YR, surv.event=hkdr_data$AF_END_29YR, method="noether")
cindex2$n
cindex2$c.index
cindex2$se
cindex2$lower
cindex2$upper
cindex2$p.value

cindex3 <- concordance.index(x=pred3, surv.time=hkdr_data$AF_TIME_29YR, surv.event=hkdr_data$AF_END_29YR, method="noether")
cindex3$n
cindex3$c.index
cindex3$se
cindex3$lower
cindex3$upper
cindex3$p.value


#### count the sample sizes according to AF event status in HKDR

n_case <- nrow(hkdr_data %>% filter(AF_END_29YR == 1))
n_ctrl <- nrow(hkdr_data %>% filter(AF_END_29YR == 0))
n_case
n_ctrl


### compare the c-index between Model 3 and Model 2 in HKDR

cindex.comp(cindex3, cindex2)





####################################################################################################################################
#### Compute 95% CI for c-index in HKDR using bootstrap resampling approach
#### Validate Models 1-3 in HKDR for predicting 29-year AF risk
####################################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "AF_TIME_29YR" = time to event (censoring time is 29 years)
#### "AF_END_29YR" = indictor of AF event during the 29-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","AF_TIME_29YR","AF_END_29YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### select the HKDR cohort

hkdr_data <- mydata %>% filter(COHORT == 0)


#### fit cox regression model in HKDB (training cohort)
#### Model 1 (m1): PRS-AF only
#### Model 2 (m2): CHARGE-AF score only
#### Model 3 (m3): PRS-AF and CHARGE-AF score

cox_m1 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ PRS_AF, data=hkdb_data, x=TRUE)
cox_m2 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE, data=hkdb_data, x=TRUE)
cox_m3 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE + PRS_AF, data=hkdb_data, x=TRUE)


#### Define the bootstrapping function for the validation dataset

c_index_boot_valid <- function(data, indices) {
  # Create a bootstrapped dataset from the validation data
  boot_data <- data[indices, ]
  boot_data <- boot_data[complete.cases(boot_data[,c("AF_TIME_29YR","AF_END_29YR","CHARGE_AF_SCORE","PRS_AF")]), ]
  
  time <- boot_data$AF_TIME_29YR
  status <- boot_data$AF_END_29YR
  x1 <- boot_data$PRS_AF
  x2 <- boot_data$CHARGE_AF_SCORE
   
  # Predict risk scores using the Cox model fitted on the validation data
  pred1 <- predict(cox_m1, newdata = boot_data)
  pred2 <- predict(cox_m2, newdata = boot_data)
  pred3 <- predict(cox_m3, newdata = boot_data)
  
  # Calculate the C-index
  c_index_1 <- concordance.index(x=pred1, surv.time=time, surv.event=status, method="noether")$c.index
  c_index_2 <- concordance.index(x=pred2, surv.time=time, surv.event=status, method="noether")$c.index
  c_index_3 <- concordance.index(x=pred3, surv.time=time, surv.event=status, method="noether")$c.index
  
  return(c(c_index_1, c_index_2, c_index_3))
}


# Perform bootstrapping
set.seed(123)
boot_results_valid <- boot(data = hkdr_data, statistic = c_index_boot_valid, R = 10000)


# Calculate 95% CI for the C-index
c_index_ci_valid_m1 <- quantile(boot_results_valid$t[,1], probs = c(0.025, 0.975))
c_index_ci_valid_m2 <- quantile(boot_results_valid$t[,2], probs = c(0.025, 0.975))
c_index_ci_valid_m3 <- quantile(boot_results_valid$t[,3], probs = c(0.025, 0.975))


# View bootstrapped C-index and 95% CI
print(boot_results_valid)
print(c_index_ci_valid_m1)
print(c_index_ci_valid_m2)
print(c_index_ci_valid_m3)




####################################################################################################################################
#### compare two models in HKDR using NRI and IDI (for predicting 29-year AF risk) 
#### Model 1 (m1): CHARGE-AF score only
#### Model 2 (m2): PRS-AF and CHARGE-AF score
####################################################################################################################################

#### read the data
#### "AF_TIME_10YR" = time to event (censoring time is 10 years)
#### "AF_END_10YR" = indictor of AF event during the 10-year follow-up (0: no; 1: yes)
#### "AF_TIME_29YR" = time to event (censoring time is 29 years)
#### "AF_END_29YR" = indictor of AF event during the 29-year follow-up (0: no; 1: yes)
#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)
#### "PRS_AF" = polygenic risk score for AF
#### "CHARGE_AF_SCORE" = CHARGE AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")


#### Select subjects for whom complete data is available for all variables in the Cox regression model

mydata <- mydata[complete.cases(mydata[,c("AF_TIME_10YR","AF_END_10YR","AF_TIME_29YR","AF_END_29YR","CHARGE_AF_SCORE","PRS_AF")]), ]


#### select the HKDB cohort

hkdb_data <- mydata %>% filter(COHORT == 1)


#### select the HKDR cohort

hkdr_data <- mydata %>% filter(COHORT == 0)


#### fit cox regression model in HKDB (training cohort)
#### Model 1 (m1): CHARGE-AF score only
#### Model 2 (m2): PRS-AF and CHARGE-AF score

cox_m1 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE, data=hkdb_data, x=TRUE)
cox_m2 <- coxph(Surv(AF_TIME_10YR, AF_END_10YR) ~ CHARGE_AF_SCORE + PRS_AF, data=hkdb_data, x=TRUE)


#### in HKDR (validation cohort), compute the predicted value for each model estimated in HKDB

pred1 <- predict(cox_m1, hkdr_data, type="lp")
pred2 <- predict(cox_m2, hkdr_data, type="lp") 


#### create the array for subsequent analysis

hkdr_time <- hkdr_data$AF_TIME_29YR
hkdr_status <- hkdr_data$AF_END_29YR
pheno <- hkdr_data[,c("AF_TIME_29YR","AF_END_29YR")]

covs0 <- as.matrix(pred1)    # old model
covs1 <- as.matrix(pred2)  # new model


#### compute IDI and NRI
#### 29-year risk
#### no. of iteration

t29_n10000 <- IDI.INF(pheno, covs0, covs1, t0 = 29, npert=10000)


#### results

IDI.INF.OUT(t29_n10000)






