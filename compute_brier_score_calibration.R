#!/bin/bash

########################################
#### Computation of Brier score
########################################

#### load the required library

library(survival)


#### read the data
#### "AF_TIME" = time to event (year)
#### "AF_END" = indictor of AF event (0: no; 1: yes)
#### "PRS_CHARGE_AF_SCORE" = a novel risk score that integrates the PRS-AF into the CHARGE-AF score

mydata <- read.table(file="mydata.txt", header=TRUE, sep="\t")

mydata <- mydata[complete.cases(mydata[,c("AF_TIME","AF_END","PRS_CHARGE_AF_SCORE")]), ]


#### "COHORT" = 0 (HKDR - validation cohort)
#### "COHORT" = 1 (HKDB - model training cohort)    

hkdr_data <- mydata %>% filter(COHORT == 0)

hkdb_data <- mydata %>% filter(COHORT == 1)



# Fit a Cox model in HKDB cohort

cox_model <- coxph(Surv(AFIB_TIME, AFIB_END) ~ PRS_CHARGE_AF_SCORE, data=hkdb_data)



# compute predicted risk from Cox model in HKDB cohort

hkdb_data$pred_risk <- 1 - predict(cox_model, type = "survival", newdata = hkdb_data)



# compute predicted risk from Cox model in HKDR cohort

hkdr_data$pred_risk <- 1 - predict(cox_model, type = "survival", newdata = hkdr_data)


# compute the Brier score in HKDB cohort

brier_score1 <- mean((hkdb_data$pred_risk - hkdb_data$AFIB_END)^2)
print(brier_score1)

# compute the Brier score in HKDR cohort

brier_score2 <- mean((hkdr_data$pred_risk - hkdr_data$AFIB_END)^2)
print(brier_score2)









