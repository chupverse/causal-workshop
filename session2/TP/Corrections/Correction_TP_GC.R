#########################################################################################################
### R Code correction TP causality GC
#########################################################################################################

### PART 1
#########################################################################################################
### 1) a)
#########################################################################################################
# remotes::install_github("chupverse/gcomputation")
library(gcomputation)
data(dataCOHORT)

#########################################################################################################
### 1) b)
#########################################################################################################
dim(dataCOHORT)
dataCOHORT <- subset(dataCOHORT, AGE > 30)
dim(dataCOHORT) #81 patients were removed

#########################################################################################################
### 1) c)
#########################################################################################################
vars <- c("AGE", "BMI", "SEX", "GLASGOW", "INJURY", "PAO2FIO2")

#########################################################################################################
### 2) a)
#########################################################################################################
dataCOHORT$VAP_num <- ifelse(dataCOHORT$VAP == "Yes", 1, 0)
#0 for non-events and 1 for events.
dataCOHORT$GROUP_num <- ifelse(dataCOHORT$GROUP == "Untreated", 0, 1)
#0 for the untreated/unexposed patients and 1 for the treated/exposed ones.

#########################################################################################################
### 2) b)
#########################################################################################################
dataCOHORT_CC <- dataCOHORT[complete.cases(dataCOHORT[,vars]),]
dim(dataCOHORT_CC) #removed 74 patients due to missing values
fit <- glm(formula = VAP_num ~ GROUP_num + AGE + BMI + SEX + GLASGOW + INJURY + PAO2FIO2,
           data = dataCOHORT_CC, family = binomial(link="logit"))

#########################################################################################################
### 3) a)
#########################################################################################################
data0 <- dataCOHORT_CC
data1 <- dataCOHORT_CC
data0$GROUP_num <- 0 #Set everyone to Untreated
data1$GROUP_num <- 1 #Set everyone to Treated

#########################################################################################################
### 3) b)
#########################################################################################################
p0i <- predict(fit, newdata = data0, type="response")
p1i <- predict(fit, newdata = data1, type="response")

#########################################################################################################
### 4) a)
#########################################################################################################
p0 <- mean(p0i)
p1 <- mean(p1i)

#########################################################################################################
### 4) b)
#########################################################################################################
ATE <- p1 - p0

#########################################################################################################
### 5)
#########################################################################################################
# The estimated ATE was -0.15
# Meaning that approximately 15.2 VAP cases were prevented for every 100 treated patients
#  by using Ceftriaxone compared to the untreated group.

#########################################################################################################

### PART 2
#########################################################################################################
### 6) a)
#########################################################################################################
gc_ATE <- function(data, formula, group) {
  fit <- glm(formula = formula, data = data, family = binomial(link="logit"))
  
  data0 <- data
  data1 <- data
  data0[,group] <- 0
  data1[,group] <- 1
  
  p0i <- predict(fit, newdata = data0, type="response")
  p1i <- predict(fit, newdata = data1, type="response")
  
  ATE <- mean(p1i) - mean(p0i)
  return(ATE)
}
gc_ATE(data = dataCOHORT_CC,
       formula = formula(VAP_num ~ GROUP_num + AGE + BMI + SEX +
                           GLASGOW + INJURY + PAO2FIO2),
       group = "GROUP_num")

#########################################################################################################
### 6) b)
#########################################################################################################
ATE_boot <- NULL
for (i in 1:100) {
  #Sample with replacement
  id <- sample(1:nrow(dataCOHORT_CC), size = nrow(dataCOHORT_CC), replace = TRUE)
  data_boot <- dataCOHORT_CC[id,]
  
  ATE_boot[i] <- gc_ATE(data_boot,
                        formula = formula(VAP_num ~ GROUP_num + AGE + BMI + SEX +
                                            GLASGOW + INJURY + PAO2FIO2),
                        group = "GROUP_num")
}

#########################################################################################################
### 6) c)
#########################################################################################################
mean(ATE_boot)
quantile(ATE_boot, probs = 0.975, na.rm = TRUE)
quantile(ATE_boot, probs = 0.025, na.rm = TRUE)
#The ATE is statistically significant since the 95% confidence interval does not include zero

#########################################################################################################

### Part 3
#########################################################################################################
### 7) a) and b)
#########################################################################################################
.f1 <- formula(VAP_num ~ GROUP_num + AGE + BMI + SEX +
                GLASGOW + INJURY + PAO2FIO2)

gc_bin_all <- gc_binary(formula = .f1, data = dataCOHORT_CC, 
                        group = "GROUP_num", model = "all",
                        boot.number = 100, effect = "ATE",
                        progress = TRUE, seed = 5186)

#########################################################################################################
### 8) a) and b) and c)
#########################################################################################################
gc_bin_all
summary(gc_bin_all, ci.type = "perc", ci.level = 0.95)

cat("ATE manual =", round(mean(ATE_boot), 3),
    " (95% CI:", 
    round(quantile(ATE_boot, probs = 0.025, na.rm = TRUE), 3), "to",
    round(quantile(ATE_boot, probs = 0.975, na.rm = TRUE), 3), ")")
summary(gc_bin_all, ci.type = "perc", ci.level = 0.95, unadjusted = FALSE)

#########################################################################################################
### 9) a) and b)
#########################################################################################################
data_LEUKO <- subset(dataCOHORT_CC, LEUKO == ">=20000")
dim(data_LEUKO) #n=319 patients

.f2 <- VAP_num ~ GROUP_num * DIABETES + PAO2FIO2 + GLASGOW + TIME_INTUBATION
.f3 <- VAP_num ~ GROUP_num + AGE + BMI + GLASGOW + INJURY + PAO2FIO2

gc_bin_2 <- gc_binary(formula = .f2, data = data_LEUKO,
                       group = "GROUP_num", model = "all",
                       boot.number = 100, effect = "ATE",
                       progress = FALSE, seed = 5186)
#Patients with missing data are automatically removed from the analysis

gc_bin_3 <- gc_binary(formula = .f3, data = data_LEUKO,
                       group = "GROUP_num", model = "all",
                       boot.number = 100, effect = "ATE",
                       progress = FALSE, seed = 5186)

#########################################################################################################
### 9) c)
#########################################################################################################
plot(gc_bin_2, method = "calibration", 
     main = "Calibration Plot for .f2 (LEUKO >= 20000)")
plot(gc_bin_3, method = "calibration", 
     main = "Calibration Plot for .f3 (LEUKO >= 20000)")
#The first model shows poor calibration with predicted probabilities deviating from
# the diagonal (over and underestimation of risk).
#The second model is well-calibrated, as predictions overlap the diagonal line.

### BONUS
#########################################################################################################
### 10)
#########################################################################################################
dataCOHORT$DEATH_num <- ifelse(dataCOHORT$DEATH == "Yes", 1, 0)

.fsurv <- formula(Surv(TIME_DEATH, DEATH_num) ~ GROUP_num + AGE + BMI + SEX +
                                                GLASGOW + INJURY + PAO2FIO2)

gc_surv <- gc_times(formula = .fsurv, model = "all", data = dataCOHORT, 
                    group = "GROUP_num", boot.number = 100, effect = "ATE", 
                    pro.time = 60, #Time in days to estimate marginal survival
                    seed = 5186)
summary(gc_surv) #RMST = Restricted mean survival time up to pro.time (here = 60 days)
plot(gc_surv)

#########################################################################################################
### 11)
#########################################################################################################
.f4 <- VAP_num ~ GROUP_num * (AGE + SEX + BMI + DIABETES + ALCOOL +
                                SMOKING + INJURY + GLASGOW + PAO2FIO2 +
                                LEUKO + TIME_INTUBATION)

gc_all <- gc_binary(formula = .f4, model = "all", data = dataCOHORT,
                   group = "GROUP_num", boot.number = 100, effect = "ATE",
                   seed = 5186) 

gc_mi <- gc_binary(formula = .f4, model = "all", data = dataCOHORT,
                   group = "GROUP_num", boot.number = 100, effect = "ATE",
                   boot.mi = TRUE, #Apply Multiple Imputation before GC
                   m = 5, #Number of imputed datasets
                   seed = 5186) 

summary(gc_all, ci.type = "perc", unadjusted = FALSE)
summary(gc_mi, ci.type = "perc", unadjusted = FALSE)
#Note that the population is not the same, n=400 patients were excluded in the first model
#Excluding patients can introduce selection bias, while imputing their values adds variability.

plot(gc_mi, method = "calibration") #One plot per imputed dataset
plot(gc_mi, method = "calibration", smooth = TRUE) #Plots smoothed calibration curve

#########################################################################################################
### 11)
#########################################################################################################
#Using penalized regression (Lasso) for variable selection
gc_lasso <- gc_binary(formula = .f4, #complete formula, could even add splines with bs()
                      data = dataCOHORT, group = "GROUP_num", model = "lasso",
                      boot.number = 100, effect = "ATE", progress = TRUE, 
                      cv = 10, #number off cross-validation to use to estimate the hyperparameter
                      boot.tune = FALSE, #whether to re-estimate the hyperparamter in each bootstrap iteration
                      seed = 5186)

summary(gc_all, ci.type = "perc", unadjusted = FALSE)
summary(gc_lasso, ci.type = "perc", unadjusted = FALSE)

gc_lasso$tuning.parameters #the hyperparameter lambda obtained by cross-validation on the total population

#The hyperparameters can be directly given or a grid can be given to search on with param.tune :
.tune.lambda <- 10^seq(-4, 1, length.out = 100)
gc_lasso2 <- gc_binary(formula = .f4, #complete formula, could even add splines with bs()
                      data = dataCOHORT, group = "GROUP_num", model = "lasso",
                      boot.number = 100, effect = "ATE", progress = TRUE, 
                      cv = 10, #number off cross-validation to use to estimate the hyperparameter
                      boot.tune = FALSE, #whether to re-estimate the hyperparamter in each bootstrap iteration
                      param.tune = .tune.lambda, #search grid for lambda
                      seed = 5186)

gc_lasso2$tuning.parameters

#########################################################################################################