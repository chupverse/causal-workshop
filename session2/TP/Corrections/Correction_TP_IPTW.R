######################################################
# ecrit par : Lisa Durocher 
# version : 1
# description : TP IPTW formation causalite
######################################################


library(lmtest)
library(sandwich)
library(RISCA)
library(jtools)
library(glmnet)
library(mice)
library(cobalt)

#### 1. Load the database dataCOHORT #### 
load("dataCOHORT.RData")
dataCOHORT.complet <- dataCOHORT

#### 2. Looking at the descriptive table, check the positivity assumption #### 

dataCOHORT$GROUP <- as.character(dataCOHORT$GROUP)
dataCOHORT$GROUP <- ifelse(dataCOHORT$GROUP == "Untreated", 0, 1)

# The port algorithm from the RISCA package allows to identify potential positivity violations.
port(group = "GROUP", cov.quanti = c("BMI","AGE"), 
     cov.quali = c("SEX","DIABETES","ALCOOL","SMOKING","INJURY","GLASGOW","PAO2FIO2","LEUKO"),
     data = dataCOHORT, alpha = 0.03, beta = 0.02, gamma = 2)

# We choose to exclude patients under the age of 30.
dataCOHORT <- dataCOHORT[-which(dataCOHORT$AGE < 30),]

#### 3. Compare the outcome (VAP) observed in each group  ####

length(dataCOHORT$VAP[which(dataCOHORT$GROUP == 0 & dataCOHORT$VAP == "Yes")]) / length(dataCOHORT$VAP[which(dataCOHORT$GROUP == 0)]) *100
length(dataCOHORT$VAP[which(dataCOHORT$GROUP == 1 & dataCOHORT$VAP == "Yes")]) / length(dataCOHORT$VAP[which(dataCOHORT$GROUP == 1)]) *100

chisq.test(dataCOHORT$VAP, dataCOHORT$GROUP)

#### 5. Construct the exposure model using multivariate logistic regression ####

# We have to work in complete-case
dataCOHORT <- dataCOHORT[-which(is.na(dataCOHORT$AGE) | is.na(dataCOHORT$BMI) | is.na(dataCOHORT$GLASGOW) | 
                                  is.na(dataCOHORT$PAO2FIO2) | is.na(dataCOHORT$SEX) | is.na(dataCOHORT$INJURY)),]

# We choose we will choose to adjust on the age, the BMI, the gender, the Glasgow Coma Scale, the injury and the PAO2/FIO2.
ps <- glm(GROUP ~ AGE + BMI + GLASGOW + PAO2FIO2 + SEX + INJURY,
          data = dataCOHORT, family = binomial(link = "logit"))

#### 6. For each subject, calculate his propensity score, i.e. the probability of being treated according to the individual characteristics #### 

Pr1 = fitted(ps)
dataCOHORT$Pr1 <- Pr1

#### 7. Check the positivity assumption by plotting the distribution of the propensity score in the two separate groups. ####

H1 <- hist(x = Pr1[dataCOHORT$GROUP == 1], breaks = 40, plot = F)
H0 <- hist(x = Pr1[dataCOHORT$GROUP == 0], breaks = 40, plot = F)
hist(x = Pr1[dataCOHORT$GROUP == 1], xlim = c(0,1), ylim = c(max(max(H1$counts),max(H0$counts)),
                                                             -max(max(H1$counts),max(H0$counts))),
     breaks = 40, col = "gray60", .5, axes = FALSE, ylab = "", xlab = "", main = "")
par(new = T)
hist(x = Pr1[dataCOHORT$GROUP == 0], xlim = c(0,1), ylim = c(-max(max(H1$counts),max(H0$counts)),
                                                             max(max(H1$counts),max(H0$counts))),
     breaks = 40, col = "gray86", .5, ylab = "Frequency", xlab = "Propensity score", main = "")
text(0.95, -max(H1$counts),"Ceftriaxone", font = 2, pos = 2)
text(0.95, max(H0$counts),"Untreated", font = 2, pos = 2)

#### 8. Calculate weights ####

# Probability of being treated. 
Pr0 = glm(GROUP~1,data=dataCOHORT, family = binomial(link = "logit"))$fitted.values

# Weights for ATE 
dataCOHORT$W <- (dataCOHORT[,"GROUP"]==1) *Pr0/Pr1 + (dataCOHORT[,"GROUP"] == 0) * (1-Pr0)/(1-Pr1) 

#### 9. Check the correct specifications for your model by calculating the standardized mean difference3 (SMD) of the covariates AGE and SEX in the counterfactual population ####

smd <- function(data,var,TRT, WEIGHT){
  nbvar <- length(var)
  group0 <- c()
  group1 <- c()
  SMD <- c()
  varnames <- c()
  for (i in 1:nbvar) {
    if (is.numeric(data[,var[i]])){
      z1 <- weighted.mean(data[,var[i]][which(data[,TRT] == 1)], data[,WEIGHT][which(data[,TRT] == 1)], na.rm = T)
      z0 <- weighted.mean(data[,var[i]][which(data[,TRT] == 0)], data[,WEIGHT][which(data[,TRT] == 0)], na.rm = T)
      s1 <- wtd.sd(data[,var[i]][which(data[,TRT] == 1)], data[,WEIGHT][which(data[,TRT] == 1)])**2
      s0 <- wtd.sd(data[,var[i]][which(data[,TRT] == 0)], data[,WEIGHT][which(data[,TRT] == 0)])**2
      group1 <- c(group1, format(round(z1, 2), nsmall = 2))
      group0 <- c(group0, format(round(z0, 2), nsmall = 2))
      SMD <- c(SMD,format(round((abs(z1-z0)/ sqrt((s1+s0)/2)) *100,2),nsmall=2))
      varnames <- c(varnames,var[i])
    }
    if(is.factor(data[,var[i]])){
      if (length(levels(data[,var[i]])) == 2) {
        tot1 <- sum(dataCOHORT[,WEIGHT][which(dataCOHORT[,TRT] == 1)])
        tot0 <- sum(dataCOHORT[,WEIGHT][which(dataCOHORT[,TRT] == 0)])
        p1 <- sum(dataCOHORT[,WEIGHT][which(dataCOHORT[,TRT] == 1 &
                                              dataCOHORT[,var[i]] == levels(dataCOHORT[,var[i]])[1])])/tot1 
        p0 <- sum(dataCOHORT[,WEIGHT][which(dataCOHORT[,TRT] == 0 &
                                              dataCOHORT[,var[i]] == levels(dataCOHORT[,var[i]])[1])])/tot0
        group1 <- c(group1, format(round(p1*100, 2), nsmall = 2))
        group0 <- c(group0, format(round(p0*100, 2), nsmall = 2))
        SMD <- c(SMD,format(round((abs(p1-p0)/sqrt((p1*(1-p1)+p0*(1-p0))/2)) *100,2),nsmall=2))
        varnames <- c(varnames, paste(var[i],":",levels(dataCOHORT[,var[i]])[1]))
      
      }
      
      if (length(levels(data[,var[i]])) > 2) {
        tot1 <- sum(dataCOHORT[,WEIGHT][which(dataCOHORT[,TRT] == 1)])
        tot0 <- sum(dataCOHORT[,WEIGHT][which(dataCOHORT[,TRT] == 0)])
        for (j in 1:length(levels(data[,var[i]]))) {
          p1 <- sum(dataCOHORT[,WEIGHT][which(dataCOHORT[,TRT] == 1 &
                                                dataCOHORT[,var[i]] == levels(dataCOHORT[,var[i]])[j])])/tot1
          p0 <- sum(dataCOHORT[,WEIGHT][which(dataCOHORT[,TRT] == 0 &
                                                dataCOHORT[,var[i]] == levels(dataCOHORT[,var[i]])[j])])/tot0
          group1 <- c(group1, format(round(p1*100, 2), nsmall = 2))
          group0 <- c(group0, format(round(p0*100, 2), nsmall = 2))
          SMD <- c(SMD,format(round((abs(p1-p0)/sqrt((p1*(1-p1)+p0*(1-p0))/2)) *100,2),nsmall=2))
          varnames <- c(varnames, paste(var[i],":",levels(dataCOHORT[,var[i]])[j]))
        }
      }
    }
  }
  data.frame("variable"=varnames, "Group 1" = group1, "Group 0" = group0, "SMD" = SMD)
}
# SMD for covariate we choose to adjust on
smd(data = dataCOHORT, var =c("AGE","BMI","GLASGOW","PAO2FIO2","SEX","INJURY"), TRT = "GROUP", WEIGHT = "W")
# SMD for the other covariate
smd(data = dataCOHORT, var =c("DIABETES","ALCOOL","SMOKING","LEUKO"), TRT = "GROUP", WEIGHT = "W")


# love plot 
dataCOHORT$Pr1 <- Pr1

balIPTW <- bal.tab(x = GROUP~.,
                   data = dataCOHORT[,-which(colnames(dataCOHORT) %in% c("VAP","TIME_DEATH","DEATH"))],
                   distance = "Pr1",
                   weights = "W",
                   method = "weighting",
                   disp.v.ratio = TRUE,
                   un = TRUE,
                   s.d.denom = "treated")

love.plot(x = balIPTW,
          stat = "mean.diffs",
          abs = TRUE,
          var.order = "unadjusted",
          threshold = 0.1,
          shape=23)


#### 10. Calculate the average treatment effect on the entire population (ATE) by computing the difference of proportion of VAP between treated and untreated ####

glm.ipw <- glm(VAP ~ GROUP, weights = W, family = quasibinomial, data=dataCOHORT)
summary(glm.ipw)$coefficients

# The coeftest function from the lmtest package is used to obtain robust standard errors.

# estimated rate of VAP in the group untreated with 95% confidence interval
glm.rob <- coeftest(glm.ipw, vcov=sandwich)
P0 <- exp(glm.rob[1,1])/(1+exp(glm.rob[1,1]))
P0_inf <- exp(glm.rob[1,1] - 1.96*glm.rob[1,2])/(1+exp(glm.rob[1,1] - 1.96*glm.rob[1,2]))
P0_sup <- exp(glm.rob[1,1] + 1.96*glm.rob[1,2])/(1+exp(glm.rob[1,1] + 1.96*glm.rob[1,2]))
sum(dataCOHORT$W[which(dataCOHORT$GROUP == 0 & dataCOHORT$VAP == "Yes")]) / sum(dataCOHORT$W[which(dataCOHORT$GROUP == 0)]) *100

# estimated rate of VAP in the group treated with 95% confidence interval
dataCOHORT$GROUP <- abs(dataCOHORT$GROUP - 1)
glm.ipw <- glm(VAP ~ GROUP, weights = W, family = quasibinomial, data=dataCOHORT)
glm.rob <- coeftest(glm.ipw, vcov=sandwich)
P1 <- exp(glm.rob[1,1])/(1+exp(glm.rob[1,1]))
P1_inf <- exp(glm.rob[1,1] - 1.96*glm.rob[1,2])/(1+exp(glm.rob[1,1] - 1.96*glm.rob[1,2]))
P1_sup <- exp(glm.rob[1,1] + 1.96*glm.rob[1,2])/(1+exp(glm.rob[1,1] + 1.96*glm.rob[1,2]))
sum(dataCOHORT$W[which(dataCOHORT$GROUP == 1 & dataCOHORT$VAP == "Yes")]) / sum(dataCOHORT$W[which(dataCOHORT$GROUP == 1)]) *100

# estimated difference and pvalue
dataCOHORT$GROUP <- abs(dataCOHORT$GROUP - 1)
glm.ipw <- glm(VAP ~ GROUP, weights = W, family = quasibinomial, data=dataCOHORT)
glm.rob <- coeftest(glm.ipw, vcov=sandwich)
diff <- P0 - P1
pvalue <- glm.rob[2,4]

data.frame("P0" = paste0(round(P0*100,2)," [95% CI: ", round(P0_inf*100,2), "-", round(P0_sup*100,2), "]"),
           "P1" = paste0(round(P1*100,2)," [95% CI: ", round(P1_inf*100,2), "-", round(P1_sup*100,2), "]"),
           "difference" = round(diff*100,2),
           "pvalue" = pvalue)

#### 11. Calculate the marginal effect of the group on the death ####

dataCOHORT$failures <- ifelse(dataCOHORT$DEATH == "No", 0, 1)

res.akm <-ipw.survival(times=dataCOHORT$TIME_DEATH, failures=dataCOHORT$failures,
                       variable=dataCOHORT$GROUP, weights=dataCOHORT$W)
plot(res.akm, ylab="Confounder-adjusted survival",
     xlab="Time (days)", col=c(1,2), grid.lty=1)

coxph(Surv(TIME_DEATH, failures) ~ GROUP, weights = W, data = dataCOHORT)

#################################################################################
#### Question 12, 13, 14 ####
#################################################################################

#### MICE ####

dataCOHORT.complet$GROUP <- as.character(dataCOHORT.complet$GROUP)
dataCOHORT.complet$GROUP <- ifelse(dataCOHORT.complet$GROUP == "Untreated", 0, 1)

sapply(dataCOHORT.complet, function(x) sum(is.na(x))) *100/nrow(dataCOHORT.complet)

dm = mice(dataCOHORT.complet,maxit = 0)
pred <- dm$predictorMatrix
meth <- dm$method
# The group and the endpoints are not used to impute covariates.
pred[, c("GROUP","VAP","TIME_DEATH","DEATH","TIME_INTUBATION")]=0

# We create 5 imputed databases. 
dm = mice(dataCOHORT.complet, m = 5, print = FALSE, maxit = 20, predictorMatrix = pred, method = meth, vis =  "monotone")
summary(dm)
sapply(complete(dm,1), function(x) sum(is.na(x))) *100/nrow(complete(dm,1))

#### Transform categorical variable as indicator variable.  ####

categorial <- function(data, without = NULL){
  
  if(is.null(without)){
    name <- colnames(data)
  }else{
    name <- colnames(data[,-which(colnames(data) %in% without)])
  }
  
  nb <- length(name)
  var<- NULL
  for (col in 1:nb) {
    i <- name[col]
    if (is.factor(data[,i])){
      if (length(levels(data[,i]))>2){
        
        for (j in 2:length(levels(data[,i]))) {
          data$a <- as.factor(ifelse(data[,i] == sort(levels(data[,i]))[j], 1, 0))
          colnames(data)[ncol(data)] <- gsub(" ","",paste0(i,"_",sort(levels(data[,i]))[j]))
          var <- c(var, i)
        }
      
      } 
    }
  }
  data <- data[,-which(colnames(data) %in% var)]
  return(data)
}

categorial(dataCOHORT.complet)


#### Lasso ####

data <- complete(dm,1)
data <- categorial(data)
data <- data[,c("AGE","SEX","BMI","DIABETES","ALCOOL","SMOKING","LEUKO","INJURY_Ischemicstroke",
                "INJURY_Subarachnoidhemorrhage","INJURY_Traumabraininjury","GLASGOW_3","GLASGOW_4-8",
                "PAO2FIO2_>=200","PAO2FIO2_100-199","GROUP","VAP")]

data.b<- data

# The tuning parameter is estimated by minimizing the 20-fold cross-validation mean square errors
set.seed(1234)
trainX <- model.matrix(VAP~., data = data.b[,-which(colnames(data.b)%in% c("GROUP","TIME_INTUBATION","TIME_DEATH","DEATH"))])[,-1]

lasso.mod = cv.glmnet(trainX,data.b[,"VAP"],alpha=1, family = binomial(link = "logit") ,nfolds = 20)
lambda.min <- lasso.mod$lambda.min

# Multivariate linear regression with Lasso penalization
lasso.mod2 = glmnet(trainX,data.b[,"VAP"],lambda = lambda.min, alpha=1, family = binomial(link = "logit"))
lasso.mod2$beta

# all covariates are selected as risk factors

#### IPTW ####

ps <- glm(GROUP ~ .,
          data = data.b[,-which(colnames(data.b) == "VAP")], family = binomial(link = "logit"))

Pr1 = fitted(ps)

Pr0 = glm(GROUP~1,data=data.b, family = binomial(link = "logit"))$fitted.values

data.b$W <- (data.b[,"GROUP"]==1) *Pr0/Pr1 + (data.b[,"GROUP"] == 0) * (1-Pr0)/(1-Pr1) 


glm.ipw <- glm(VAP ~ GROUP, weights = W, family = quasibinomial, data=data.b)
summary(glm.ipw)


p0 <- sum(data.b$W[which(data.b$GROUP == 0 & data.b$VAP == "Yes")]) / sum(data.b$W[which(data.b$GROUP == 0)]) *100
p1 <- sum(data.b$W[which(data.b$GROUP == 1 & data.b$VAP == "Yes")]) / sum(data.b$W[which(data.b$GROUP == 1)]) *100


variablePS <- as.data.frame(matrix(data=NA, 1, 16))
colnames(variablePS) <- c(colnames(data.b[1:14]),"P0","P1")
variablePS[,1:14] <- lasso.mod2$beta
variablePS[,"P0"] <- p0
variablePS[,"P1"] <- p1


#### Function to bootstrap the previous step. ####

MI.boot <- function(data, GROUP, OUTCOME){
  # creation of the bootstrap sample
  bootstrap <- sample(1:nrow(data), size=nrow(data), replace=TRUE)
  data.b<- data[bootstrap,] 
  
  # estimation of the tunning parameter
  trainX <- model.matrix(as.formula(paste0(OUTCOME,"~.")), data = data.b[,-which(colnames(data.b)%in% c(GROUP,"TIME_INTUBATION","TIME_DEATH","DEATH"))])[,-1]
  
  lasso.mod = cv.glmnet(trainX,data.b[,OUTCOME],alpha=1, family = binomial(link = "logit") ,nfolds = 20)
  lambda.min <- lasso.mod$lambda.min
  
  # Multivariate linear regression with Lasso penalization
  lasso.mod2 = glmnet(trainX,data.b[,OUTCOME],lambda = lambda.min, alpha=1, family = binomial(link = "logit"))

  #  selection of the risk factors to adjust on
  if (length(which(lasso.mod2$beta != 0)) >0 ) {
    if (length(which(lasso.mod2$beta == 0)) !=0 ){
      d <- data.b[,-which(lasso.mod2$beta == 0)]
    }else { 
      d <- data.b
    }
  }
  
  # IPTW 
  ps <- glm(as.formula(paste0(GROUP, "~ .")),
            data = d[,-which(colnames(d) == OUTCOME)], family = binomial(link = "logit"))
  Pr1 = fitted(ps)
  Pr0 = glm(as.formula(paste0(GROUP,"~1")), data=d, family = binomial(link = "logit"))$fitted.values
  d$W <- (d[,GROUP]==1) *Pr0/Pr1 + (d[,GROUP] == 0) * (1-Pr0)/(1-Pr1) 
  glm.ipw <- glm(as.formula(paste0(OUTCOME, "~", GROUP)), weights = W, family = quasibinomial, data=d)
  summary(glm.ipw)
  p0 <- sum(d$W[which(d[,GROUP] == 0 & d[,OUTCOME] == "Yes")]) / sum(d$W[which(d[,GROUP] == 0)]) *100
  p1 <- sum(d$W[which(d[,GROUP] == 1 & d[,OUTCOME] == "Yes")]) / sum(d$W[which(d[,GROUP] == 1)]) *100
  
  # we save the results on a data frame
  variablePS <- as.data.frame(matrix(data=NA, 1, 16))
  colnames(variablePS) <- c(colnames(data.b[1:14]),"P0","P1")
  variablePS[,1:14] <- lasso.mod2$beta
  variablePS[,"P0"] <- p0
  variablePS[,"P1"] <- p1
  
  return(variablePS)
}

IPTW.result <- as.data.frame(matrix(data=NA, 1, 16))
colnames(IPTW.result) <- c("AGE","SEX","BMI","DIABETES","ALCOOL","SMOKING","LEUKO","INJURY_Ischemicstroke","INJURY_Subarachnoidhemorrhage",
                           "INJURY_Traumabraininjury","GLASGOW_3","GLASGOW_4-8","PAO2FIO2_>=200","PAO2FIO2_100-199","P0","P1")
for (i in 1:5) {
  
  data <- complete(dm,i)
  data <- categorial(data)
  data <- data[,c("AGE","SEX","BMI","DIABETES","ALCOOL","SMOKING","LEUKO","INJURY_Ischemicstroke",
                  "INJURY_Subarachnoidhemorrhage","INJURY_Traumabraininjury","GLASGOW_3","GLASGOW_4-8",
                  "PAO2FIO2_>=200","PAO2FIO2_100-199","GROUP","VAP")]
  
  for (j in 1:100) {
    a <- MI.boot(data = data, GROUP = "GROUP", OUTCOME = "VAP")
    IPTW.result <- rbind(IPTW.result, a)
  }
  
  
}

# 95% confidence interval
round(quantile(IPTW.result$P0, probs = c(0.025, 0.5, 0.975), na.rm = T),4)
round(quantile(IPTW.result$P1, probs = c(0.025, 0.5,0.975), na.rm = T),4)
IPTW.result$DIFF <- IPTW.result$P1 - IPTW.result$P0
round(quantile(data_var$DIFF, probs = c(0.025, 0.975), na.rm = T),4)


