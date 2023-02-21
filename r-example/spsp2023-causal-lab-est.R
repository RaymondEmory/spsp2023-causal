rm(list=ls())

# download data from GitHub
D <- read.csv("https://raw.githubusercontent.com/wwloh/spsp2023-causal/main/spsp2023-causal-lab-Data1.csv")
head(D)

# vector to save estimates
est <- NULL

# no confounding adjustment ###################################################
fit.noL <- lm(Y~TRT, data=D)
summary(fit.noL)
est["noL"] <- coef(fit.noL)["TRT"]

# outcome regression ##########################################################
fit.regL <- lm(Y~TRT+L1+L2, data=D)
summary(fit.regL)
est["regL"] <- coef(fit.regL)["TRT"]

round(est,2)

# propensity score methods ####################################################
# fit propensity score model using logistic regression
fit.ps <- glm(TRT~L1+L2, data=D, family=binomial('logit'))
# calculate predicted propensity scores
pred.ps <- predict.glm(fit.ps, type="response")

# add propensity scores as additional covariate
D[,"PS"] <- pred.ps
fit.addPS <- lm(Y~TRT+PS, data=D)
est["addPS"] <- coef(fit.addPS)["TRT"]

# use inverse probability weights
D[,"W"] <- D[,"TRT"]*(1/pred.ps) + (1-D[,"TRT"])*(1/(1-pred.ps))
fit.regIPW <- lm(Y~TRT, data=D, weights=W)
est["IPW"] <- coef(fit.regIPW)["TRT"]

round(est,2)

# bootstrap SEs and CIs #######################################################
## function to generate a single bootstrap sample and calculate estimates
OneBoot <- function(DATA) {
  # resample observations
  boot.id <- sort(unique(sample(nrow(DATA),nrow(DATA),replace = TRUE)))
  bootD <- DATA[boot.id,]
  # estimates of effect
  bootest <- NULL
  fit.noL <- lm(Y~TRT, data=bootD)
  bootest["noL"] <- coef(fit.noL)["TRT"]
  fit.regL <- lm(Y~TRT+L1+L2, data=bootD)
  bootest["regL"] <- coef(fit.regL)["TRT"]
  fit.ps <- glm(TRT~L1+L2, data=bootD, family=binomial('logit'))
  pred.ps <- predict.glm(fit.ps, type="response")
  bootD[,"PS"] <- pred.ps
  fit.addPS <- lm(Y~TRT+PS, data=bootD)
  bootest["addPS"] <- coef(fit.addPS)["TRT"]
  bootD[,"W"] <- bootD[,"TRT"]*(1/pred.ps) + (1-bootD[,"TRT"])*(1/(1-pred.ps))
  fit.regIPW <- lm(Y~TRT, data=bootD, weights=W)
  bootest["IPW"] <- coef(fit.regIPW)["TRT"]
  return(bootest)
}
est.boot <- replicate(n=500,OneBoot(DATA=D))
# bootstrap CIs
round(apply(est.boot,1,quantile, probs=c(0.025, 0.975)),2)
# bootstrap SEs
round(apply(est.boot,1,sd),2)


# visual diagnostics: overlap in propensity scores between the treatment groups
boxplot(pred.ps~D[,"TRT"], horizontal=TRUE,
        main="Propensity score distributions", 
        ylab="Propensity score estimate", xlab="Treatment group")
