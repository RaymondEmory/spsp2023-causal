rm(list=ls())
set.seed(202302)

# helper expit function (also termed inverse logit)
expit <- function(x) exp(x)/(1+exp(x))

# function to generate a single observed dataset
OneD <- function(n=1000) {
  # baseline covariates L
  L1 <- rexp(n=n, rate=2)
  L2 <- rpois(n=n, lambda=0.5)
  L <- cbind(L1,L2)

  # true propensity score
  PS <- expit(rowSums(L))
  # observed treatment assignment or selection based on propensity score
  TRT <- rbinom(n=n,size=1, prob=PS)
  
  # potential outcomes under no treatment
  Y0 <- rnorm(n=n, mean=2*rowSums(L), sd=1)
  # no individual causal effect
  Y1 <- Y0
  Y <- TRT*Y1 + (1-TRT)*Y0
  return(data.frame(L,TRT,Y))
}

write.csv(OneD(),file=paste0("spsp2023-causal-lab-Data",1,".csv"),
          row.names=FALSE)
