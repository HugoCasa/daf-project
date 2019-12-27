# rm(list=ls())
# 
# # library for dataframe manipulation
# library(fGarch)
# library("MASS")
# 
# ## Parameters
# 
# # Portfolio directory in Data folder
# # portfolio 1
# #pf_path <- 'Data/Portfolio1/'
# # portfolio 2
# pf_path <- 'Data/Portfolio2/'
# 
# # Number of days ahead the VaR is calculated
# VaR_days <- 5
# 
# # VaR alpha
# VaR_alpha <- 0.05
# 
# # Monte Carlo sim
# MC_n <- 1000
# 
# # Conditional distribution GARCH: Students t distribution
# GARCHcondDist <- "std"
# 
# # GARCH model
# #GARCH_model <- 'GARCH'
# GARCH_model <- 'TGARCH'
# 
# # sample or fit t distribution on standardized returns
# #ret_method <- "sample"
# ret_method <- "fit"
# 
# ## Run data handling file 
# source('DataHandling.R')


# Record start time
time_start <- Sys.time()


# Precalculate some numbers for higher calculation performance
pf_days <- length(pf_log)

## Get mean returns
pf_log_mean <- mean(pf_log)

# Deduct the mean from every stock, for standardized news process in GARCH
pf_log_mean0 <- pf_log - pf_log_mean


# fit garch model

if (GARCH_model == "GARCH") {
  gfit01  <- garchFit(formula = ~ garch(1, 1), data=pf_log_mean0, cond.dist=GARCHcondDist)
  
  mu <- coef(gfit01)[["mu"]]
  omega <- coef(gfit01)[["omega"]]
  alpha <- coef(gfit01)[["alpha1"]]
  beta <- coef(gfit01)[["beta1"]]
}
if (GARCH_model == "TGARCH") {
  gfit01 <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=pf_log_mean0, cond.dist=GARCHcondDist)
  
  mu <- coef(gfit01)[["mu"]]
  omega <- coef(gfit01)[["omega"]]
  alpha <- coef(gfit01)[["alpha1"]]
  beta <- coef(gfit01)[["beta1"]]
  
  gamma <- coef(gfit01)[["gamma1"]]
}

## Functions

# GARCH

GARCH_ht_function <- function(omega,alpha,beta,hPrevious,zPrevious){
  epsilon <- sqrt(hPrevious)*zPrevious
  h <- omega + alpha*epsilon^2 + beta*hPrevious
  return(h)
}

# TGARCH 

TGARCH_ht_function <- function(omega,alpha,beta,gamma,hPrevious,zPrevious){
  h <- omega + alpha*hPrevious*(abs(zPrevious)-gamma*zPrevious)^2 + beta*hPrevious
  return(h)
}

# conditional volatility
h.t <- gfit01@h.t
# standardized returns
Z <- gfit01@residuals/gfit01@sigma.t 

if (ret_method == "sample") {
  # sample from standardized residuals (3 dimensions: for each trading day, each nday VaR and each simulation per day)
  MC_Z <- array(sample(Z,pf_days*VaR_days*MC_n,replace = TRUE), dim=c(pf_days,VaR_days,MC_n))
}
if (ret_method == "fit") {
  
  # fit student t distribution
  Z_fit <- fitdistr(Z,"t")
  Z_mu <- Z_fit$estimate[["m"]]
  Z_s <- Z_fit$estimate[["s"]]
  Z_df <- Z_fit$estimate[["df"]]
  
  # plot fitting
  hist(Z,breaks=80,main='Standardised Residuals',freq=F,col='cyan')
  lines(seq(-10,10,0.01),dt((seq(-10,10,0.01)-Z_mu)/Z_s,Z_df)/Z_s,col='red',lwd=2)
  legend('bottomright',c('Fitted t-distribution'),col=c('red'), lwd=2)
  
  
  # sample from fitted student t distribution (3 dimensions: for each trading day, each nday VaR and each simulation per day)
  MC_Z <- array(rt(pf_days*VaR_days*MC_n,df=Z_df), dim=c(pf_days,VaR_days,MC_n))
  # Scale back all standardized residuals
  MC_Z <- MC_Z * Z_s + Z_mu
}

## Create array to store ht of simulations
MC_h <- array(dim=c(pf_days,VaR_days,MC_n))

# The variance of the first simulated day is the same for all
MC_h[,1,] <- matrix(h.t,nrow=pf_days,ncol=MC_n)

# Apply GARCH function
for(i in 2:VaR_days){
  if (GARCH_model == "GARCH") {
    MC_h[,i,] <- GARCH_ht_function(omega,alpha,beta,MC_h[,i-1,],MC_Z[,i-1,])
  }
  if (GARCH_model == "TGARCH") {
    MC_h[,i,] <- TGARCH_ht_function(omega,alpha,beta,gamma,MC_h[,i-1,],MC_Z[,i-1,])
  }
}

## Simulate returns
MC_log <- MC_Z*sqrt(MC_h)

# Add back the mean to every stock
MC_log <- MC_log + pf_log_mean

# Get cumulative returns of the n days
MC_log <- apply(MC_log,c(1,3),sum)


# Apply quantile function to get VaR
VaR <- apply(MC_log,1,quantile,VaR_alpha)


### Check model
hitSeq       <- pf_log_nday < VaR[1:length(pf_log_nday)]
numberOfHits <- sum(hitSeq)
exRatio      <- numberOfHits/length(pf_log_nday)

plot(pf_log_nday)
lines(VaR,col='green')


# Higher precision is needed, otherwise numerator and denumerator are treated as 0
N <- mpfr(length(pf_log),precBits= 128)
exRatio <- mpfr(exRatio,precBits = 128)
numberOfHits <- mpfr(numberOfHits,precBits = 128)
VaR_alpha <- mpfr(VaR_alpha,precBits = 128)
num <- (exRatio^numberOfHits)*(1-exRatio)^(N-numberOfHits)
den <- (VaR_alpha^numberOfHits )*(1-VaR_alpha)^(N-numberOfHits)
K   <- as.numeric(2*log(num/den))
VaR_alpha <- as.numeric(VaR_alpha)
p <- 0.99

if(K < qchisq(p,1)){
  print("VaR model is accurate at 99% level")
}else{
  print("VaR model is not accurate at 99% level")
}


#=================================================================================================
## Calculate some unconditional summary statistics for simulated portfolio returns

# For simulated n day portfolio returns
MC_log_mean <- mean(MC_log)
MC_log_std <- sqrt(var(as.vector(MC_log)))
MC_log_skewness <- skewness(as.vector(MC_log))
MC_log_kurtosis <- kurtosis(as.vector(MC_log),method = 'moment')


time_end <- Sys.time()
time_model <- time_end- time_start

printFile(1)
writeFile(1)
