rm(list=ls())

# library for dataframe manipulation
library(dplyr)
library(fGarch)

## Parameters

# Portfolio directory in Data folder
pf_path <- 'Data/Portfolio1/'

# Number of days ahead the VaR is calculated
VaR_days <- 5

# VaR alpha
VaR_alpha <- 0.05

# Monte Carlo sim
MC_n <- 1000

# Conditional distribution GARCH: Students t distribution
GARCHcondDist <- "std"

## Run data handling file 
source('DataHandling.R')

# Precalculate some numbers for higher calculation performance
pf_days <- length(pf_log)

## Get mean returns
pf_log_mean <- mean(pf_log)

# Deduct the mean from every stock, for standardized news process in GARCH
pf_log_mean0 <- pf_log - pf_log_mean


plot(pf_log)

# garch
gfit01  <- garchFit(formula = ~ garch(1, 1), data=pf_log_mean0, cond.dist=GARCHcondDist)
# tgarch
#gfit01 <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=pf_log_mean0, cond.dist=GARCHcondDist)

# Functions

# GARCH

GARCH_ht_function <- function(omega,alpha,beta,hPrevious,zPrevious){
  epsilon <- sqrt(hPrevious)*zPrevious
  h <- omega + alpha*epsilon^2 + beta*hPrevious
  return(h)
}


# TGARCH 

TGARCH_ht_function <- function(omega,alpha,beta,gamma,hPrevious,zPrevious){
  epsilon <- sqrt(hPrevious)*zPrevious
  h <- omega + alpha*epsilon^2 + gamma*pmin(epsilon,0)^2 +beta*hPrevious
  return(h)
}


mu <- coef(gfit01)[1]
omega <- coef(gfit01)[2]
alpha <- coef(gfit01)[3]

#garch
beta <- coef(gfit01)[4]
#tgarch
#gamma <- coef(gfit01)[4]
#beta <- coef(gfit01)[5]

h.t <- gfit01@h.t
Z <- gfit01@residuals/gfit01@sigma.t 


MC_Z <- array(sample(Z,pf_days*VaR_days*MC_n,replace = TRUE), dim=c(pf_days,VaR_days,MC_n))

MC_h <- array(dim=c(pf_days,VaR_days,MC_n))

MC_h[,1,] <- matrix(h.t,nrow=pf_days,ncol=MC_n)

for(i in 2:VaR_days){
  # garch
  MC_h[,i,] <- GARCH_ht_function(omega,alpha,beta,MC_h[,i-1,],MC_Z[,i-1,])
  # tgarch
  #MC_h[,i,] <- TGARCH_ht_function(omega,alpha,beta,gamma,MC_h[,i-1,],MC_Z[,i-1,])
}

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

