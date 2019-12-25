rm(list=ls())

# library for dataframe manipulation
library(fGarch)
library("MASS")

## Parameters

# Portfolio directory in Data folder
pf_path <- 'Data/Portfolio1/'

# Number of days ahead the VaR is calculated
VaR_days <- 5

# VaR alpha
VaR_alpha <- 0.05

# Monte Carlo sim
MC_n <- 10000

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
# gfit01  <- garchFit(formula = ~ garch(1, 1), data=pf_log_mean0, cond.dist=GARCHcondDist)
# tgarch
gfit01 <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=pf_log_mean0, cond.dist=GARCHcondDist)

# Functions

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


mu <- coef(gfit01)["mu"]
omega <- coef(gfit01)["omega"]
alpha <- coef(gfit01)["alpha1"]
beta <- coef(gfit01)["beta1"]

#tgarch
#gamma <- coef(gfit01)["gamma1"]

h.t <- gfit01@h.t
Z <- gfit01@residuals/gfit01@sigma.t 

# fit student t dist
Z_fit <- fitdistr(Z,"t")
Z_mu <- Z_fit$estimate[["m"]]
Z_s <- Z_fit$estimate[["s"]]
Z_df <- Z_fit$estimate[["df"]]


# hist(pf_ret,breaks=80,main='Cree returns',freq=F,density=30,col='cyan',ylim=c(0,20),xlim=c(-0.2,0.3))
# lines(seq(-0.5,0.5,0.001),dt((seq(-0.5,0.5,0.001)-Z_mu)/Z_s,Z_df)/Z_s,col='blue',lwd=2)
# legend('topright',c('Fitted normal', 'Fitted t-distribution'),col=c('red', 'blue'),lwd=2)


#MC_Z <- array(sample(Z,pf_days*VaR_days*MC_n,replace = TRUE), dim=c(pf_days,VaR_days,MC_n))
MC_Z <- array(rt(pf_days*VaR_days*MC_n,df=Z_df), dim=c(pf_days,VaR_days,MC_n))
MC_Z <- MC_Z * Z_s + Z_mu


MC_h <- array(dim=c(pf_days,VaR_days,MC_n))

MC_h[,1,] <- matrix(h.t,nrow=pf_days,ncol=MC_n)

for(i in 2:VaR_days){
  # garch
  # MC_h[,i,] <- GARCH_ht_function(omega,alpha,beta,MC_h[,i-1,],MC_Z[,i-1,])
  # tgarch
  MC_h[,i,] <- TGARCH_ht_function(omega,alpha,beta,gamma,MC_h[,i-1,],MC_Z[,i-1,])
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

# Kupiec test
library(Rmpfr)

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

