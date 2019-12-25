library("copula")
library("fGarch")
library("MASS")

rm(list=ls())

## Parameters

# Portfolio directory in Data folder
pf_path <- 'Data/Portfolio1/'
#pf_path <- 'Data/Portfolio2/'

# Number of days ahead the VaR is calculated
VaR_days <- 5

# VaR alpha
VaR_alpha <- 0.05

# Monte Carlo sim
MC_n <- 500

# Conditional distribution GARCH: Students t distribution
GARCHcondDist <- "std"

## Run data handling file 
source('DataHandling.R')

# Precalculate some numbers for higher calculation performance
stock_n <- length(stockList)
stock_days <- nrow(stock_log)
pf_days <- length(pf_ret_nday)

## Get mean returns
stock_log_mean <- apply(stock_log,2,mean)

# Deduct the mean from every stock, for standardized news process in GARCH
stock_log_mean0 <- sweep(stock_log,2,stock_log_mean)

## Fit GARCH model on all marginal returns. Store all models in list. Also create vectors with GARCH parameters
GARCH_mu <- vector(mode = "list", length = stock_n)
GARCH_omega <- vector(mode = "list", length = stock_n)
GARCH_alpha <- vector(mode = "list", length = stock_n)
GARCH_beta <- vector(mode = "list", length = stock_n)
GARCH_residuals <- matrix(nrow = stock_days,ncol = stock_n)
GARCH_h.t <- matrix(nrow = stock_days,ncol = stock_n)
GARCH_sigma.t <- matrix(nrow = stock_days,ncol = stock_n)
GARCH_Z <- matrix(nrow = stock_days,ncol = stock_n)

for(s in 1:stock_n){
  GARCH <- garchFit(formula = ~ garch(1, 1), data=stock_log_mean0[,s], cond.dist=GARCHcondDist)
  GARCH_mu[s] <- coef(GARCH)[1]
  GARCH_omega[s] <- coef(GARCH)[2]
  GARCH_alpha[s] <- coef(GARCH)[3]
  GARCH_beta[s] <- coef(GARCH)[4]
  GARCH_residuals[,s] <- GARCH@residuals
  GARCH_h.t[,s] <- GARCH@h.t
  GARCH_sigma.t[,s] <- GARCH@sigma.t
}

GARCH_Z <- GARCH_residuals/GARCH_sigma.t

# Delete not needed variables
# rm(list='GARCH')

##### Copula

# Define the used copula
copula_obj <- claytonCopula(dim=length(stockList))

# Matrix of marginals, make universal
copula_pobs <- pobs(GARCH_Z)
copula_fit <- fitCopula(copula_obj,copula_pobs,method='ml')
copula_theta <- coef(copula_fit)

# Estimate marginals
copula_Z_mu <- vector(mode = "numeric", length = stock_n)
copula_Z_s <- vector(mode = "numeric", length = stock_n)
copula_Z_df <- vector(mode = "numeric", length = stock_n)

for(s in 1:stock_n){
  copula_Z_fit <- fitdistr(GARCH_Z[,s],"t")
  copula_Z_mu[s] <- copula_Z_fit$estimate[["m"]]
  copula_Z_s[s] <- copula_Z_fit$estimate[["s"]]
  copula_Z_df[s] <- copula_Z_fit$estimate[["df"]]
}

# Delete not needed variables
rm(list='copula_Z_fit')

# Parameters for copula function
copula_margins_list <- matrix(data='t',nrow=1,ncol=stock_n)[1,]

copula_paramMargins_list <- vector(mode="list", length = stock_n)
copula_paramMargins_list_names <- matrix(data='df',nrow=1,ncol=stock_n)[1,]
names(copula_paramMargins_list) <- copula_paramMargins_list_names

for(s in 1:stock_n){
  copula_paramMargins_list[s] <- copula_Z_df[s]
}

# Copula object
copula_dist <- mvdc(copula=claytonCopula(copula_theta, dim = length(stockList)), margins=copula_margins_list,
                    paramMargins = copula_paramMargins_list)

# Functions
GARCH_ht_function <- function(omega,alpha,beta,hPrevious,zPrevious){
  epsilon <- sqrt(hPrevious)*zPrevious
  h <- omega + alpha*epsilon^2 + beta*hPrevious
  return(h)
}


## For high calculating performance 4 dimensional arrays are used, to reduce the number of loops

# Dimensions: 1: the different stocks, 2: the total days of the model, 3: days of VaR, 4: Monte Carlo simulations for each day

# Draw all random residuals at once with correct dependence between stocks. In array the dimensions are filled consecutively with data.
# In the residuals matrix the first dimension is the different stocks, so that they can be filled with the correct dependence.

# Normal that this line takes a lot of time
MC_Z <- array(as.vector(t(rMvdc(copula_dist, n=MC_n*VaR_days*stock_days))), dim=c(stock_n,stock_days,VaR_days,MC_n))


# Scale back all standardized residuals
for(s in 1:stock_n){
  MC_Z[s,,,] <- MC_Z[s,,,]*copula_Z_s[s]+copula_Z_mu[s]
}  

## Check for dependence between modeled residuals

# Pairplot library
library(psych)

# Pairplot
pairs.panels(t(MC_Z[,1,1,]))

## Create array to store ht of simulations
MC_h <- array(dim=c(stock_n,stock_days,VaR_days,MC_n))

# The variance of the first simulated day is the same for all
for(s in 1:stock_n){
  MC_h[s,,1,] <- matrix(GARCH_h.t[,s],nrow=stock_days,ncol=MC_n)
}

# Apply GARCH function
for(s in 1:stock_n){
  for(i in 2:VaR_days){
    MC_h[s,,i,] <- GARCH_ht_function(GARCH_omega[[s]],GARCH_alpha[[s]],GARCH_beta[[s]],MC_h[s,,i-1,],MC_Z[s,,i-1,])
  }
}

## Simulate returns
MC_stock_log <- MC_Z*sqrt(MC_h)

# Add back the mean
for(s in 1:stock_n){
  MC_stock_log[s,,,] <- MC_stock_log[s,,,] + stock_log_mean[s]
}

# Get cumulative returns of the n days
MC_stock_log <- apply(MC_stock_log,c(1,2,4),sum)

# Change to simple returns
MC_stock_ret <- exp(MC_stock_log)-1

# Get portfolio returns
MC_ret <- apply(MC_stock_ret,c(2,3),mean)

# Change back to log returns
MC_log <- log(MC_ret+1)

# Apply quantile function to get VaR
VaR <- apply(MC_log,1,quantile,VaR_alpha)



### Check model
hitSeq       <- pf_ret_nday < VaR[1:length(pf_ret_nday)]
numberOfHits <- sum(hitSeq)
exRatio      <- numberOfHits/length(pf_ret_nday)


# Plot from class
index = index[1:length(pf_ret_nday)]
plot(index,pf_ret_nday, type="p")
lines(index,VaR[1:length(pf_ret_nday)], col="red" )
time = c(1:length(pf_ret_nday))
points(index[hitSeq ], pf_ret_nday[hitSeq], pch="+", col="green")


# Kupiec test
library(Rmpfr)

# Higher precision is needed, otherwise numerator and denumerator are treated as 0
N <- mpfr(length(pf_ret_nday),precBits= 128)
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



# GARCH test
GARCH <- garchFit(formula = ~ garch(1, 1), data=stock_log_mean0[,s], cond.dist=GARCHcondDist)
plot(GARCH@sigma.t,type='l')
test_h.t <- vector(mode = "numeric", length = length(GARCH@h.t))
h <- 0.02^2
test_h.t[1] <- h


GARCH_mu <- coef(GARCH)[1]
GARCH_omega <- coef(GARCH)[2]
GARCH_alpha <- coef(GARCH)[3]
GARCH_beta <- coef(GARCH)[4]
z <- GARCH@residuals/GARCH@sigma.t

for(i in 2:length(test_h.t)){
  test_h.t[i] <- GARCH_ht_function(GARCH_omega,GARCH_alpha,GARCH_beta,test_h.t[i-1],z[i])
}
lines(sqrt(test_h.t),col='green')

# TGARCH test
GARCH <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=stock_log_mean0[,5], cond.dist=GARCHcondDist)
plot(GARCH@sigma.t,type='l')
test_h.t <- vector(mode = "numeric", length = length(GARCH@h.t))
h <- 0.01^2
test_h.t[1] <- h


TGARCH_ht_function <- function(omega,alpha,beta,gamma,hPrevious,zPrevious){
  epsilon <- sqrt(hPrevious)*zPrevious
  h <- omega + alpha*epsilon^2 + gamma*pmin(epsilon,0)^2 +beta*hPrevious
  return(h)
}

GARCH_mu <- coef(GARCH)[1]
GARCH_omega <- coef(GARCH)[2]
GARCH_alpha <- coef(GARCH)[3]
GARCH_gamma <- coef(GARCH)[4]
GARCH_beta <- coef(GARCH)[5]

z <- GARCH@residuals/GARCH@sigma.t

for(i in 2:length(test_h.t)){
  test_h.t[i] <- TGARCH_ht_function(GARCH_omega,GARCH_alpha,GARCH_beta,GARCH_gamma,test_h.t[i-1],z[i])
}
lines(sqrt(test_h.t),col='green')

unconditionalVariance <- GARCH_omega/(1-GARCH_alpha-0.5*GARCH_gamma-GARCH_beta)



