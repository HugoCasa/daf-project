## Libraries

library("copula")
library("fGarch")
library("MASS")
library("psych")
library("Rmpfr")

rm(list=ls())

# Record start time
time_start <- Sys.time()

# ================================================================================================
## Parameters

# Portfolio number (1: 5 stocks, 2: 10 stocks)
pf_n <- 1

# Number of days ahead the VaR is calculated
VaR_days <- 5

# VaR alpha
VaR_alpha <- 0.05

# Monte Carlo simulation number
MC_n <- 500

# GARCH model
#GARCH_model <- 'GARCH'
GARCH_model <- 'TGARCH'

# Conditional distribution GARCH: Students t distribution
GARCHcondDist <- "std"

#=================================================================================================
## Data handling

pf_path <- paste('Data/Portfolio',toString(pf_n),'/',sep='')

# Run data handling file 
source('DataHandling.R')

# Precalculate some numbers
stock_n <- length(stockList)
stock_days <- nrow(stock_log)
pf_days <- length(pf_ret_nday)

# Get mean returns of the stocks
stock_log_mean <- apply(stock_log,2,mean)

# Deduct the mean from every stock, for standardized news process in GARCH process
stock_log_mean0 <- sweep(stock_log,2,stock_log_mean)

#=================================================================================================
## GARCH model: Fit GARCH model on all marginal returns

# Store all relevant parameters in variables:

# GARCH model parameters:
GARCH_mu <- vector(mode = "list", length = stock_n)
GARCH_omega <- vector(mode = "list", length = stock_n)
GARCH_alpha <- vector(mode = "list", length = stock_n)
GARCH_beta <- vector(mode = "list", length = stock_n)
GARCH_gamma <- vector(mode = "list", length = stock_n)

# GARCH residuals, condtional volatility and standardized residuals
GARCH_residuals <- matrix(nrow = stock_days,ncol = stock_n)
GARCH_h.t <- matrix(nrow = stock_days,ncol = stock_n)
GARCH_sigma.t <- matrix(nrow = stock_days,ncol = stock_n)
GARCH_Z <- matrix(nrow = stock_days,ncol = stock_n)

# For GARCH model
if(GARCH_model == 'GARCH'){
  for(s in 1:stock_n){
    # Estimate GARCH for every individual stocks returns, adjusted for mean 0
    # Use fGARCH package. Normal GARCH model is a special case of fGARCH model.
    GARCH <- garchFit(formula = ~ garch(1, 1), data=stock_log_mean0[,s], cond.dist=GARCHcondDist)
    
    # Store estiated parameters in lists
    GARCH_mu[s] <- coef(GARCH)[1]
    GARCH_omega[s] <- coef(GARCH)[2]
    GARCH_alpha[s] <- coef(GARCH)[3]
    GARCH_beta[s] <- coef(GARCH)[4]
    GARCH_residuals[,s] <- GARCH@residuals
    GARCH_h.t[,s] <- GARCH@h.t
    GARCH_sigma.t[,s] <- GARCH@sigma.t
  }
  # Standardized residuals are residuals divided by conditional volatility
  GARCH_Z <- GARCH_residuals/GARCH_sigma.t
}

# For TGARCH model
if(GARCH_model == 'TGARCH'){
  for(s in 1:stock_n){
    # Estimate TGARCH for every individual stocks returns, adjusted for mean 0
    # Use fGARCH package. GJR-GARCH model is a special case of fGARCH model with delta = 2 and leverage
    GARCH <- garchFit(formula = ~ garch(1, 1), delta = 2, leverage = TRUE, data=stock_log_mean0[,s], cond.dist=GARCHcondDist)
    
    # Store estiated parameters in lists
    GARCH_mu[s] <- coef(GARCH)[1]
    GARCH_omega[s] <- coef(GARCH)[2]
    GARCH_alpha[s] <- coef(GARCH)[3]
    GARCH_gamma[s] <- coef(GARCH)[4]
    GARCH_beta[s] <- coef(GARCH)[5]
    GARCH_residuals[,s] <- GARCH@residuals
    GARCH_h.t[,s] <- GARCH@h.t
    GARCH_sigma.t[,s] <- GARCH@sigma.t
  }
  # Standardized residuals are residuals divided by conditional volatility
  GARCH_Z <- GARCH_residuals/GARCH_sigma.t
}

# Delete not needed variables to free working memory
rm(list='GARCH')

#=================================================================================================
## Copula

# This section creates a copula model in order to get the multivariate distribution of the marginal 
# returns. The copula package is used.

# Define the used copula
copula_obj <- claytonCopula(dim=length(stockList))

# Matrix of marginal returns as pseudo observations [0,1]
copula_pobs <- pobs(GARCH_Z)

# Fit the copula
copula_fit <- fitCopula(copula_obj,copula_pobs,method='ml')

# Get the fitted copula parameter
copula_theta <- coef(copula_fit)

# Fit a t distribution on standardised redisuals of marginal returns and store parameters in vectors
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

## Create lists with parameters for use of copula function

# T distribution for all standardized residuals
copula_margins_list <- matrix(data='t',nrow=1,ncol=stock_n)[1,]

# Make a list with the degrees of freedom of the standardized residuals
copula_paramMargins_list <- vector(mode="list", length = stock_n)
copula_paramMargins_list_names <- matrix(data='df',nrow=1,ncol=stock_n)[1,]
names(copula_paramMargins_list) <- copula_paramMargins_list_names
for(s in 1:stock_n){
  copula_paramMargins_list[s] <- copula_Z_df[s]
}

# Create a copula object with all correct parameters, that can later be used to generate an array with returns
copula_dist <- mvdc(copula=claytonCopula(copula_theta, dim = length(stockList)), margins=copula_margins_list,
                    paramMargins = copula_paramMargins_list)

#=================================================================================================
## GARCH functions to simulate conditional variance h(t) as function of h(t-1) and z(t-1) and the estimated parameters

# GARCH
GARCH_ht_function <- function(omega,alpha,beta,hPrevious,zPrevious){
  epsilon <- sqrt(hPrevious)*zPrevious
  h <- omega + alpha*epsilon^2 + beta*hPrevious
  return(h)
}

# GJR GARCH as special version of fGARCH. From introduction to the rugarch package, p.9
TGARCH_ht_function <- function(omega,alpha,beta,gamma,hPrevious,zPrevious){
  h <- omega + alpha*hPrevious*(abs(zPrevious)-gamma*zPrevious)^2 + beta*hPrevious
  return(h)
}

#=================================================================================================
## Monte Carlo simulation

# For high calculating performance 4 dimensional arrays are used, to reduce the number of loops

# Dimensions: 1: the different stocks, 2: the total days of the model, 3: days of VaR, 4: Monte Carlo simulations for each day

#==============================================
## Standardized residuals

# Draw all random residuals at once with correct dependence between stocks. Use the copula object generated before.
# In array the dimensions are filled consecutively with data. In the MC_Z array the first dimension are 
# the different stocks, so that they can be filled with the correct dependence. The residuals in other 3 dimensions are i.i.d.

# Normal that this line takes a lot of time
MC_Z <- array(as.vector(t(rMvdc(copula_dist, n=MC_n*VaR_days*stock_days))), dim=c(stock_n,stock_days,VaR_days,MC_n))

# Scale back all standardized residuals
for(s in 1:stock_n){
  MC_Z[s,,,] <- MC_Z[s,,,]*copula_Z_s[s]+copula_Z_mu[s]
}  

## Check correct sampling with dependence between modeled standardised residuals

# Pairplot
MC_residuals_plot <- t(MC_Z[,1,1,])
colnames(MC_residuals_plot) <- stockList
pairs.panels(t(MC_Z[,1,1,]))

#==============================================
## Conditional variance

# Create array to store h(t)
MC_h <- array(dim=c(stock_n,stock_days,VaR_days,MC_n))

# The variance of the first simulated day is taken from GARCH model and is the same for all simulated days
for(s in 1:stock_n){
  MC_h[s,,1,] <- matrix(GARCH_h.t[,s],nrow=stock_days,ncol=MC_n)
}

# Apply GARCH or TGARCH function. Function works over the 2 dimensions all days of the model 
# and all simulated days at the same time
if(GARCH_model == 'GARCH'){
  for(s in 1:stock_n){
    for(i in 2:VaR_days){
      MC_h[s,,i,] <- GARCH_ht_function(GARCH_omega[[s]],GARCH_alpha[[s]],GARCH_beta[[s]],MC_h[s,,i-1,],MC_Z[s,,i-1,])
    }
  }
}
if(GARCH_model == 'TGARCH'){
  for(s in 1:stock_n){
    for(i in 2:VaR_days){
      MC_h[s,,i,] <- TGARCH_ht_function(GARCH_omega[[s]],GARCH_alpha[[s]],GARCH_beta[[s]],GARCH_gamma[[s]],MC_h[s,,i-1,],MC_Z[s,,i-1,])
    }
  }
}

#=================================================================================================
## Simulate returns and get VaR

# Multiply sampled standardized returns with conditional volatility. Works over all 4 dimensions
MC_stock_log <- MC_Z*sqrt(MC_h)

# Add back the mean
for(s in 1:stock_n){
  MC_stock_log[s,,,] <- MC_stock_log[s,,,] + stock_log_mean[s]
}

# Get cumulative returns of the n days. Apply function on the other 3 dimensions
MC_log <- apply(MC_stock_log,c(1,2,4),sum)

# Change to simple returns because they aggregate over assets
MC_ret <- exp(MC_log)-1

# Get portfolio returns as mean of the returns of the different stocks
MC_ret <- apply(MC_ret,c(2,3),mean)

# Change back to log returns
MC_log <- log(MC_ret+1)

# Apply quantile function to get VaR
VaR <- apply(MC_log,1,quantile,VaR_alpha)


#=================================================================================================
## Check model

# Exceedance ratio, code as shown in class
hitSeq       <- pf_ret_nday < VaR[1:length(pf_ret_nday)]
numberOfHits <- sum(hitSeq)
exRatio      <- numberOfHits/length(pf_ret_nday)


# Plot returns and VaR, code as shown in class
index = index[1:length(pf_ret_nday)]
plot(index,pf_ret_nday, type="p")
lines(index,VaR[1:length(pf_ret_nday)], col="red" )
time = c(1:length(pf_ret_nday))
points(index[hitSeq], pf_ret_nday[hitSeq], pch="+", col="green")


## Kupiec test

# Higher precision is needed, otherwise numerator and denumerator are treated as 0
N <- mpfr(length(pf_ret_nday),precBits= 128)
exRatio <- mpfr(exRatio,precBits = 128)
numberOfHits <- mpfr(numberOfHits,precBits = 128)
VaR_alpha <- mpfr(VaR_alpha,precBits = 128)
num <- (exRatio^numberOfHits)*(1-exRatio)^(N-numberOfHits)
den <- (VaR_alpha^numberOfHits )*(1-VaR_alpha)^(N-numberOfHits)
K   <- as.numeric(2*log(num/den))
VaR_alpha <- as.numeric(VaR_alpha)
exRatio <- as.numeric(exRatio)
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
hist(as.vector(MC_log))

#=================================================================================================
## Output file
time_end <- Sys.time()
time_model <- time_end- time_start

outputFile_name <- paste('output/VaR_model2_pf',toString(pf_n),'_',toString(VaR_days),'day_',toString(VaR_alpha*100),'%VaR_',MC_n,'simulations_',strftime(Sys.time(),format = "%Y-%m-%d--%H-%M-%S"),'.txt',sep='')
source('OutputFile.R')

