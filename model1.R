
# library for dataframe manipulation
library(dplyr)
library(fGarch)

N = nrow(get(stockList[1]))
M = length(stockList)
returns = as.data.frame(mget(stockList))
returns = select(returns, -contains("logReturns"))
colnames(returns) <- stockList


ret_mat = matrix(as.numeric(unlist(returns)), nrow=N, ncol=M)

pf_ret = apply(ret_mat, 1, mean)
pf_log = log(pf_ret + 1)


condDist <- 'std'
gfit01  <- garchFit(formula = ~ garch(1, 1), data=pf_log, cond.dist=condDist)

alpha = 0.05

sqrtht <- gfit01@sigma.t

nu     <- gfit01@fit$coef["shape"]
VaR    <- sqrtht*qt(alpha, nu)/sqrt(nu/(nu-2))


plot(pf_log)
lines(VaR, col='red')
hit = pf_log<VaR

nb_hits = sum(hit)
exc = nb_hits/length(hit)
var_mean = mean(VaR)

###### 5 days VaR

# Functions
GARCH_ht <- function(omega,alpha,beta,hPrevious,zPrevious){
  epsilon <- sqrt(hPrevious)*zPrevious
  h <- omega + alpha1*epsilon^2 + beta1*hPrevious
  return(h)
}

MC_VaR_OneDay <- function(h,mu,omega,alpha,beta,standardResid,days,VaR_alpha){
  # MC number simulations
  n <- 1000
  
  # Create a random residual matrix by sampling from standardized residuals
  residMat <- matrix(data=sample(standardResid,n,replace = TRUE),ncol = days, nrow=n)
  
  # Create a matrix to store ht of simulations
  hMat <- matrix(ncol = days,nrow = n)
  
  # The variance of the first simulated day is the same for all
  hMat[,1] <- h
  
  # Apply the GARCH_ht function to get variance for each simulated day, depending on variance and residual of t-1
  if(days > 1){
    for(j in 2:days){
      hMat[,j] <- GARCH_ht(omega,alpha,beta,hMat[,j-1],residMat[,j-1])
    }
  }
  # Simulate daily returns by multiplying the random residuals with the conditional volatility
  retMat <- residMat*sqrt(hMat) + mu
  
  # Cumulative log returns of all days
  retMat <- rowSums(retMat)
  
  # Get the VaR of the given risk level with the quantile function
  VaR <- quantile(retMat,VaR_alpha)
  
  return(VaR)
}



MC_VaR_AllDays <- function(ht,mu,omega,alpha,beta,standardResid,days,VaR_alpha){
  h <- c(1:length(ht))
  h <- 0
  for(i in 1:length(ht)){
    h[i] <- MC_VaR_OneDay(ht[i],mu,omega,alpha1,beta1,Z,days,VaR_alpha)
  }
  return(h)
}

# Try 5 day MC



mu <- coef(gfit01)[1]
omega <- coef(gfit01)[2]
alpha1 <- coef(gfit01)[3]
beta1 <- coef(gfit01)[4]
Z <- (gfit01@residuals-mu)/gfit01@sigma.t 

VaR5 <- MC_VaR_AllDays(gfit01@h.t[1:(length(gfit01@h.t)-1)],mu,omega,alpha1,beta1,Z,1,0.05)*sqrt(5)

plot(pf_log)
lines(VaR5,col='green')

hitSeq       <- pf_log< VaR5
numberOfHits <- sum(hitSeq)
exRatio      <- numberOfHits/length(pf_log)






