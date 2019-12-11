library("VineCopula")
library("copula")
library("fGarch")
library("MASS")

## Parameters

# Students t distribution
condDist <- "std"


## Fit GARCH on all marginal returns
for(s in stockList){
  assign(paste(s,'gfit',sep = ''),garchFit(formula = ~ garch(1, 1), data=get(s)$logReturns, cond.dist=condDist))
}

# Coefficients of all marginal returns
for(s in stockList){
  assign(paste(s,'mu',sep = ''),coef(get(paste(s,'gfit',sep = '')))[1])
  assign(paste(s,'omega',sep = ''),coef(get(paste(s,'gfit',sep = '')))[2])
  assign(paste(s,'alpha1',sep = ''),coef(get(paste(s,'gfit',sep = '')))[3])
  assign(paste(s,'beta',sep = ''),coef(get(paste(s,'gfit',sep = '')))[4])
}

# Standardized residuals and unconditional volatility for all marginal return distributions
for(s in stockList){
  assign(paste(s,'Z',sep = ''),(get(paste(s,'gfit',sep=''))@residuals-get(paste(s,'mu',sep='')))/get(paste(s,'gfit',sep=''))@sigma.t)
}

##### Copula
claytonCop <- claytonCopula(dim=length(stockList))

# m, make universal
m <- pobs(as.matrix(cbind(AMDZ,BAZ,MCDZ,WMTZ,XOMZ)))
fit <- fitCopula(claytonCop,m,method='ml')
theta <- coef(fit)

AMDZfit <- fitdistr(AMDZ,"t")
AMDZmu <- AMDZfit$estimate[["m"]]
AMDZs <- AMDZfit$estimate[["s"]]
AMDZdf <- AMDZfit$estimate[["df"]]

BAZfit <- fitdistr(BAZ,"t")
BAZmu <- BAZfit$estimate[["m"]]
BAZs <- BAZfit$estimate[["s"]]
BAZdf <- BAZfit$estimate[["df"]]

MCDZfit <- fitdistr(MCDZ,"t")
MCDZmu <- MCDZfit$estimate[["m"]]
MCDZs <- MCDZfit$estimate[["s"]]
MCDZdf <- MCDZfit$estimate[["df"]]

WMTZfit <- fitdistr(WMTZ,"t")
WMTZmu <- WMTZfit$estimate[["m"]]
WMTZs <- WMTZfit$estimate[["s"]]
WMTZdf <- WMTZfit$estimate[["df"]]

XOMZfit <- fitdistr(XOMZ,"t")
XOMZmu <- XOMZfit$estimate[["m"]]
XOMZs <- XOMZfit$estimate[["s"]]
XOMZdf <- XOMZfit$estimate[["df"]]



copula_dist <- mvdc(copula=claytonCopula(theta, dim = length(stockList)), margins=c('t','t','t','t','t'),
                    paramMargins = list(df=AMDZdf,df=BAZdf,df=MCDZdf,df=WMTZdf,df=XOMZdf))

simStandardizedResiduals <- rMvdc(copula_dist, n=1000)

simStandardizedResiduals[,1] <- simStandardizedResiduals[,1]*AMDZs+AMDZmu
simStandardizedResiduals[,2] <- simStandardizedResiduals[,2]*BAZs+BAZmu
simStandardizedResiduals[,3] <- simStandardizedResiduals[,3]*MCDZs+MCDZmu
simStandardizedResiduals[,4] <- simStandardizedResiduals[,4]*WMTZs+WMTZmu
simStandardizedResiduals[,5] <- simStandardizedResiduals[,5]*XOMZs+XOMZmu

colnames(sim) <- stockList

## Combine data for functions
htMat <- matrix(nrow=nrow(m),ncol=length(stockList))
muVec <- vector(mode='numeric',length(stockList))
omegaVec <- vector(mode='numeric',length(stockList))
alpha1Vec <- vector(mode='numeric',length(stockList))
betaVec <- vector(mode='numeric',length(stockList))

for(s in 1:length(stockList)){
  htMat[,s] <- get(paste(stockList[s],'gfit',sep=''))@h.t
  muVec[s] <- get(paste(stockList[s],'mu', sep=''))
  omegaVec[s] <- get(paste(stockList[s],'omega', sep=''))
  alpha1Vec[s] <- get(paste(stockList[s],'alpha1', sep=''))
  betaVec[s] <- get(paste(stockList[s],'beta', sep=''))
}
colnames(htMat) <- stockList
names(muVec) <- stockList
names(omegaVec) <- stockList
names(alpha1Vec) <- stockList
names(betaVec) <- stockList

## 5 day VaR

# Functions
GARCH_ht <- function(omega,alpha,beta,hPrevious,zPrevious){
  epsilon <- sqrt(hPrevious)*zPrevious
  h <- omega + alpha*epsilon^2 + beta*hPrevious
  return(h)
}

MC_VaR_OneDay <- function(h,mu,omega,alpha,beta,standardResid,days,VaR_alpha){
  # MC number simulations
  n <- 1000
  
  # Create a random residual matrix by sampling from standardized residuals
  #residMat <- matrix(data=sample(standardResid,n*days,replace = TRUE),ncol = days, nrow=n)
  residMat <- array(dim = c(n,days,length(stockList)))
  for(s in 1:length(stockList)){
    residMat[,,s] <- matrix(data=sample(standardResid[,s],n*days,replace = TRUE),ncol = days, nrow=n)
    # assign(paste(stockList[s],'residMat',sep = ''), matrix(data=sample(standardResid[,s],n*days,replace = TRUE),ncol = days, nrow=n))
  }
  
  # Create a matrix to store ht of simulations
  #hMat <- matrix(ncol = days,nrow = n)
  hMat <- array(dim = c(n,days,length(stockList)))
  
  # The variance of the first simulated day is the same for all
  for(i in 1:length(stockList)){
    hMat[,1,i] <- h[i]
  }
  
  
  # Apply the GARCH_ht function to get variance for each simulated day, depending on variance and residual of t-1
  if(days > 1){
    for(i in 1:length(stockList)){
      for(j in 2:days){
        hMat[,j,i] <- GARCH_ht(omega,alpha,beta,hMat[,j-1,i],residMat[,j-1,i])
      }
    }
  }
  
  # Simulate daily returns by multiplying the random residuals with the conditional volatility
  # Check if it runs without loop
  retMat <- array(dim = c(n,days,length(stockList)))
  for(s in 1:length(stockList)){
    for(j in 1:days){
      for(i in 1:n){
        retMat[i,j,s] <- residMat[i,j,s]*sqrt(hMat[i,j,s]) + mu[s]
      }
    }
  }
  # retMat <- residMat*sqrt(hMat) #+ mu
  
  # Cumulative log returns of simulated assets over all days
  cumulativeRet <- matrix(nrow = n,ncol=length(stockList))
  for(s in 1:length(stockList)){
    cumulativeRet[,s] <- apply(retMat[,,s],1,sum)
  }
  cumulativeRet <- exp(cumulativeRet)-1
  # Portfolio returns
  # change to simple returns
  
  portReturns <- apply(cumulativeRet,1,mean)

  #retMat <- rowSums(retMat)
  
  # Get the VaR of the given risk level with the quantile function
  VaR <- quantile(portReturns,VaR_alpha)
  
  return(VaR)
}

MC_VaR_AllDays <- function(htMat,muVec,omegaVec,alpha1Vec,betaVec,simStandardizedResiduals,days,VaR_alpha){
  VaR <- c(1:nrow(htMat))
  VaR <- 0
  for(i in 1:nrow(htMat)){
    VaR[i] <- MC_VaR_OneDay(htMat[i,],muVec,omegaVec,alpha1Vec,betaVec,simStandardizedResiduals,days,VaR_alpha)
  }
  return(VaR)
}

### 5 day VaR
#VaR <- MC_VaR_OneDay(htMat[1,],muVec,omegaVec,alpha1Vec,betaVec,simStandardizedResiduals,5,0.05)
VaR <- MC_VaR_AllDays(htMat,muVec,omegaVec,alpha1Vec,betaVec,simStandardizedResiduals,5,0.05)

# 5 Day portfolio returns
portReturns5days <- matrix(nrow=length(AMD5day$simpleReturns),ncol=length(stockList))
i = 1
for(s in stockList){
  portReturns5days[,i] <- get(paste(s,'5day',sep=''))$simpleReturns
  i = i+1
}
portReturns5days <- apply(portReturns5days,1,mean)

# Check
plot(portReturns5days)
lines(VaR[1:length(portReturns5days)],col='green')

hitSeq       <- portReturns5days < VaR[1:length(portReturns5days)]
numberOfHits <- sum(hitSeq)
exRatio      <- numberOfHits/length(AMD$logReturns)


### Nice plots

# Compare observed and simulated standardized returns
plot(AMDZ,XOMZ,main='Returns')
points(sim[,1],sim[,5],col='red')
legend('bottomright',c('Observed','Simulated'),col=c('black','red'),pch=21)

#### Code below for trial purpose

# Plot all standardized residuals
for(s in stockList){
  plot(get(paste(s,'Z',sep='')),type='l')
}

# Residuals vs. returns

plot(AMD$logReturns,type = 'l',col='green')
lines(AMDgfit@residuals)



### Test GARCH ht function

# TGARCH
#gfit01  <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=AMD$logReturns, cond.dist=condDist)

gfit01  <- garchFit(formula = ~ garch(1, 1), data=AMD$logReturns, cond.dist=condDist)

mu <- coef(gfit01)[1]
omega <- coef(gfit01)[2]
alpha1 <- coef(gfit01)[3]
beta1 <- coef(gfit01)[4]

# Standardize residuals
Z <- (gfit01@residuals-mu)/gfit01@sigma.t 
sigma <- sqrt(omega/(1-alpha1-beta1))

# Check residuals
plot(AMD$logReturns,type='l')
lines(gfit01@residuals,col='green')

plot(gfit01@residuals/sigma,type='l')
lines(Z,col='green')

# ht 0
ht0 <- omega/(1-alpha1-beta1)
ht0 <- 0.2
epsilon <- sqrt(ht0)*Z[1]

h <- c(NA)
length(h) <- length(Z)
h[1] <- ht0

for(i in 2:length(h)){
  h[i] <- omega + alpha1*epsilon^2 + beta1*h[i-1]
  epsilon <- sqrt(h[i])*Z[i]
}
s <- sqrt(h)


plot(gfit01@sigma.t ,type='l')
lines(s,col='green')



# One day VaR


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
      hMat[,j] <- GARCH_ht(omega,alpha,beta,hMat[,j-1],residMat[,j-1])/sqrt(days)
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

VaR <- MC_VaR_AllDays(gfit01@h.t[1:(length(gfit01@h.t))],mu,omega,alpha1,beta1,Z,1,0.05)*sqrt(5)

plot(AMD5day$logReturns)
lines(VaR,col='green')

hitSeq       <- AMD5day$logReturns< VaR
numberOfHits <- sum(hitSeq)
exRatio      <- numberOfHits/length(AMD$logReturns)






### Delete below


# TGARCH
#gfit01  <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=AMD$logReturns, cond.dist=condDist)

# Level of VaR
alpha  <- 0.05

# VaR
sqrtht <- volatility(gfit01)
nu     <- gfit01@fit$coef["shape"]
VaR    <- sqrtht*qt(alpha, nu)/sqrt(nu/(nu-2))

hitSeq       <- AMD$logReturns< VaR 
numberOfHits <- sum(hitSeq)
exRatio      <- numberOfHits/length(AMD$logReturns)

dev.new()
plot(AMD$logReturns, type="p")
lines(VaR, col="red" )
points(AMD$logReturns[hitSeq], pch="+", col="green")

# Get daily std from residuals
stdResiduals <- sqrt(var(gfit01@residuals))
muResiduals <- mean(gfit01@residuals)

# Standardize residuals
standardResiduals <- (gfit01@residuals-muResiduals)/stdResiduals
# epsilon <- gfit01@sigma.t*standardResiduals
test <- AMD$logReturns/gfit01@sigma.t



h <- omega/(1-alpha1-beta1)

hVector <- c(NA)
length(hVector) <- length(standardResiduals)
hVector[1] <- h

#for (i in 2:length(hVector)){
#  hVector[i] <- omega + alpha1*(epsilon[i-1]^2)+beta1*h
#  h <- hVector[i]
#}

for (i in 2:length(hVector)){
  epsilon <- sqrt(h)*standardResiduals[i-1]
  hVector[i] <- omega + alpha1*(epsilon^2)+beta1*h
  h <- hVector[i]
}
sigmaVector <- sqrt(hVector)

plot(sqrtht,type='l')
lines(sigmaVector, col='green')






