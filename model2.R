libraries = c("fGarch")

lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x, repos = "http://cran.uni-muenster.de/")
})

lapply(c(libraries), require, character.only = TRUE)

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

## Copula




#### Code below for trial purpose

# Plot all standardized residuals
for(s in stockList){
  plot(get(paste(s,'Z',sep='')),type='l')
}


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

VaR <- MC_VaR_AllDays(gfit01@h.t[1:(length(gfit01@h.t)-1)],mu,omega,alpha1,beta1,Z,1,0.05)*sqrt(5)

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






