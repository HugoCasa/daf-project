libraries = c("fGarch")

lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x, repos = "http://cran.uni-muenster.de/")
})

lapply(c(libraries), require, character.only = TRUE)

# Students t distribution
condDist <- "std"

## GARCH
gfit01  <- garchFit(formula = ~ garch(1, 1), data=AMD$logReturns, cond.dist=condDist)

# TGARCH
#gfit01  <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=AMD$logReturns, cond.dist=condDist)

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






