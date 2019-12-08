libraries = c("fGarch")

lapply(libraries, function(x) if (!(x %in% installed.packages())) {
  install.packages(x, repos = "http://cran.uni-muenster.de/")
})

lapply(c(libraries), require, character.only = TRUE)

# Students t distribution
condDist <- "std"

## GARCH
gfit01  <- garchFit(formula = ~ garch(1, 1), data=AMD$logReturns, cond.dist=condDist)

## TGARCH
#gfit01  <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=AMD$logReturns, cond.dist=condDist)

# Level of VaR
alpha  <- 0.05

# VaR
sqrtht <- gfit01@sigma.t
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
epsilon <- gfit01@sigma.t*gfit01@residuals
test <- AMD$logReturns/gfit01@sigma.t

mu <- gfit01@fit$matcoef[1,1]
omega <- gfit01@fit$matcoef[2,1]
alpha1 <- gfit01@fit$matcoef[3,1]
beta1 <- gfit01@fit$matcoef[4,1]

h <- sqrt(omega/(1-alpha1-beta1))

hVector <- c(NA)
length(hVector) <- length(epsilon)
hVector[1] <- h

for (i in 2:length(hVector)){
  hVector[i] <- omega + alpha1*(gfit01@sigma.t[i-1])+beta1*h
  h <- hVector[i]
}

