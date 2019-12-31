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

# Precalculate some numbers
pf_days <- length(pf_log)

# Get mean returns
pf_log_mean <- mean(pf_log)

# Deduct the mean from every stock, for standardized news process in GARCH
pf_log_mean0 <- pf_log - pf_log_mean


#=================================================================================================
## GARCH model: Fit GARCH model on all marginal returns


# For GARCH model
if (GARCH_model == "GARCH") {
  # Estimate GARCH for every individual stocks returns, adjusted for mean 0
  gfit01  <- garchFit(formula = ~ garch(1, 1), data=pf_log_mean0, cond.dist=GARCHcondDist)
  
  # store estimated parameters
  mu <- coef(gfit01)[["mu"]]
  omega <- coef(gfit01)[["omega"]]
  alpha <- coef(gfit01)[["alpha1"]]
  beta <- coef(gfit01)[["beta1"]]
}


# For TGARCH model
if (GARCH_model == "TGARCH") {
  # Estimate GARCH for every individual stocks returns, adjusted for mean 0
  gfit01 <- garchFit(formula = ~ garch(1, 1), delta =2, leverage = TRUE, data=pf_log_mean0, cond.dist=GARCHcondDist)
  
  # store estimated parameters
  mu <- coef(gfit01)[["mu"]]
  omega <- coef(gfit01)[["omega"]]
  alpha <- coef(gfit01)[["alpha1"]]
  beta <- coef(gfit01)[["beta1"]]
  gamma <- coef(gfit01)[["gamma1"]]
}

#=================================================================================================
## GARCH functions to simulate conditional variance h(t) as function of h(t-1) and z(t-1) and the estimated parameters

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

#=================================================================================================
## Monte Carlo simulation

# For high calculating performance 3 dimensional arrays are used, to reduce the number of loops

# Dimensions: 1: the total days of the model, 3: days of VaR, 4: Monte Carlo simulations for each day

#==============================================
## Standardized residuals

# This section either sample from the standardized resiudals or fit a student t distribution and generate values

# Standardized residuals are residuals divided by conditional volatility
Z <- gfit01@residuals/gfit01@sigma.t 

# store conditional volatility
h.t <- gfit01@h.t

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
  
  # sample from fitted student t distribution (3 dimensions: for each trading day, each nday VaR and each simulation per day)
  MC_Z <- array(rt(pf_days*VaR_days*MC_n,df=Z_df), dim=c(pf_days,VaR_days,MC_n))
  # Scale back all standardized residuals
  MC_Z <- MC_Z * Z_s + Z_mu
}

#==============================================
## Conditional variance

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


#=================================================================================================
## Simulate returns and get VaR

# Simulate returns
MC_log <- MC_Z*sqrt(MC_h)

# Add back the mean to every stock
MC_log <- MC_log + pf_log_mean

# Get cumulative returns of the n days
MC_log <- apply(MC_log,c(1,3),sum)


# Apply quantile function to get VaR
VaR <- apply(MC_log,1,quantile,VaR_alpha)

#=================================================================================================
## Check model

hitSeq       <- pf_log_nday < VaR[1:length(pf_log_nday)]
numberOfHits <- sum(hitSeq)
exRatio      <- numberOfHits/length(pf_log_nday)

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
## Plots

index = index[1:length(pf_log_nday)]

# Conditional volatility
plot_name <- 'plots/VaR_model1_volatility.pdf'
plot_xlab <- 'Year'
plot_ylab <- 'Std'
plot_main <- 'annual stdev.'
pdf(plot_name)
plot(index,sqrt(252)*gfit01@sigma.t[1:length(index)],type='l',main = plot_main,xlab = plot_xlab,ylab=plot_ylab)
plot_stock_uncond_vol <- vector(mode='numeric',length = length(index))
# Unconditional variance as mean of conditional variance
plot_stock_uncond_vol[] <- mean(sqrt(252)*gfit01@sigma.t[1:length(index)])
lines(index,plot_stock_uncond_vol,col='green',lwd = 2)
legend('topright',c(expression(paste(sqrt('h'['t']),': GARCH cond. stdev.')), expression(paste(sigma,':    uncond. stdev.'))),col=c('black', 'green'), lwd=2)
dev.off()

# Standardized residuals
plot_name <- 'plots/VaR_model1_standardized_residuals.pdf'
plot_xlab <- 'Year'
plot_ylab <- 'Std'
plot_main <- 'standardized residuals'
pdf(plot_name)
plot(index,Z[1:length(index)],type='l',main = plot_main,xlab = plot_xlab,ylab=plot_ylab)
dev.off()

# Student t fitting on the standardized residuals
if (ret_method == "fit") {
  plot_name <- 'plots/VaR_model1_hist_standardized_residuals.pdf'
  pdf(plot_name)
  hist(Z,breaks=80,main='Standardised Residuals',freq=F,col='blue')
  lines(seq(-10,10,0.01),dt((seq(-10,10,0.01)-Z_mu)/Z_s,Z_df)/Z_s,col='red',lwd=2)
  legend('topright',c('Standardized \nresiduals','Fitted t-distribution'),col=c('blue','red'), lwd=2)
  dev.off()
}

# VaR and portfolio returns
plot_name <- paste('plots/VaR_model1_pf',toString(pf_n),'_',toString(VaR_days),'day_',toString(VaR_alpha),'VaR_',MC_n,'simulations_',strftime(Sys.time(),format = "%Y-%m-%d--%H-%M-%S"),'.pdf',sep='')
plot_main <- paste('Portfolio ',toString(pf_n),' - ',toString(VaR_days),' day ','- ',toString(VaR_alpha*100),'% VaR',sep='')
plot_ylab <- paste(toString(VaR_days),'day log returns')
plot_xlab <- 'Year'
pdf(plot_name)
plot(index, pf_log_nday, main = plot_main,xlab = plot_xlab,ylab = plot_ylab)
lines(index, VaR[1:length(pf_log_nday)], col="red" )
points(index[hitSeq], pf_log_nday[hitSeq], pch="+", col="green")
legend('topright',c('VaR Model 1','Exceedances'),col=c('red','green'), pch=c('-','+'))
dev.off()

#=================================================================================================
## Calculate some unconditional summary statistics for simulated portfolio returns

# For simulated n day portfolio returns
MC_log_mean <- mean(MC_log)
MC_log_std <- sqrt(var(as.vector(MC_log)))
MC_log_skewness <- skewness(as.vector(MC_log))
MC_log_kurtosis <- kurtosis(as.vector(MC_log),method = 'moment')

#=================================================================================================
## Output file and console print

time_end <- Sys.time()
time_model <- time_end- time_start

writeFile(1)
results <- fillData(1)
