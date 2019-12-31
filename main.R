# The following code creates the output file

# clean
rm(list=ls())

# libraries
library("fGarch")
library("MASS")
library("Rmpfr")
library("psych")
library("copula")

### Parameters

# Portfolio number (1: 5 stocks, 2: 10 stocks)
pf_n <- 1

# Number of days ahead the VaR is calculated
VaR_days <- 5

# VaR alpha
VaR_alpha <- 0.01

# number of Monte Carlo simulations by day
MC_n <- 1000

# Conditional distribution GARCH: Students t distribution
GARCHcondDist <- "std"

# GARCH model
#GARCH_model <- 'GARCH'
GARCH_model <- 'TGARCH'

# sample from observed returns or fit t distribution on standardized returns (for model 1)
#ret_method <- "sample"
ret_method <- "fit"


### Files

## Run data handling file 

source('DataHandling.R')

## Output Functions

source('OutputFile.R')

# Record start time
overall_time_start <- Sys.time()

# plot returns
plot(time, pf_log_nday, type="p")

## MODEL 1

source('model1.R')

# save VaR of model 1 for plotting
VaR1 <- VaR[1:length(time)]

## MODEL 2

source('model2.R')

# save VaR of model 2 for plotting
VaR2 <- VaR[1:length(time)]

## Results

# Record end time
overall_time_end <- Sys.time()
overall_time <- overall_time_end - overall_time_start

## Comparison Plots

# VaR and portfolio returns
plot_name <- paste('plots/VaR_model1_vs_model2_pf',toString(pf_n),'_',toString(VaR_days),'day_',toString(VaR_alpha),'VaR_',MC_n,'simulations_',strftime(Sys.time(),format = "%Y-%m-%d--%H-%M-%S"),'.pdf',sep='')
plot_main <- paste('Portfolio ',toString(pf_n),' - ',toString(VaR_days),' day ','- ',toString(VaR_alpha*100),'% VaR',sep='')
plot_ylab <- paste(toString(VaR_days),'day log returns')
plot_xlab <- 'Year'

hitSeq1 <- pf_log_nday < VaR1
hitSeq2 <- pf_log_nday < VaR2
differentExceedances <- hitSeq1 != hitSeq2

pdf(plot_name)
plot(index, pf_log_nday, type="p", main = plot_main,xlab = plot_xlab,ylab = plot_ylab)
lines(index, VaR1, col="red" )
lines(index, VaR2, col="blue" )
points(index[differentExceedances], pf_log_nday[differentExceedances], pch="+", col="green",cex=1.5)
legend('topright',c('VaR Model 1', 'VaR Model 2','Differences exceedance'),col=c('red', 'blue','green'), pch=c('-','-','+'))
dev.off()

# Histogram VaR
plot_name <- paste('plots/Histogram_VaR_model1_vs_model2_pf',toString(pf_n),'_',toString(VaR_days),'day_',toString(VaR_alpha),'VaR_',MC_n,'simulations_',strftime(Sys.time(),format = "%Y-%m-%d--%H-%M-%S"),'.pdf',sep='')
h1 <- hist(-VaR1,breaks=30,freq=FALSE)
h2 <- hist(-VaR2,breaks=30,freq=FALSE)
plot_ymax <- max(max(unlist(h1['density'])),max(unlist(h2['density'])))
plot_xmin <- min(min(unlist(h1['breaks'])),min(unlist(h2['breaks'])))
plot_xmax <- max(max(unlist(h1['breaks'])),max(unlist(h2['breaks'])))
pdf(plot_name)
plot(h1,col=rgb(1,0,0,1/2),ylim=c(0,plot_ymax),xlim=c(plot_xmin,plot_xmax),freq=FALSE,xlab='VaR',main='Histogram VaR')
plot(h2,col=rgb(0,0,1,1/2),ylim=c(0,plot_ymax),xlim=c(plot_xmin,plot_xmax),freq=FALSE,add=TRUE)
legend('topright',c('VaR Model 1', 'VaR Model 2'),col=c('red', 'blue'), lwd=2)
dev.off()

# Show Results
printFile()


