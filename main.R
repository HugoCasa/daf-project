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
VaR_alpha <- 0.05

# number of Monte Carlo simulations by day
MC_n <- 1000

# Conditional distribution GARCH: Students t distribution
GARCHcondDist <- "std"

# GARCH model
#GARCH_model <- 'GARCH'
GARCH_model <- 'TGARCH'

# sample or fit t distribution on standardized returns (for model 1)
ret_method <- "sample"
#ret_method <- "fit"


### Files

## Run data handling file 

source('DataHandling.R')

## Output Functions

source('OutputFile.R')

## MODEL 1

source('model1.R')

## MODEL 2

source('model2.R')