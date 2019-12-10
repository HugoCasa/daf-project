# This file reads all csv files in the Data folder and prepares them for further analysis. 
# The files are taken from Wharton Research Data Services. The data points Adjustment factor and 
# "Price - Close - Daily" are needed. 

rm(list=ls())

# List all files in Data folder
allStocks <- list.files('Data/')

# List to store all stock names
stockList <- c(NA)
length(stockList) <- length(allStocks)
i <- 1

for(s in allStocks){
  
  fileDir <- paste('Data/',s,sep = '')

  t <- read.csv(fileDir)
  name <- levels(t$tic)
  
  # Add name of stock to list
  stockList[i] <- name
  i <- i + 1
  
  # Get the adjusted closing price with the adjustment factor
  priceAdjusted = t$prccd/t$ajexdi
  
  # Get returns
  logReturns = diff(log(priceAdjusted))
  simpleReturns = priceAdjusted[2:length(priceAdjusted)]/priceAdjusted[1:length(priceAdjusted)-1] - 1
  
  # Assign to dataframe
  df <- data.frame(simpleReturns,logReturns)
  colnames(df) <- c('simpleReturns','logReturns')
  rownames(df) <- t$datadate[2:length(priceAdjusted)]
  
  # Create a dataframe named after the stock
  assign(name,df)
  
  ## Same for 5 day returns
  # Get returns
  logReturns = diff(log(priceAdjusted),lag=5)
  simpleReturns = priceAdjusted[6:length(priceAdjusted)]/priceAdjusted[1:(length(priceAdjusted)-5)] - 1
  
  # Assign to dataframe
  df <- data.frame(simpleReturns,logReturns)
  colnames(df) <- c('simpleReturns','logReturns')
  rownames(df) <- t$datadate[6:length(priceAdjusted)]
  
  # Create a dataframe named after the stock
  name <- paste(name,'5day',sep='')
  assign(name,df)
  
}

# Remove not needed variables
rm(list=c('df','t','allStocks','fileDir','name','priceAdjusted','s', 'i','logReturns','simpleReturns'))

