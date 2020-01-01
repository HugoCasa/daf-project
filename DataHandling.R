# This file reads all csv files in the Data folder and prepares them for further analysis. 
# The files are taken from Wharton Research Data Services. The data points Adjustment factor and 
# "Price - Close - Daily" are needed. 

# Get the right path for the chosen portfolio
pf_path <- paste('Data/Portfolio',toString(pf_n),'/',sep='')

# List all stocks in path of portfolio
allStocks <- list.files(pf_path)

# List to store all stock names
stockList <- c(NA)
length(stockList) <- length(allStocks)
i <- 1

for(s in allStocks){
  fileDir <- paste(pf_path,s,sep='')

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
  
  ## Same for n day returns
  # Get returns
  logReturns = diff(log(priceAdjusted),lag=VaR_days)
  simpleReturns = priceAdjusted[(VaR_days+1):length(priceAdjusted)]/priceAdjusted[1:(length(priceAdjusted)-(VaR_days))] - 1
  
  # Assign to dataframe
  df <- data.frame(simpleReturns,logReturns)
  colnames(df) <- c('simpleReturns','logReturns')
  rownames(df) <- t$datadate[(VaR_days+1):length(priceAdjusted)]
  
  # Create a dataframe named after the stock
  name <- paste(name,'nday',sep='')
  assign(name,df)
  
}

# Create dataframe of the indiviual stock returns
stock_ret <- matrix(nrow=length(get(stockList[1])$simpleReturns),ncol=length(stockList))
stock_log <- matrix(nrow=length(get(stockList[1])$logReturns),ncol=length(stockList))
i = 1
for(s in stockList){
  stock_ret[,i] <- get(s)$simpleReturns
  stock_log[,i] <- get(s)$logReturns
  i = i+1
}
colnames(stock_ret) <- stockList
colnames(stock_log) <- stockList

index <- as.Date(rownames(get(stockList[1])),"%Y%m%d")

### TEST
index <- 1:length(index)

stock_ret = as.data.frame(stock_ret, row.names = index)
stock_log = as.data.frame(stock_log, row.names = index)



## Portfolio returns

pf_ret <- as.numeric(apply(stock_ret,1,mean))
pf_log <- log(pf_ret+1)

# n Day portfolio returns
pf_ret_nday <- matrix(nrow=length(get(paste(stockList[1],'nday',sep=''))$simpleReturns),ncol=length(stockList))
i = 1
for(s in stockList){
  pf_ret_nday[,i] <- get(paste(s,'nday',sep=''))$simpleReturns
  i = i+1
}

# Can take n day simple returns to get portfolio returns because they were calculated end - beginning, not 1 day aggregated

# Mean of Assets
pf_ret_nday <- apply(pf_ret_nday,1,mean)

# Change to log returns
pf_log_nday <- log(pf_ret_nday+1)

# Remove not needed variables
for(s in stockList){
  rm(list=paste(s,'nday',sep=''))
  rm(list=s)
}

#

time <- as.Date(row.names(stock_ret))
time <- time[VaR_days:length(time)]


rm(list=c('df','t','allStocks','fileDir','name','priceAdjusted','s', 'i','logReturns','simpleReturns'))


