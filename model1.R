
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
pf_log = log(portfolio_ret + 1)


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






