writeFile = function(model_nb) {
  
  outputFile_name <- paste('output/VaR_model',model_nb,'_pf',toString(pf_n),'_',toString(VaR_days),'day_',toString(VaR_alpha*100),'%VaR_',MC_n,'simulations_',strftime(Sys.time(),format = "%Y-%m-%d--%H-%M-%S"),'.txt',sep='')
  
  write("=====================================================",outputFile_name)
  write("",outputFile_name,append=TRUE)
  write("Parameters:",outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  
  write(paste("Portfolio number:              ",pf_n),outputFile_name,append=TRUE)
  write(paste("VaR model:                     ",toString(model_nb)),outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  
  write(paste("VaR Days:                      ",VaR_days),outputFile_name,append=TRUE)
  write(paste("VaR Alpha:                     ",VaR_alpha),outputFile_name,append=TRUE)
  write(paste("Number simulations:            ",MC_n),outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  
  write(paste("GARCH model:                   ",GARCH_model),outputFile_name,append=TRUE)
  write(paste("GARCH cond. dist:              ",GARCHcondDist),outputFile_name,append=TRUE)
  if (model_nb == 1) {
    write("Only for model 1:",outputFile_name, append= TRUE)
    write(paste("Returns generation method:     ",ret_method),outputFile_name,append=TRUE)
  }
  
  write("=====================================================",outputFile_name,append = TRUE)
  write("",outputFile_name,append=TRUE)
  write("Results:",outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  write(paste("Mean VaR:                      ",round(mean(-VaR),digits = 4)),outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  write(paste("Exceedance ratio:              ",round(as.numeric(exRatio),digits = 4)),outputFile_name,append=TRUE)
  write(paste("Kupiec K:                      ",round(K,digits = 4)),outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  if(K < qchisq(p,1)){
    write("VaR model is accurate at 99% level",outputFile_name,append=TRUE)
  }else{
    write("VaR model is not accurate at 99% level",outputFile_name,append=TRUE)
  }
  
  write("=====================================================",outputFile_name,append = TRUE)
  write("",outputFile_name,append=TRUE)
  write("Unconditional summary statistics:",outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  write("Simulated nday portfolio returns:",outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  write(paste("Mean:                          ",round(MC_log_mean,digits = 4)),outputFile_name,append=TRUE)
  write(paste("Std:                           ",round(MC_log_std,digits = 4)),outputFile_name,append=TRUE)
  write(paste("Skewness                       ",round(MC_log_skewness,digits = 4)),outputFile_name,append=TRUE)
  write(paste("Kurtosis:                      ",round(MC_log_kurtosis,digits = 4)),outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  write("Observed nday portfolio returns:",outputFile_name,append=TRUE)
  write("",outputFile_name,append=TRUE)
  write(paste("Mean:                          ",round(mean(pf_log_nday),digits = 4)),outputFile_name,append=TRUE)
  write(paste("Std:                           ",round(sqrt(var(pf_log_nday)),digits = 4)),outputFile_name,append=TRUE)
  write(paste("Skewness                       ",round(skewness(pf_log_nday),digits = 4)),outputFile_name,append=TRUE)
  write(paste("Kurtosis:                      ",round(kurtosis(pf_log_nday,method = 'moment'),digits = 4)),outputFile_name,append=TRUE)
  
  write("=====================================================",outputFile_name,append = TRUE)
  write(paste("Time passed:                   ",round(difftime(time_end,time_start,units='mins'),2),' minutes'),outputFile_name,append=TRUE)
}

printFile = function() {
  
  cat(paste('Portfolio',toString(pf_n),'for a',toString(VaR_days),'day',toString(VaR_alpha*100),'%VaR with',MC_n,'simulations.'),'\n')
  
  
  cat("=====================================================",'\n')
  cat("\n")
  cat("Parameters:",'\n')
  cat("\n")
  
  cat(paste("Portfolio number:              ",pf_n),'\n')
  cat("\n")
  
  
  cat(paste("VaR Days:                      ",VaR_days),'\n')
  cat(paste("VaR Alpha:                     ",VaR_alpha),'\n')
  cat(paste("Number simulations:            ",MC_n),'\n')
  cat("\n")
  
  cat(paste("GARCH model:                   ",GARCH_model),'\n')
  cat(paste("GARCH cond. dist:              ",GARCHcondDist),'\n')
  cat("Only for model 1:",'\n')
  cat(paste("Returns generation method:     ",ret_method),'\n')
  
  cat("=====================================================",'\n')
  cat("\n")
  cat("Results:",'\n')
  cat("\n")
  print(format(results[1:4,]))
  
  cat("=====================================================",'\n')
  cat("\n")
  cat("Unconditional summary statistics:",'\n')
  cat("\n")
  cat("Simulated nday portfolio returns:",'\n')
  cat("\n")
  print(format(results[5:8,]))
  cat("\n")
  cat("Observed nday portfolio returns:",'\n')
  cat("\n")
  cat(paste("Mean:                          ",round(mean(pf_log_nday),digits = 4)),'\n')
  cat(paste("Std:                           ",round(sqrt(var(pf_log_nday)),digits = 4)),'\n')
  cat(paste("Skewness                       ",round(skewness(pf_log_nday),digits = 4)),'\n')
  cat(paste("Kurtosis:                      ",round(kurtosis(pf_log_nday,method = 'moment'),digits = 4)),'\n')
  
  cat("=====================================================",'\n')
  cat(paste("Time passed:                   ",round(difftime(overall_time_end,overall_time_start,units='mins'),2),' minutes'),'\n')
}


results = data.frame(row.names=c("Mean VaR:","Exceedance ratio:", "Kupiec K:", "Chi-squared test:", "Sim Mean:", "Sim Std:", "Sim Skewness:", "Sim Kurtosis"))
# parameters = c(VaR_days, VaR_alpha, MC_n, GARCH_model, GARCHcondDist, ret_method)
# names(parameters)=c("VaR Days:", "VaR Alpha:", "Number of simulations:", "GARCH model:", "GARCH cond. dist:", "Returns generation method:")
# observed = c(round(mean(pf_log_nday),digits = 4), round(sqrt(var(pf_log_nday)),digits = 4),round(skewness(pf_log_nday),digits = 4),round(kurtosis(pf_log_nday,method = 'moment'),digits = 4))
# names(observed) = c("Obs Mean:", "Obs Std:", "Obs Skewness:", "Obs Kurtosis:")

fillData = function(model_nb) {
  
  if(K < qchisq(p,1)){
    chisqtest = "VaR model is accurate at 99% level"
  } else {
    chisqtest = "VaR model is not accurate at 99% level"
  }
  
  results[[paste("Model", model_nb)]] <- 
    c(
      round(mean(-VaR),digits = 4),
      round(as.numeric(exRatio),digits = 4),
      round(K,digits = 4),
      chisqtest,
      
      round(MC_log_mean,digits = 4),
      round(MC_log_std,digits = 4),
      round(MC_log_skewness,digits = 4),
      round(MC_log_kurtosis,digits = 4)
      )
  
  return(results)
}

