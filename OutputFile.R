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
  write(paste("Mean VaR:                      ",round(mean(VaR),digits = 4)),outputFile_name,append=TRUE)
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
