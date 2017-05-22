#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 5/22/2017
#Seasonal Kendall test for monotonic trend in time series

rm(list=ls()); cat('\014')

#setup ####

#set working directory
setwd('C:/Users/Mike/git/stream_nuts_DFA/data/')
setwd('~/git/puget_sound_rivers_DFA/data')

#load raw temperature data
temps = readRDS('../saved_structures/raw_waterTemp.rds')

#install (if necessary) and load packages
package_list <- c('imputeTS','EnvStats')
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos="http://cran.rstudio.com/")
for(i in c(package_list)) library(i, character.only=TRUE)

#preliminary ####

#convert data to time series objects
for(i in 2:ncol(temps)){
    temps[,i] = ts(temps[,i], frequency=12, start=c(1978,1))
}

#interpolate NAs (turns out to be useful only for acf) and plot with loess
ints = temps
par(mfrow=c(4,3))
for(i in 2:ncol(temps)){
    ints[,i] = na.seasplit(temps[,i], 'interpolation')
    plot(ints[,i], col='red', main=colnames(ints)[i], ylab='temp C')
    lines(temps[,i], col='gray')
    lines(lowess(time(ints[,i]), ints[,i]), col='blue', lwd=1)
}

#test for autocorrelation (there's plenty, so we'll have to use block bootstrapping)
# par(mfrow=c(4,3))
# for(i in 2:ncol(temps)){
#     acf(ints[,i])
# }

#block-bootstrapped seasonal mann-kendall test (crappy approach)
# SMKtau <- function(x) SeasonalMannKendall(x)$tau
# for(i in 2:ncol(ints)){
#     boot.out <- tsboot(ints[,i], MKtau, R=500, l=nrow(ints)/3, sim="fixed")
#     boot.ci(boot.out, type="perc")
#     # print(SeasonalMannKendall(ints[,i]))
# }

#Seasonal Kendall Test for Trend Modified for Serial Correlation (good approach; interpolated data)
# output = matrix(NA, nrow=24, ncol=10, dimnames=list(colnames(ints)[2:25],
#                 c('tau','slope','LCL','UCL','int','chisq','chisq_p','z','z_p','df')))
# for(i in 2:ncol(ints)){
#     results = kendallSeasonalTrendTest(ints[,i], season=rep(1:12, times=38), 
#                              year=rep(1978:2015, each=12), independent.obs=FALSE)
#     output[(i-1),c(6,8)] = results$statistic
#     output[(i-1),10] = results$parameters
#     output[(i-1),c(7,9)] = results$p.value
#     output[(i-1),c(1,2,5)] = results$estimate
#     output[(i-1),3:4] = results$interval$limits
# }

#Seasonal Kendall Test for Trend Modified for Serial Correlation (good approach; UNinterpolated data)
output = matrix(NA, nrow=24, ncol=10, dimnames=list(colnames(temps)[2:25],
                                                    c('Kendall Tau','Slope','Lower_95','Upper_95','int','chisq','chisq_p','z','p_value','df')))
for(i in 2:ncol(temps)){
    results = kendallSeasonalTrendTest(temps[,i], season=rep(1:12, times=38), 
                                       year=rep(1978:2015, each=12), independent.obs=FALSE)
    output[(i-1),c(6,8)] = results$statistic
    output[(i-1),10] = results$parameters
    output[(i-1),c(7,9)] = results$p.value
    output[(i-1),c(1,2,5)] = results$estimate
    output[(i-1),3:4] = results$interval$limits
}

#write output table
options(scipen=100)
output2 = as.data.frame(round(output, 2))
output2 = as.data.frame(cbind(row.names(output2), output2$Slope, output2$Lower_95, output2$Upper_95, 
                output2$`Kendall Tau`, output2$z, output2$p_value))
colnames(output2) = c('Site','Slope','Lower 95','Upper 95','Kendall Tau','z','p value')
write.csv(output2, '../manuscript/figures/mann_kendall.csv', row.names=FALSE)

#plot final results (chisq for testing heterosked, z for testing slope?)
# par(mfrow=c(4,3))
# for(i in 2:ncol(ints)){
#     plot(ints[,i], col='gray80', main=colnames(ints)[i], ylab='temp C')
#     lines(temps[,i], col='gray10')
    # abline(a=output[(i-1),5], b=output[(i-1),2], col='red', lty=1) #???
    # abline(a=output[(i-1),5], b=output[(i-1),3], col='red', lty=2)
    # abline(a=output[(i-1),5], b=output[(i-1),4], col='red', lty=2)
# }
