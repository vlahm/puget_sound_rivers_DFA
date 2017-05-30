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

#seasonal Mann-Kendall tests on water temp series ####

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
# write.csv(output2, '../manuscript/figures/mann_kendall.csv', row.names=FALSE)

#plot final results (chisq for testing heterosked, z for testing slope?)
# par(mfrow=c(4,3))
# for(i in 2:ncol(ints)){
#     plot(ints[,i], col='gray80', main=colnames(ints)[i], ylab='temp C')
#     lines(temps[,i], col='gray10')
    # abline(a=output[(i-1),5], b=output[(i-1),2], col='red', lty=1) #???
    # abline(a=output[(i-1),5], b=output[(i-1),3], col='red', lty=2)
    # abline(a=output[(i-1),5], b=output[(i-1),4], col='red', lty=2)
# }

#MK tests on decomposed water temp series (stopped before finishing?) ####

trnds = ints
for(i in 2:ncol(ints)){
    trnds[,i] = decompose(ints[,i])$trend
}

output = matrix(NA, nrow=24, ncol=7, dimnames=list(colnames(temps)[2:25],
                                                    c('Kendall Tau','Slope','Lower_95','Upper_95','int','z','p_value')))
for(i in 2:ncol(trnds)){
    results = kendallTrendTest(trnds[7:450,i], season=c(7:12,rep(1:12, times=36),1:6), 
                               year=c(rep(1978,6),rep(1979:2014, each=12),rep(2015,6)),
                               independent.obs=FALSE)
    output[(i-1),c(6)] = results$statistic
    output[(i-1),c(7)] = results$p.value
    output[(i-1),c(1,2,5)] = results$estimate
    output[(i-1),3:4] = results$interval$limits
}

#write output table
# rownames(output) = colnames(temps)[-1]
chili = output
for(j in 1:ncol(output)){
    for(i in 1:nrow(output)){
        chili[i,j] = sprintf('%+1.3f', output[i,j])
    }
} #this doesnt render automatically in excel. gotta redo the 3-digit, +/- thing
# output2 = as.data.frame()
# output2 = as.data.frame(cbind(row.names(output2), output2$Slope, output2$Lower_95, output2$Upper_95,
                              # output2$`Kendall Tau`, output2$z, output2$p_value))
# colnames(output2) = c('Site','Slope','Lower 95','Upper 95','Kendall Tau','z','p value')
write.csv(chili, '../manuscript/figures/mann_kendall.csv', row.names=TRUE)

#MK tests on predictors (used in paper) ####

awd <- readRDS('../saved_structures/air_water_discharge2.rds')
pcsn <- readRDS('../saved_structures/precip_snowmelt.rds')
airT = decompose(ts(awd[[1]], frequency=12, start=c(1978,1)))$trend
precip = decompose(ts(pcsn[,1], frequency=12, start=c(1978,1)))$trend
Smelt = decompose(ts(pcsn[,2], frequency=12, start=c(1978,1)))$trend

# kendallSeasonalTrendTest(as.vector(awd[[1]]), season=rep(1:12, times=38), 
#                          year=rep(1978:2015, each=12), independent.obs=FALSE)
plot(airT)
airMK = kendallTrendTest(airT[7:450], season=c(7:12,rep(1:12, times=36),1:6), 
                         year=c(rep(1978,6),rep(1979:2014, each=12),rep(2015,6)),
                         independent.obs=FALSE)
plot(precip)
precipMK = kendallTrendTest(precip[7:450], season=c(7:12,rep(1:12, times=36),1:6), 
                         year=c(rep(1978,6),rep(1979:2014, each=12),rep(2015,6)),
                         independent.obs=FALSE)
plot(Smelt)
snowMK = kendallTrendTest(Smelt[7:450], season=c(7:12,rep(1:12, times=36),1:6), 
                         year=c(rep(1978,6),rep(1979:2014, each=12),rep(2015,6)),
                         independent.obs=FALSE)

covs_out = matrix(NA, nrow=7, ncol=3, 
                dimnames=list(c('Kendall Tau','Slope','Lower_95','Upper_95','int','z','p_value'),
                c('T_air','precip','snowmelt')))

covs_out[6,1] = airMK$statistic
covs_out[7,1] = airMK$p.value
covs_out[c(1,2,5),1] = airMK$estimate
covs_out[3:4,1] = airMK$interval$limits
covs_out[6,2] = precipMK$statistic
covs_out[7,2] = precipMK$p.value
covs_out[c(1,2,5),2] = precipMK$estimate
covs_out[3:4,2] = precipMK$interval$limits
covs_out[6,3] = snowMK$statistic
covs_out[7,3] = snowMK$p.value
covs_out[c(1,2,5),3] = snowMK$estimate
covs_out[3:4,3] = snowMK$interval$limits

#write output table
# options(scipen=100)
chili2 = covs_out
for(j in 1:ncol(covs_out)){
    for(i in 1:nrow(covs_out)){
        chili2[i,j] = sprintf('%+1.3f', covs_out[i,j])
    }
} 
# covs_out2 = as.data.frame(round(covs_out, 3))
write.csv(t(chili2), '../manuscript/figures/mann_kendall_covs.csv', row.names=TRUE)

#write combined output table
# comb = rbind(t(covs_out2), output2)
