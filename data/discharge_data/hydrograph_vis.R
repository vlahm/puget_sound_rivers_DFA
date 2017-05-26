#visualize hydrographs

rm(list=ls()); cat('\014')


#setup ####
setwd('C:/Users/Mike/git/stream_nuts_DFA/data/discharge_data')

if(!require('lubridate')) install.packages('lubridate')
library(lubridate)

dis = read.csv('discharge.csv', colClasses=c('date'='Date'))

#preprocess ####

#subset, make class labels
dis = dis[year(dis$date) >= 1978 & year(dis$date) <= 2015,]
dis = subset(dis, select=-c(X))

classes = c('RS','RD','SD','RD','SD','RS','SD','SD','RS','RD','RS',
            'RS','RS','RS','RS','RS','SD','SD','RD')

cols = as.character(factor(classes, labels=c('darkolivegreen3','brown','lightblue')))

#aggregate by monthly discharge
dis_agg = aggregate(dis, by=list(substr(dis$date,6,7)), FUN=mean, na.rm=TRUE)

#scale discharge
dis_agg[3:21] = scale(dis_agg[,3:21])

#plot ####

#plot all together
plot(dis_agg[,3], type='l', col=cols[1], ylim=c(min(dis_agg[,3:21]), max(dis_agg[,3:21])))
for(i in 4:21){
    lines(dis_agg[,i], col=cols[i-1])
}

#indiv plots
defpar = par(mfrow=c(4,5))
for(i in 3:21){
plot(dis_agg[,i], type='l', col=cols[i-2], 
     ylim=c(min(dis_agg[,3:21]), max(dis_agg[,3:21])),
     main=colnames(dis_agg)[i], lwd=2)
}
par(defpar)
