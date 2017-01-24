#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 1/23/17

cat('\014'); rm(list=ls())

# 0 - setup ####
setwd('C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/snow_data/raw')

if (!require("devtools")) install.packages("devtools")
library(devtools)
if (!require("manipulateR")) install_github('vlahm/manipulateR')
if (!require("stringr")) install.packages("stringr")
if (!require("imputeTS")) install.packages("imputeTS")

library(manipulateR)
library(stringr)
library(imputeTS)

#read raw snowtel data
x <- ultimate_reader(dir_args=list(path='./', pattern='[^00]\\.csv'),
                read_args=list(sep=',', quote="\"", header=TRUE, stringsAsFactors=FALSE,
                               fill=TRUE, comment.char='', skip=58),
                merge=TRUE, by='Date', all=TRUE)

# 1 - trim and sort dataset

#grab just the relevant measurement
x <- x[,c(1, which(grepl('Change.In.Snow.Water.Equivalent', colnames(x))))]

#rename cols by site
colnames(x) <- c('date', str_match(dir(pattern='[^00]\\.csv'), '(.*).csv')[,2])

#format date column
year <- substr(x$date, 5, 8)
month_num <- as.numeric(as.character(factor(substr(x$date, 1, 3), labels=order(month.abb))))
month_num <- sprintf('%02.f', month_num)
x$date <- paste0(year, month_num)

#sort by date
x <- x[order(substr(x$date,1,4), substr(x$date,5,6)),]
rownames(x) <- 1:nrow(x)

# 2 - get mean and sd snowfall

#means and sds by row
swe <- apply(x[,2:ncol(x)], 1, mean, na.rm=TRUE)
sds <- apply(x[,2:ncol(x)], 1, sd, na.rm=TRUE)

#number of NAs per row (max 6)
NAs <- apply(x[,2:ncol(x)], 1, function(i) sum(is.na(i)))
sum(NAs > 3) #not bad coverage, also:

#decide whether to aggregate across all sites
sum(abs(range(swe, na.rm=TRUE))) #range of means
sum(abs(range(sds, na.rm=TRUE))) #range of sds
mean(abs(swe)) #mean
abs_mean <- apply(x[,2:ncol(x)], 1, function(i) mean(abs(i), na.rm=TRUE))
sum(sds>=abs_mean, na.rm=TRUE) #sds higher than abs(means)

# 3 - impute 1978 values and remove snow accumulation values
swe <- c(rep(NA, 9), swe[-1])
swe_ts <- ts(swe, start=1, frequency=12)

#black shows the imputed section
swe_tsi <- na.seasplit(swe_ts, 'interpolation')
plot(swe_tsi, col='black', type='l')
lines(swe_ts, col='gray')

#set all positive values (accumulations) to zero, and flip the sign of snowmelt
swe <- as.vector(swe_tsi)
swe[swe>0] <- 0
swe <- abs(swe)

# 4 - build output dataset and export
out <- data.frame(date=c(paste0('1978', sprintf('%02.f', 1:8)), x$date),
                  snowmelt=swe)

write.csv(out, '../snowmelt.csv', row.names=FALSE)
