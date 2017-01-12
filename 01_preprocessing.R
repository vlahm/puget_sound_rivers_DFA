
rm(list=ls())

library(stringr)

nutsfull <- read.csv("C:/Users/Mike/git/stream_nuts_DFA/data/chemPhys_data/nutsdata.csv",
                     stringsAsFactors=FALSE)
# nutsfull <- read.csv("~/Desktop/grad/Projects/Thesis/stream_nuts_DFA/data/nutsdata.csv",
#                      stringsAsFactors=FALSE)
nutsfull[,3:(ncol(nutsfull)-2)] <- apply(nutsfull[,3:(ncol(nutsfull)-2)], 2, as.numeric)

#fix a few date errors ####
for(i in c(5084:5531,11063:11447)){
yr <- as.numeric(str_match(nutsfull$dateTime[i], "\\d+/\\d+/(\\d+)")[2])
mo <- as.numeric(str_match(nutsfull$dateTime[i], "(\\d+)/\\d+/\\d+")[2])
da <- as.numeric(str_match(nutsfull$dateTime[i], "\\d+/(\\d+)/\\d+")[2])

  if(yr > 50){
      if(nchar(mo) == 1 & nchar(da) == 1){
        nutsfull$dateTime[i] <- gsub('(\\d+)/(\\d+)/(\\d+) (\\d+):(\\d+)',
                                             '19\\3-0\\1-0\\2 \\4:\\5:00', nutsfull$dateTime[i])
      } else {
          if(nchar(mo) == 1 & nchar(da) == 2){
              nutsfull$dateTime[i] <- gsub('(\\d+)/(\\d+)/(\\d+) (\\d+):(\\d+)',
                                           '19\\3-0\\1-\\2 \\4:\\5:00', nutsfull$dateTime[i])
          } else {
              if(nchar(mo) == 2 & nchar(da) == 1){
                  nutsfull$dateTime[i] <- gsub('(\\d+)/(\\d+)/(\\d+) (\\d+):(\\d+)',
                                               '19\\3-\\1-0\\2 \\4:\\5:00', nutsfull$dateTime[i])
              } else {
                  nutsfull$dateTime[i] <- gsub('(\\d+)/(\\d+)/(\\d+) (\\d+):(\\d+)',
                                               '19\\3-\\1-\\2 \\4:\\5:00', nutsfull$dateTime[i])
              }
          }
      }
  } else {
      if(nchar(mo) == 1 & nchar(da) == 1){
          nutsfull$dateTime[i] <- gsub('(\\d+)/(\\d+)/(\\d+) (\\d+):(\\d+)',
                                       '20\\3-0\\1-0\\2 \\4:\\5:00', nutsfull$dateTime[i])
      } else {
          if(nchar(mo) == 1 & nchar(da) == 2){
              nutsfull$dateTime[i] <- gsub('(\\d+)/(\\d+)/(\\d+) (\\d+):(\\d+)',
                                           '20\\3-0\\1-\\2 \\4:\\5:00', nutsfull$dateTime[i])
          } else {
              if(nchar(mo) == 2 & nchar(da) == 1){
                  nutsfull$dateTime[i] <- gsub('(\\d+)/(\\d+)/(\\d+) (\\d+):(\\d+)',
                                               '20\\3-\\1-0\\2 \\4:\\5:00', nutsfull$dateTime[i])
              } else {
                  nutsfull$dateTime[i] <- gsub('(\\d+)/(\\d+)/(\\d+) (\\d+):(\\d+)',
                                               '20\\3-\\1-\\2 \\4:\\5:00', nutsfull$dateTime[i])
              }
          }
      }
  }
}

# 1 CHOOSE whether to aggregate by month or year (and then mean or max) ####

# 1a aggregate by month
nuts <- aggregate(nutsfull[,seq(3,35,2)], by=list(nutsfull$siteCode, substr(nutsfull$dateTime, 1,7)),
          FUN=mean, na.rm=TRUE, na.action=NULL)
# 1b aggregate by year (means)
nuts <- aggregate(nutsfull[,seq(3,35,2)],
                  by=list(nutsfull$siteCode, substr(nutsfull$dateTime, 1,4)),
                  FUN=mean, na.rm=TRUE, na.action=NULL)
# 1c aggregate by year (maxes)
nuts <- aggregate(nutsfull[,seq(3,35,2)],
                  by=list(nutsfull$siteCode, substr(nutsfull$dateTime, 1,4)),
                  FUN=max, na.rm=TRUE, na.action=NULL)


#replace NaNs and -Infs with NA
nuts[,3:ncol(nuts)][apply(nuts[,3:ncol(nuts)], 2, is.nan)] <- NA
nuts[,3:ncol(nuts)][apply(nuts[,3:ncol(nuts)], 2, is.infinite)] <- NA

# 2 find out proportions of NAs in each column####
checkNA <- function(){
    na_prop <- rep(NA, ncol(nuts))
    for(i in 1:ncol(nuts)){
        na_prop[i] <- sum(is.na(nuts[,i]))/nrow(nuts)
    }
    out <- cbind(colnames(nuts), na_prop)
    return(out)
}
checkNA()

# 3 throw away columns with lots of NAs####
nuts <- nuts[,-c(5,15,16,17,19)] #works for both year and month
colnames(nuts)[1:2] <- c('site', 'date')

# 4 assemble response dataframes####

#this will help with sorting
temp <- factor(nuts$site)
levels(temp)[2:5] <- c('ZA','ZC','ZD','ZE')
nuts$site <- as.vector(temp)
moreLETTERS <- c(LETTERS,'ZA','ZC','ZD','ZE')

#this builds properly structured dataframes for each chem metric of interest
assembler <- function(colname){
    colnum <- which(colnames(nuts)==colname)
    for(i in unique(nuts$site)){
        assign(paste0('site', i), nuts[which(nuts$site == i), c(2,colnum)])
    }

    #some day find a way to do this elegantly. a separate 'assign' call runs, but doesn't work
    colnames(siteA)[2] <- 'A'; colnames(siteB)[2] <- 'B'; colnames(siteC)[2] <- 'C'
    colnames(siteD)[2] <- 'D'; colnames(siteE)[2] <- 'E'; colnames(siteF)[2] <- 'F'
    colnames(siteG)[2] <- 'G'; colnames(siteH)[2] <- 'H'; colnames(siteI)[2] <- 'I'
    colnames(siteJ)[2] <- 'J'; colnames(siteK)[2] <- 'K'; colnames(siteL)[2] <- 'L'
    colnames(siteM)[2] <- 'M'; colnames(siteN)[2] <- 'N'; colnames(siteO)[2] <- 'O'
    colnames(siteP)[2] <- 'P'; colnames(siteQ)[2] <- 'Q'; colnames(siteR)[2] <- 'R'
    colnames(siteS)[2] <- 'S'; colnames(siteT)[2] <- 'T'; colnames(siteU)[2] <- 'U'
    colnames(siteV)[2] <- 'V'; colnames(siteW)[2] <- 'W'; colnames(siteX)[2] <- 'X'
    colnames(siteY)[2] <- 'Y'; colnames(siteZ)[2] <- 'Z'; colnames(siteZA)[2] <- 'ZA'
    colnames(siteZC)[2] <- 'ZC'; colnames(siteZD)[2] <- 'ZD'; colnames(siteZE)[2] <- 'ZE'

    for(i in 2:length(moreLETTERS)){
        siteA <- merge(siteA, eval(parse(text=paste0('site', moreLETTERS[i]))), by='date', all=TRUE)
    }

    #format dates and set to middle of the year for plotting
    siteA$date <- as.Date(paste0(siteA$date, '-06-15'))

    return(siteA)
}

COND <- assembler('COND')
FC <- assembler('FC')
NH3_N <- assembler('NH3_N')
NO2_NO3 <- assembler('NO2_NO3')
OP_DIS <- assembler('OP_DIS')
OXYGEN <- assembler('OXYGEN')
PH <- assembler('PH')
PRESS <- assembler('PRESS')
SUSSOL <- assembler('SUSSOL')
TEMP <- assembler('TEMP')
TP_P <- assembler('TP_P')
TURB <- assembler('TURB')

# 5 dump to disk (change filename accordingly)####
save(COND, FC, NH3_N, NO2_NO3, OP_DIS, OXYGEN, PH, PRESS, SUSSOL, TEMP, TP_P, TURB,
     list=c('COND', 'FC', 'NH3_N', 'NO2_NO3', 'OP_DIS', 'OXYGEN', 'PH', 'PRESS', 'SUSSOL', 'TEMP', 'TP_P', 'TURB'),
     file="C:/Users/Mike/git/stream_nuts_DFA/data/chemPhys_data/yys_by___.rda")

# 6 add temp, precip, drought, heating degree days ####
#data from http://www.ncdc.noaa.gov/cag/time-series/us/45/3/zndx/ytd/12/1970-2014?base_prd=true&firstbaseyear=1970&lastbaseyear=2014
#hydro drought is PHDI, meteor drought is PMDI, shortTerm is Z-index

#by year

#read in temp and convert to C
# meantemp <- read.csv("C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_year/ hdd3 .csv",
#          skip=3)
# meantemp$Date <- 1970:2014
# meantemp$Value <- round((meantemp$Value - 32) * (5/9), 2)
# meantemp$Anomaly <- meantemp$Value - mean(meantemp$Value)
# colnames(meantemp)[2:3] <- c('meantemp', 'meantemp_anom')
#
# #read in precip (unit = inches)
# precip <- read.csv("C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\stream_nuts_DFA\\data\\precip_byYear.csv",
#          skip=3)
# precip$Date <- 1970:2014
# colnames(precip)[2:3] <- c('precip', 'precip_anom')
#
# #hydrological drought (PHDI)
# hydroDrought <- read.csv("C:\\Users\\Mike\\Desktop\\Grad\\Projects\\Thesis\\stream_nuts_DFA\\data\\hydrol_drought_byYear.csv",
#          skip=2)
# hydroDrought$Date <- 1970:2014
# colnames(hydroDrought)[2:3] <- c('hydroDrought', 'hydroDrought_anom')
#
# #combine
# climate <- Reduce(function(x,y) merge(x,y, by='Date'),
#                   x=list(meantemp, precip, hydroDrought, meteoDrought, ZDrought))
#
# save(climate, list='climate', file="C:/Users/Mike/Desktop/Grad/Projects/Thesis/stream_nuts_DFA/data/climate.rda")


