
#compose data from individual sites ####
setwd("C:\\Users\\Mike\\Desktop\\DOE_WaterChemistry_HistoricalData")

names <- read.csv("tempNames.csv", header=T)
names$fileNames <- paste(names$DoE_StationCode,"_",names$siteCode, sep="")

compiledData <- read.csv(paste(names$fileNames[1], ".csv", sep=""),
                         header=T, stringsAsFactors=F)#, colClasses=c(dateTime="POSIXct"))
compiledData$siteCode <- names$siteCode[1]
compiledData$DoE_StationCode <- names$DoE_StationCode[1]

for (i in 2:length(names$fileNames)){
  temp <- read.csv(paste(names$fileNames[i], ".csv", sep=""),
                   header=T, stringsAsFactors=F)#, colClasses=c(dateTime="POSIXct"))
  temp$siteCode <- names$siteCode[i]
  temp$DoE_StationCode <- names$DoE_StationCode[i]
  compiledData <- rbind(compiledData,temp)
}

#compose data from 2015-16 ####
more <- read.csv('DOE_WaterChemistry_2015WY_16May2016.csv', stringsAsFactors=F)
library(stringr)

#format dates and times
datevec <- as.character(as.Date(more$Date, format='%m/%d/%Y'))
timevec <- str_pad(paste0(more$Time, ':00'), 8, 'left', 0)
dateTime <- paste(datevec, timevec)

# library(qpcR)
# qpcR:::cbind.na(colnames(compiledData), colnames(more))

#sort and assemble columns
new <- more[,c(9:12,4,4,13:30,4,4,4,4,31:34,4,4,2,1)]
new[,substr(colnames(new), 1, 4) == 'Year'] <- rep(NA, nrow(new))
new <- cbind(dateTime, new)
# qpcR:::cbind.na(colnames(compiledData), colnames(new))
colnames(new) <- colnames(compiledData)

#bind to historical data
compiledData <- rbind(compiledData, new)

#write output
write.csv(compiledData,
          "C:/Users/Mike/git/stream_nuts_DFA/data/chemPhys_data/nutsdata.csv")



