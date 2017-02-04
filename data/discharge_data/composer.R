#discharge data composer
#all data expressed in CFS

rm(list=ls()); cat('\014')
setwd('~/git/puget_sound_rivers_DFA/data/discharge_data')
setwd('C:/Users/Mike/git/stream_nuts_dfa/data/discharge_data')

library(devtools)
# install_github('vlahm/manipulateR')
library(manipulateR)

# compose daily data ####

#read in all files
library(manipulateR)
ultimate_reader(dir_args=list(pattern='.txt', path='./daily/'),
                read_args=list(sep=',', quote='\"', header=FALSE, fill=TRUE,
                               comment.char="", skip=2, stringsAsFactors=FALSE))

#slice, aggregate, sort
files <- ls()
discharge <- data.frame(date=character())
for(i in 1:length(files)){
    print(paste('processing daily -', files[i]))
    x <- eval(parse(text=files[i]))
    x <- x[,1:3]
    colnames(x)<-c('station', 'date', 'discharge')
    x <- aggregate(x$discharge, list(year=substr(x$date,8,11), month=substr(x$date,2,3)), mean, na.rm=TRUE)
    date <- paste(x$year, x$month, '15', sep='-')
    out <- data.frame('date'=date, 'temp'=x$x)[order(substr(date,1,4), substr(date,6,7)),]
    colnames(out)[2] <- files[i]
    row.names(out) <- 1:nrow(out)
    discharge <- merge(discharge, out, by='date', all=TRUE)
}

#re-sort
discharge$date <- as.Date(discharge$date)
discharge <- discharge[order(discharge$date),]
# discharge <- discharge[order(substr(discharge$date,1,4), substr(discharge$date,6,7)),]

#clean up
discharge <- discharge[-(1:27),]
rownames(discharge) <- 1:nrow(discharge)
rm(list=ls()[which(ls() != 'discharge')])

#compose 15 minute data ####

ultimate_reader(dir_args=list(pattern='.txt', path='./15_min/continuations/'),
                read_args=list(sep='\t', quote='\"', header=FALSE, fill=TRUE,
                               comment.char="", skip=29, stringsAsFactors=FALSE))

#slice, aggregate, sort
files <- ls(pattern='[^discharge]')
discharge2 <- data.frame(date=character())
for(i in 1:length(files)){
    print(paste('processing 15_min -', files[i]))
    x <- eval(parse(text=files[i]))
    x <- x[,c(3,5)]
    colnames(x)<-c('date', 'discharge')
    if(any(x$discharge %in% c('Ice', 'Eqp'))){
        x$discharge[x$discharge %in% c('Ice', 'Eqp')] <- NA
        x$discharge <- as.numeric(x$discharge)
    }
    x <- aggregate(x$discharge, list(year=substr(x$date,1,4), month=substr(x$date,6,7)), mean, na.rm=TRUE)
    date <- paste(x$year, x$month, '15', sep='-')
    out <- data.frame('date'=date, 'temp'=x$x)[order(substr(date,1,4), substr(date,6,7)),]
    colnames(out)[2] <- files[i]
    row.names(out) <- 1:nrow(out)
    discharge2 <- merge(discharge2, out, by='date', all=TRUE)
}

#re-sort
discharge2$date <- as.Date(discharge2$date)
# discharge2 <- discharge2[order(substr(discharge2$date,1,4), substr(discharge2$date,6,7)),]
discharge2 <- discharge2[order(discharge2$date),]
discharge2 <- discharge2[-(1:3),]
rownames(discharge2) <- 1:nrow(discharge2)

#merge
disch_2008 <- discharge[1171:1179,]
disch2_2008 <- discharge2[1:9,]
best_2008 <- disch2_2008
best_2008[is.na(best_2008)] <- as.numeric(disch_2008[is.na(best_2008)])
discharge <- rbind(discharge[-(1171:1179),], best_2008, discharge2[-(1:9),])
rownames(discharge) <- 1:nrow(discharge)

#clean up
# rm(list=ls()[!grepl('discharge', ls())])
rm(list=ls()[which(ls() != 'discharge')])

#one more consolidated 15 minute site (not worth it)
# x <- read.table('15_min/stillaguamish_nf_hi.txt', header=FALSE, skip=29, stringsAsFactors=FALSE)
# x <- x[,c(3,6)]
# colnames(x) <- c('date', 'discharge')
# x <- aggregate(x$discharge, list(year=substr(x$date,1,4), month=substr(x$date,6,7)), mean, na.rm=TRUE)
# date <- paste(x$year, x$month, '15', sep='-')
# out <- data.frame('date'=date, 'stillaguamish_nf_hi'=x$x)[order(substr(date,1,4), substr(date,6,7)),]
# row.names(out) <- 1:nrow(out)
# discharge <- merge(discharge, out, by='date', all=TRUE)

#compose 15 minute by-year sites (not worth it. too little data)####

# #stallaguamish south fork
# files <- dir('15_min/by_year/stillaguamish_sf')
# discharge3 <- data.frame(date=character())
# for(i in 1:length(files)){
#     print(files[i])
#     # x <- read.fwf(paste0('15_min/by_year/stillaguamish_sf/2005.TXT'),
#     x <- read.fwf(paste0('15_min/by_year/stillaguamish_sf/', files[i]),
#                   skip=12, widths=c(11,9,15,20), stringsAsFactors=FALSE, strip.white=TRUE)
#     x <- x[,1:3]
#     colnames(x)<-c('date', 'time', 'discharge')
#     # unique(x$discharge[!grepl('0',x$discharge)])
#     # which(x$discharge=='h other station')
#     x$discharge[x$discharge %in% c('Ice', 'Eqp', 'IEWED data', 'ut within 2x', 'h other station',
#                                    'ation across ga', 'ut within 2x', 'another station',
#                                    '- unreliable', 'visional data', 'visional data -', '',
#                                    'xceeded (data w', 'No Data')] <- NA
#     x$discharge <- as.numeric(x$discharge)
#     x <- aggregate(x$discharge, list(year=substr(x$date,7,10), month=substr(x$date,1,2)), mean, na.rm=TRUE)
#     date <- paste(x$year, x$month, '15', sep='-')
#     out <- data.frame('date'=date, 'stillaguamish_sf'=x$x)[order(substr(date,1,4), substr(date,6,7)),]
#     row.names(out) <- 1:nrow(out)
#     discharge3 <- rbind(discharge3, out)
# }
#
# #merge and clean
# discharge <- merge(discharge, discharge3, by='date', all=TRUE)
# rm(list=ls()[which(ls() != 'discharge')])
#
# #stallaguamish mainstem
# files <- dir('15_min/by_year/stillaguamish_main')
# discharge3 <- data.frame(date=character())
# for(i in 1:length(files)){
#     print(files[i])
#     # x <- read.fwf(paste0('15_min/by_year/stillaguamish_sf/2005.TXT'),
#     x <- read.fwf(paste0('15_min/by_year/stillaguamish_main/', files[i]),
#                   skip=12, widths=c(11,9,15,20), stringsAsFactors=FALSE, strip.white=TRUE)
#     x <- x[,1:3]
#     colnames(x)<-c('date', 'time', 'discharge')
#     # unique(x$discharge[!grepl('0',x$discharge)])
#     # which(x$discharge=='h other station')
#     x$discharge[x$discharge %in% c('Ice', 'Eqp', 'IEWED data', 'ut within 2x', 'h other station',
#                                    'ation across ga', 'ut within 2x', 'another station',
#                                    '- unreliable', 'visional data', 'visional data -', '',
#                                    'xceeded (data w', 'No Data', 'ded', 't rating = Fair',
#                                    'easurement - wi', 'by cooperator a', 'less than 1/2x',
#                                    'hecked', 'nt threshold')] <- NA
#     x$discharge <- as.numeric(x$discharge)
#     x <- aggregate(x$discharge, list(year=substr(x$date,7,10), month=substr(x$date,1,2)), mean, na.rm=TRUE)
#     date <- paste(x$year, x$month, '15', sep='-')
#     out <- data.frame('date'=date, 'stillaguamish_main'=x$x)[order(substr(date,1,4), substr(date,6,7)),]
#     row.names(out) <- 1:nrow(out)
#     discharge3 <- rbind(discharge3, out)
# }
#
# #merge and clean
# discharge <- merge(discharge, discharge3, by='date', all=TRUE)
# rm(list=ls()[which(ls() != 'discharge')])

# final fixes ####

#replace NaNs with NAs
# discharge$stillaguamish_sf[is.nan(discharge$stillaguamish_sf)] <- NA
# discharge$stillaguamish_main[is.nan(discharge$stillaguamish_main)] <- NA

#remove garbage rows and Nisqually column
# discharge <- discharge[-(1336:nrow(discharge)),]
discharge <- subset(discharge, select=-nisqually)

#change site names to site IDs
# site_names <- list(elwha='Z', duckabush='I', skokomish='L', deschutes='M', puyallup='J',
#                    green_lo='B', green_hi='N', cedar_lo='A', cedar_hi='O', sammamish_issaquah='ZA',
#                    snohomish='H', skykomish='Q', snoqualmie_lo='R', snoqialmie_hi='P',
#                    stillaguamish_main='G', #stillaguamish_main also equals 'T', but only data for one
#                    stillaguamish_nf_lo='U', stillaguamish_nf_hi='V', stillaguamish_sf='S', skagit_low='F',
#                    skagit_hi='W', samish='E', nooksack_low='C', nooksack_hi='X')

names(discharge)[2:ncol(discharge)] <- c('O','A','M','I','Z','N','B','X','C','J','E','ZA','W','F',
                                         'L','Q','H','P','R','U')#,'V','S','G')

#sort by colname
discharge <- cbind(discharge$date, discharge[,-1][,order(names(discharge)[-1])])
colnames(discharge)[1] <- 'date'

#add rows for missing months
discharge <- merge(discharge,
      data.frame(date=seq.Date(discharge$date[1], discharge$date[nrow(discharge)], 'month')),
      all=T)

#sort by date
discharge <- discharge[order(discharge$date),]
# discharge <- discharge[-(1:27),]
rownames(discharge) <- 1:nrow(discharge)

#write file
write.csv(discharge, 'discharge.csv', row.names=FALSE)
