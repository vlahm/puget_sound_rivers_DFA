rm(list=ls()); cat('\014')
setwd('C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month')

batch_file_reader <- function(dir_args=list(path='./', pattern='.csv'),
                              read_args=list(sep=',', quote="\"", header=TRUE,
                                             fill=TRUE, comment.char=""),
                              merge=FALSE, ...){

    if(!(substr(dir_args$path, nchar(dir_args$path), nchar(dir_args$path)) %in% c('/', '\\'))){
        stop("'path' must include trailing '/' or '\\'")
    }

    #create global objects for contents of all specified files
    files <- do.call('dir', dir_args)
    obj_names <- vector(length=length(files))
    for(i in 1:length(files)){
        obj_names[i] <- substr(files[i], 1, nchar(files[i])-4)
        full_read_args <- append(list(file=paste0(dir_args$path, files[i])), read_args)
        temp <- do.call('read.table', args=full_read_args)
        # temp <- read.table(paste0(dir_args$path, files[i]), sep=',', quote="\"",
        #                    header=TRUE, fill=TRUE, comment.char="")

        if(merge == FALSE){
            assign(obj_names[i], temp, pos='.GlobalEnv')
        } else {
            assign(obj_names[i], temp)
        }
    }

    #merge all objects and output, if specified
    if(merge == TRUE){
        merged <- Reduce(function(x, y) {merge(x, y, ...)},
                         eval(parse(text=paste('list(', paste(obj_names, collapse=','), ')'))))
        return(merged)
    }

    message(paste('files read:', paste(obj_names, collapse=', ')))
}

batch_file_reader()

maxtemp3_4 <- cbind(maxtemp3$date, rowMeans(cbind(maxtemp3$mt, maxtemp4$mt)),
                    rowMeans(cbind(maxtemp3$mt_anom_1900.99, maxtemp4$mt_anom_1900.99)))
maxtemp3_4 <- as.data.frame(maxtemp3_4)
colnames(maxtemp3_4) <- c('date', 'mt', 'mt_anom_1900.99')

meantemp3_4 <- cbind(meantemp3$date, rowMeans(cbind(meantemp3$at, meantemp4$at)),
                    rowMeans(cbind(meantemp3$at_anom_1900.99, meantemp4$at_anom_1900.99)))
meantemp3_4 <- as.data.frame(meantemp3_4)
colnames(meantemp3_4) <- c('date', 'at', 'at_anom_1900.99')

precip3_4 <- cbind(precip3$date, rowMeans(cbind(precip3$pc, precip4$pc)),
                    rowMeans(cbind(precip3$pc_anom_1900.99, precip4$pc_anom_1900.99)))
precip3_4 <- as.data.frame(precip3_4)
colnames(precip3_4) <- c('date', 'pc', 'pc_anom_1900.99')

hdd3_4 <- cbind(hdd3$date, rowMeans(cbind(hdd3$hd, hdd4$hd)),
                    rowMeans(cbind(hdd3$hd_anom_1900.99, hdd4$hd_anom_1900.99)))
hdd3_4 <- as.data.frame(hdd3_4)
colnames(hdd3_4) <- c('date', 'hd', 'hd_anom_1900.99')

hydroDrought3_4 <- cbind(hydroDrought3$date, rowMeans(cbind(hydroDrought3$hdr, hydroDrought4$hdr)),
                    rowMeans(cbind(hydroDrought3$hdr_anom_1900.99, hydroDrought4$hdr_anom_1900.99)))
hydroDrought3_4 <- as.data.frame(hydroDrought3_4)
colnames(hydroDrought3_4) <- c('date', 'hdr', 'hdr_anom_1900.99')

names <- c('hdd3_4', 'precip3_4', 'hydroDrought3_4', 'meantemp3_4', 'maxtemp3_4')
for(i in 1:length(names)){
    write.csv(eval(parse(text=names[i])), row.names=FALSE, file=
                  paste0('C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/', names[i], '.csv'))
}

####
climate3 <- cbind(meantemp3, maxtemp3, precip3, hydroDrought3, hdd3)
climate3 <- climate3[,-seq(4,13,3)]
climate4 <- cbind(meantemp4, maxtemp4, precip4, hydroDrought4, hdd4)
climate4 <- climate4[,-seq(4,13,3)]
climate3_4 <- cbind(meantemp3_4, maxtemp3_4, precip3_4, hydroDrought3_4, hdd3_4)
climate3_4 <- climate3_4[,-seq(4,13,3)]

save(climate3, file="C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/climate3.rda")
save(climate4, file="C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/climate4.rda")
save(climate3_4, file="C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/climate3_4.rda")
