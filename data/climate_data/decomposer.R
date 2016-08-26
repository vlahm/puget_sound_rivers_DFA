#for decomposing single-month datasets from:
# http://www.ncdc.noaa.gov/cag/time-series/us/45/3/tavg/ytd/12/1895-2016?base_prd=true&firstbaseyear=1901&lastbaseyear=2000




#be sure to change input dir to batch_file_reader as well as output filename
rm(list=ls()); cat('\014')

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

#by month####
#read in raw csvs
batch_file_reader(dir_args=list(path="C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/raw/",
                                pattern='hdr3.+\\.csv'), #only read in one metric and region at a time
                  merge=FALSE)

#get rid of function object for next step
rm(batch_file_reader)

#rbind all years together
binder <- function(){
    filenames <- ls(name='.GlobalEnv')[ls(name='.GlobalEnv') != 'binder']
    x <- data.frame()
    for(i in 1:length(filenames)){
        y <- eval(parse(text=filenames[i]), envir=.GlobalEnv)
        if(substr(filenames[i], 1, 3) == 'hdr'){
            y <- y[-(1:2),-4]
        } else {
            y <- y[-(1:3),-4]
        }
        x <- rbind(x, y)
    }
    return(x)
}
out <- binder()
rm(binder)

#sort by date and rename stuff
out <- as.data.frame(apply(out, 2, as.numeric))
metric <- substr(ls(name='.GlobalEnv')[2], 1, 3) #set this to 1, 3 for hdr
colnames(out) <- c('date', metric, paste(metric, 'anom_1900-99', sep='_'))
sorted <- out[order(out$date),]
rownames(sorted) <- 1:nrow(sorted)


#this part converts temps and anoms from F to C
base <- sorted[,2] - sorted[,3] #get back to base period for calculating anomalies

base <- (base - 32) * (5/9)#convert base period to C
sorted[,2] <- (sorted[,2] - 32) * (5/9) #convert temps to C

sorted[,3] <- sorted[,2] - base #get anomalies in C

#this part converts precip and anoms from in. to cm
base <- sorted[,2] - sorted[,3] #get back to base period for calculating anomalies

base <- base * 2.54 #convert base period to cm
sorted[,2] <- sorted[,2] * 2.54 #convert temps to cm

sorted[,3] <- sorted[,2] - base #get anomalies in cm

#this part converts heating degree days and anoms from F-HDD to C-HDD
base <- sorted[,2] - sorted[,3] #get back to base period for calculating anomalies

base <- base * (5/9) #convert base period to cm
sorted[,2] <- sorted[,2] * (5/9) #convert temps to cm

sorted[,3] <- sorted[,2] - base #get anomalies in cm

#write output csv
write.csv(sorted, row.names=FALSE, file=
              "C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/hydroDrought3.csv") #change name


#by year (obsolete?) ####

#read in raw csvs
batch_file_reader(dir_args=list(path="C:/Users/Mike/Desktop/Grad/Projects/Thesis/stream_nuts_DFA/data/climate_data/by_year/raw/",
                                pattern='\\.csv'),
                  merge=FALSE)
rm(batch_file_reader)

#rbind all years together
binder <- function(){
    filenames <- ls(name='.GlobalEnv')[ls(name='.GlobalEnv') != 'binder']
    for(i in 1:length(filenames)){
        y <- eval(parse(text=filenames[i]), envir=.GlobalEnv)
        y <- y[-(1:3),-4]

        y <- as.data.frame(apply(y, 2, as.numeric))
        if(substr(filenames[i],1,3) == 'hdr'){
            metric <- substr(filenames[i], 1, 3)
        } else {
            metric <- substr(filenames[i], 1, 2)
        }
        colnames(y) <- c('date', metric, paste(metric, 'anom_1900-99', sep='_'))
        rownames(y) <- 1:nrow(y)

        base <- y[,2] - y[,3] #for determining anoms
        if(substr(filenames[i],1,2) %in% c('at', 'mt')){
            base <- (base - 32) * (5/9)
            y[,2] <- (y[,2] - 32) * (5/9)
        } else {
            if(substr(filenames[i],1,2) == 'pc'){
                base <- base * 2.54
                y[,2] <- y[,2] * 2.54
            } else {
                if(substr(filenames[i],1,2) == 'hd' & substr(filenames[i],3,3) != 'r'){
                    base <- base * (5/9)
                    y[,2] <- y[,2] * (5/9)
                }
            }
        }
        y[,3] <- y[,2] - base

        assign(filenames[i], y, pos='.GlobalEnv')
    }
}
binder()
rm(binder)

outnames <- c('meantemp3', 'meantemp4', 'hdd3', 'hdd4', 'hydroDrought3', 'hydroDrought4',
              'maxtemp3', 'maxtemp4', 'precip3', 'precip4')

list <- ls()[!(ls() %in% c('x', 'outnames', 'i'))]
for(i in 1:length(list)){
    x <- eval(parse(text=list[i]))
    write.csv(x, row.names=FALSE, file=paste(
        'C:/Users/Mike/Desktop/Grad/Projects/Thesis/stream_nuts_DFA/data/climate_data/by_year/',
        outnames[i], '.csv'))
}
