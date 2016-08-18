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

#read in raw csvs
batch_file_reader(dir_args=list(path="C:/Users/Mike/Desktop/Grad/Projects/Thesis/stream_nuts_DFA/data/climate_data/by_month/raw/",
                                pattern='at3.+\\.csv'),
                  merge=FALSE)

#get rid of function object for next step
rm(batch_file_reader)

#rbind all years together
binder <- function(){
    filenames <- ls(name='.GlobalEnv')[ls(name='.GlobalEnv') != 'binder']
    x <- data.frame()
    for(i in 1:length(filenames)){
        y <- eval(parse(text=filenames[i]), envir=.GlobalEnv)
        yy <<- y
        y <- y[-(1:3),-4]
        x <- rbind(x, y)
    }
    return(x)
}
out <- binder()

#sort by date and rename stuff
out <- as.data.frame(apply(out, 2, as.numeric))
metric <- substr(ls(name='.GlobalEnv')[1], 1, 2)
colnames(out) <- c('date', metric, paste(metric, 'anom_1900-99', sep='_'))
sorted <- out[order(out$date),]
rownames(sorted) <- 1:nrow(sorted)

#skip this part if not dealing with temperature dataframes (which are stupidly in fahrenheit)####
base <- sorted$at - sorted$`at_anom_1900-99` #get back to base period for calculating anomalies

base <- (base - 32) * (5/9)#convert base period to C
sorted$at <- (sorted$at - 32) * (5/9) #convert temps to C

sorted$`at_anom_1900-99` <- sorted$at - base #get anomalies in C

#write output csv ####
write.csv(sorted, row.names=FALSE, file=
  'C:/Users/Mike/Desktop/Grad/Projects/Thesis/stream_nuts_DFA/data/climate_data/by_month/meantemp.csv')

