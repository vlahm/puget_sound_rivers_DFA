#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 8/10/2016

#NOTES:

#collapse folds with ALT+O (windows, linux) or CMD+OPT+O (Mac); might have to do it twice

#if testing models, just run through sections 0 through 4 in order. It's set up so that you don't really
#have to change anything except the options in section 1. You will want to look carefully at how much
#the "na_thresh" option filters your dataset. Once you have a model, you can check out various diagnostics
#and basic plots in sections 5 through 5.4. For more serious plots, read the next note.

#if examining a model that was fit using the parallel loop in 03_model_fitting,
#choose the options corresponding to that model here (section 1), and then run sections 0 through
#3.2. Next, set up the path in section 4.2 to read in the appropriate model object.
#Then run sections 5 through 5.4. Finally, save the whole global environment to the figures
#directory using the code in section 4.1. The rest of the instructions are in
#stream_nuts_DFA/manuscript/figures/figscript.R

#to quickly access any of the function definitions, put the cursor on the function name and hit F2

#you should have at least 4 datapoints (observations x streams) for each parameter in model.

#if R crashes when you try to use runDFA,
#use apply(dat_z, 2, function(x) sum(is.na(x))/length(x)) to see if you have any timepoints with
#no data or 1 data point. these timepoints must either be removed or imputed (see section 1.2).

rm(list=ls()); cat('\014') #clear env and console

# 0 - setup ####
setwd('C:/Users/Mike/git/stream_nuts_DFA/data/')
setwd('~/git/puget_sound_rivers_DFA/data')
setwd('Z:/stream_nuts_DFA/data/') #set to data folder
load('chemPhys_data/yys_bymonth.rda')
DISCHARGE <- read.csv('discharge_data/discharge.csv', stringsAsFactors=FALSE, colClasses=c('date'='Date'))
snowmelt <- read.csv('climate_data/snow_data/snowmelt.csv')

#install packages that aren't already installed (see https://github.com/kaskr/adcomp for TMB package)
#imputeTS, RColorBrewer, cluster, fpc
package_list <- c('MARSS','viridis','vegan', 'e1071', 'imputeTS', 'stringr',
                  'foreach', 'doParallel', 'caret', 'Matrix')
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos="http://cran.rstudio.com/")
# if (!require("manipulateR")) {
#     if (!require("devtools")) install.packages('devtools', repos="http://cran.rstudio.com/")
#     library(devtools)
#     install_github('vlahm/manipulateR') #this one is from github
#     detach('package:devtools', unload=TRUE)
# }
# for(i in c(package_list, 'manipulateR')) library(i, character.only=TRUE) #and load them all

#open new plot window unless already open; this is here to prevent issues with plot window size in section 1.3
#it's fine to comment out this block as long as you look out for errors below
if (is.null(dev.list()) == TRUE){
    if(.Platform$OS.type == "windows"){
        windows(record=TRUE, width=16, height=9)
    } else {
        x11(width=16, height=9)
    }
}

# 1 - OPTIONS ####

# response choices: COND FC NH3_N NO2_NO3 OP_DIS OXYGEN PH PRESS SUSSOL TEMP TP_P TURB
# also DISCHARGE (from USGS)
y_choice = 'TEMP'
# cov choices: meantemp meantemp_anom precip precip_anom hydroDrought hydroDrought_anom
# maxtemp maxtemp_anom hdd hdd_anom, snowmelt
#if you're looking at effects by month, or by month over time (see 'design' option below),
#then whichever covariate is listed first here will be the one that can be examined that way.
#any others will still be included, and their average effect sizes (not by month, not over time) can
#still be determined.
#(snowmelt only available 1978-2015. also I haven't actually used snowmelt
#in a model yet, so there could be bugs)
cov_choices = c('meantemp', 'precip')
#region choices: '3' (lowland), '4' (upland), '3_4' (average of 3 and 4, or each separately)
#regions 3 and 4 were discovered to be very similar early on, so most of this script will only work if
#you choose '3_4' here and 'average_regions=TRUE'.
region = '3_4' #code not set up to include snowmelt unless region='3_4' and average_regions=TRUE
#average regions 3 and 4? (if FALSE, sites from each region will be assigned their own climate covariates)
average_regions = TRUE #only used if region = '3_4'
#seasonality model choices: 'fixed_collective', 'fixed_individual', 'fourier'
#collective vs. indiv is no longer relevant. choose either "fixed_..." option to get
#an appropriate fixed-effect cc matrix, or "fourier" to get a simpler one
method = 'fixed_individual'
#which years to include?
startyr = 1978
endyr = 2015
#model params (specific values only relevant for testing, not for parameter optimization loop)
ntrends = 1
#error matrix can either hold the MARSS specifications or the TMB ones ('DE', 'DUE', 'EVCV', 'UNC')
obs_err_var_struc = 'DUE'
#UPDATE: Mark Schueurell no longer scales his response data. scaling forces the variance of
#the D matrix to be small, thus artificially diminishing the impact of the covariates.
#Scaling is still recommended by the community at large, and I didn't see a big difference.
scale = TRUE
na_thresh = 0.55 #exclude sites with >= this proportion of NA values.
#be sure to visit section 3.1, where you can add time interval factors and interaction effects
#to the covariate matrix
#transformations are 'log' and 'none' from here. can also explore 'power' and 'boxcox' in section 3.1
#run function transformables() to see whether your response needs to be transformed.
transform = 'none'
#choose the covariate matrix design here. options are 'just_effect', 'effect_and_seasonality_without_interaction'
#'effect_byMonth', 'effect_byMonth_noSeas', and 'effect_byMonth_acrossTime'.
#(see "designer" function in section 3.1 for details)
design = 'effect_byMonth_noSeas'
#sections = number of intervals to divide the time series into, if examining change over time
#(see "designer" function in section 3.1 for details)
sections <- 5 #an integer. will be ignored if not applicable
#the months to focus on for by-month effect size (1 is jan...)
#if looking at effect_byMonth_acrossTime, including all months will be too expensive
#(see "designer" function in section 3.1 for details)
focal_months <- 1:12 #a vector of integers between 1 and 12. will be ignored if not applicable.

# 1.1 - subset data according to choices, remove problematic columns ####
library(stringr)

# selects the right data
yy <- eval(parse(text=y_choice))

# subset by year and exclude columns with >= na_thresh proportion of NAs
subsetter <- function(yy, start, end, na_thresh=1){

    #extract year subset
    years <- as.numeric(format(yy$date, '%Y'))
    yy <- yy[years >= start & years <= end, ]

    #extract subset of columns with <= na_thresh proportion of NAs
    col_na_prop <- apply(yy, 2, function(x) sum(is.na(x)) / length(x))
    yy <- yy[,col_na_prop <= na_thresh]

    rownames(yy) <- 1:nrow(yy)

    return(yy)
}
yy <- subsetter(yy, start=startyr, end=endyr, na_thresh=na_thresh)

# remove site K: strong groundwater influence prior to 2005
yy = subset(yy, select=-K)

# this shows the proportion of each column that is made up of the most frequent single value.
# For some variables (e.g. SUSSOL) there are tons of identical readings due to low precision.
# Any column with > 0.1 representation by a single value should be examined and modified/removed,
# especially if the repeated value skews the distribution. this skew may remain even after log
# transforming, which can prevent convergence and/or result in a wonky model.
(max_repeats <- apply(yy[,-1], 2, function(i) table(i)[1]/length(i[!is.na(i)])))
(screwy <-  max_repeats[which(max_repeats > 0.05)]) #I'm concerned about anything over 0.05

#get lists of all integer and non-integer vals in the frame (useful for exploring the next two chunks)
# int = numeric()
# non = numeric()
# for(i in 1:nrow(yy)){
#     for(j in 2:ncol(yy)){
#         if(!is.na(yy[i,j])){
#             hasdecimal <- grepl('\\.', yy[i,j])
#             if(!hasdecimal) int = append(int, yy[i,j])
#             if(hasdecimal) non = append(non, yy[i,j])
#         }
#     }
# }

#here I'm jittering those values (all 1s) so that they vary between 0.5 and 1.5
#I don't want to drop them because they represent almost all of my upland sites.
#Site L also contains six 0s, the only ones in the whole set, so I'm treating them as 1s.
#scratch that. Adding artificial precision to all values that don't already have a tens place
#to smooth out the steps in the density function.
# if(y_choice == 'SUSSOL'){
#     for(i in 1:nrow(yy)){
#         # for(j in names(screwy)){ #uncomment these and comment the ones with '#*' to jitter only the worst
#         for(j in 2:ncol(yy)){
#             # yy[i,j] <- round(ifelse(yy[i,j] %in% 0:1, runif(1, 0.5, 1.5), yy[i,j]), 1)
#             if(!(is.na(yy[i,j]))){                                              #*
#                 if(yy[i,j] == 0) yy[i,j] <- 1                                   #*
#                 hasdecimal <- grepl('\\.', yy[i,j])                             #*
#                 if(!hasdecimal){                                                #*
#                     yy[i,j] <- round(runif(1, yy[i,j]-0.5, yy[i,j]+0.5), 1)     #*
#                 }                                                               #*
#             }                                                                   #*
#         }
#     }
# }

#TURB has a different problem. All months preceding June 1989 were measured at integer precision,
#so the same repeated-value skew problem arises, though it isn't so obvious by looking at the
#overall proportion of duplicates. Here I'm jittering everything pre-June-1989.
#scratch that. Adding artificial precision to all values that don't already have a tens place
#to smooth out the steps in the density function.
# if(y_choice == 'TURB'){
#     for(i in 1:nrow(yy)){
#         # for(i in 1:137){ #uncomment these and comment the ones with '#*' to jitter only pre-june-1989
#         for(j in 2:ncol(yy)){
#             # yy[i,j] <- ifelse(is.na(yy[i,j]), NA, round(runif(1, yy[i,j]-0.5, yy[i,j]+0.5), 1))
#             if(!(is.na(yy[i,j]))){                                              #*
#                 hasdecimal <- grepl('\\.', yy[i,j])                             #*
#                 if(!hasdecimal){                                                #*
#                     yy[i,j] <- round(runif(1, yy[i,j]-0.5, yy[i,j]+0.5), 1)     #*
#                 }                                                               #*
#             }                                                                   #*
#         }
#     }
# }

# subset by region
if(region == '3'){
    yy <- yy[,colnames(yy) %in% c('date','J','K','A','B','ZA','Q','H','T','G','F','E','C','R')]
} else {
    if(region == '4'){
        yy <- yy[,colnames(yy) %in% c('date','Z','I','L','M','N','O','P','S','U','X','V','W')]
    }
}

#convert dates to integers
months <- as.numeric(format(yy[,1], '%m'))
years <- as.numeric(format(yy[,1], '%Y'))
int_dates <- (years - startyr) * 12 + months

#create objects for response variables and covariates
obs_ts <- yy[,-1]

covdict <- list('meantemp'='at','meantemp_anom'='at_anom_1900.99','maxtemp'='mt',
                'maxtemp_anom'='mt_anom_1900.99','precip'='pc','precip_anom'='pc_anom_1900.99',
                'hydroDrought'='hdr','hydroDrought_anom'='hdr_anom_1900.99','hdd'='hd',
                'hdd_anom'='hd_anom_1900.99', 'snowmelt'='snowmelt')

if(region=='3_4' & average_regions==FALSE){
    # load(file='C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/climate3.rda')
    # load(file='C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/climate4.rda')
    load(file="climate_data/by_month/climate3.rda")
    load(file="climate_data/by_month/climate4.rda")
    covs3 <- as.matrix(climate3[substr(climate3$date,1,4) >= startyr &
                                    substr(climate3$date,1,4) <= endyr,
                                colnames(climate3) %in%
                                    covdict[names(covdict) %in% cov_choices]])
    covs4 <- as.matrix(climate4[substr(climate4$date,1,4) >= startyr &
                                    substr(climate4$date,1,4) <= endyr,
                                colnames(climate4) %in%
                                    covdict[names(covdict) %in% cov_choices]])
} else {
    #     load(file=paste0('C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_month/climate',
    #         region, '.rda'))
    load(file=paste0('climate_data/by_month/climate',
                     region, '.rda'))
    covs <- eval(parse(text=ls()[grep('climate.{1,3}', ls())]))

    #add snowmelt
    snowmelt <- read.csv('climate_data/snow_data/snowmelt.csv')
    covs <- merge(covs, snowmelt, by='date', all=TRUE)

    #subset by year and covariate choices; keep order specified in cov_choices
    temp_covs <- covs[substr(covs$date,1,4) >= startyr & substr(covs$date,1,4) <= endyr,]
    covs <- matrix(NA, nrow=nrow(temp_covs), ncol=length(cov_choices))
    for(i in 1:length(cov_choices)){
        covs[,i] <- as.matrix(temp_covs[,colnames(temp_covs) == covdict[names(covdict) == cov_choices[i]]])
        # covs <- as.matrix(covs[substr(covs$date,1,4) >= startyr &
        #                            substr(covs$date,1,4) <= endyr,
        #                        colnames(covs) %in%
        #                            covdict[names(covdict) %in% cov_choices]])
    }
    colnames(covs) <- cov_choices
}

# 1.2 - perform minimal interpolation for zero-data months (these will cause TMB to crash) ####
library(imputeTS)

#locate rows where there are no data
(emptyrows = unname(which(rowSums(obs_ts, na.rm=TRUE)==0)))

#check out the yy dataframe around those points. You may have to impute each individually if
#you dont want to fill in too many non-problematic NAs. 'start' is the month corresponding
#to the first observation that will be incorporated in the imputation. This is how i interpolated
#the three zero-data months for SUSSOL
if(y_choice %in% c('SUSSOL', 'TURB')){
    yts <- ts(obs_ts[70:189,c('J','L','M')], start=10, frequency=12)
    obs_ts[70:189,c('J','L','M')] <- data.frame(round(na.seasplit(yts, 'interpolation')))
}

#make sure it worked as expected
# defpar <- par(mfrow=c(3,1))
# for(i in 1:ncol(yts)){
#     plot(ts(obs_ts[70:189,c('J','L','M')], start=10, frequency=12)[,i], col='purple')
#     lines(yts[,i], col='orange')
# }
# par(defpar)

# 1.3 - (transform), center, (scale) ####
library(caret); library(e1071)

#visualize which responses should be transformed to normal
transformables <- function(){
    par(mfrow=c(4,4))
    for(i in c('COND', 'FC', 'NH3_N', 'NO2_NO3', 'OP_DIS', 'OXYGEN',
               'PH', 'PRESS', 'SUSSOL', 'TEMP', 'TP_P', 'TURB', 'DISCHARGE')){
        x <- as.numeric(as.matrix(eval(parse(text=i))[,-1]))
        plot(density(x, na.rm=TRUE), main=i)
        qqnorm(x, main=paste(i, 'qqnorm'))
        qqline(x, col='red', lwd=2)
    }
}
# transformables()

# Transform nonnormal responses (all but OXYGEN, PRESS, PH, TEMP will be transformed)
# backtransformation of boxcox and power does not work as intended, so our only current option
# is to log transform and then report effect size as change in log(response) per change in covariate
# also center and (scale) all responses and covariates (check inside function for details)
transformer <- function(data, transform, exp=NA, scale, plot=FALSE){
    #plot=T to see the effect of transforming
    #transform can be 'boxcox' or 'power'. if 'power', must specify an exp

    obs_ts2 <- data
    lambdas <- NULL

    if(!(y_choice %in% c('OXYGEN','PRESS','PH','TEMP'))){
        if(transform=='boxcox'){
            pre <- preProcess(as.matrix(obs_ts2), method=c('BoxCox'), fudge=0.01)
            obs_ts2 <- predict(pre, obs_ts2)

            lambdas <- sapply(pre$bc, function(x) x$lambda)
            unchanged <- setdiff(names(obs_ts2), names(lambdas))
            if(length(unchanged) > 0){
                temp <- c(lambdas, rep(NA, length(unchanged)))
                names(temp)[(length(lambdas)+1):(length(lambdas)+length(unchanged))] <-
                    unchanged #this will harmlessly error if nothing was left unchanged
                lambdas <- temp[order(names(temp))]
            }
        } else {
            if(transform=='power'){
                obs_ts2 <- obs_ts2^exp
            } else {
                if(transform=='log'){
                    obs_ts2[obs_ts2==0] <- 0.01 #zeros become -Inf when logged (only a few 0s in dataset)
                    obs_ts2 <- log(obs_ts2)
                }
            }
        }
    }

    sds <- apply(obs_ts2, 2, sd, na.rm=TRUE)
    means <- apply(obs_ts2, 2, mean, na.rm=TRUE)
    scaled <- scale(obs_ts2, scale=scale)
    # backtransformed <- (lambdas * scaled * sd(obs_ts2) + 1) ^ (1/lambdas)

    if(plot){
        par(mfrow=c(4,3))
        for(i in 1:ncol(data)){
            plot(density(data[,i], na.rm=TRUE), main=paste0('raw data (', colnames(data)[i], ')'))
            plot(density(scaled[,i], na.rm=TRUE), main='transformed')
            qqnorm(scaled[,i], main='qqnorm transformed')
            qqline(scaled[,i], col='red', lwd=2)
        }
    }

    out <- list(trans=scaled, sds=sds, lambdas=lambdas, means=means)
    return(out)
}
trans <- transformer(obs_ts, transform=transform, exp=NA, scale=scale, plot=F)

dat_z <- t(trans$trans)
# mean(dat_z[1,], na.rm=T); sd(dat_z[1,], na.rm=T) #verify

#scale and center covariate data
covs_z <- t(scale(as.matrix(covs)))

# 2 - plot response variable time series (unintelligible by month, deprecated) ####
# library(viridis)
# series_plotter <- function(){
#     colors1 <- viridis(ncol(yy)-1, end=1)
#     ymin <- min(yy[,-1], na.rm=TRUE)
#     ymax <- max(yy[,-1], na.rm=TRUE)
#
#     plot(int_dates, yy[,2], type='l', col=colors1[1], ylim=c(ymin,ymax), xlab='time', ylab=paste(y_choice), xaxt='n')
#     for(i in 3:ncol(yy)){
#         lines(int_dates, yy[,i], type='l', col=colors1[i-1])
#     }
#     axis(side=1, at=int_dates[c(T,F)], labels=int_dates[c(T,F)])
# }
# series_plotter()

# 3 - Set up input matrices to MARSS function call (automated for TMB only; see next line)####
#for TMB, only mm and cc (and the composite cov_and_seas) are used,
#and these are determined in the CHOICES section,
#so no need to touch this stuff. if experimenting with MARSS, some of these are relevant.
#note though, that MARSS has been essentially abandoned past section 4

# state equation params
mm <- ntrends #number of hidden processes (trends)
BB <- "identity" #'BB' is identity: 1's along the diagonal & 0's elsewhere (not part of DFA)
uu <- "zero" # 'uu' is a column vector of 0's (not part of DFA)
CCgen <- function(meth=method){
    if(meth %in% c('fixed_collective', 'fourier_collective')){
        out <- matrix(month.abb, mm, 12, byrow=TRUE)
        rownames(out) <- NULL
    } else {
        if(meth %in% c('fixed_individual', 'fourier_individual')){
            out <- paste0(month.abb, 1)
            if(mm > 1){
                for(i in 2:mm){
                    out <- rbind(out, paste0(month.abb, i))
                }
            }
            rownames(out) <- NULL
            # out <- 'unconstrained'
        }
    }
    return(out)
}
CC <- CCgen() #coeffs for covariates in process equation (cc)
# CC <- "zero" #use zero if not including seasonality
ccgen <- function(meth=method){
    if(meth %in% c('fixed_collective', 'fixed_individual')){
        year_block <- diag(12)
        nyears <- endyr-startyr+1

        cc <- Reduce(function(x,y) {cbind(x,y)},
                     eval(parse(text=paste0('list(',
                                            paste(rep('year_block', nyears,), sep=', ', collapse=', '),
                                            ')'))))
        rownames(cc) <- month.abb
    } else { #fourier method
        cos_t = cos(2 * pi * seq(dim(dat_z)[2]) / 12)
        sin_t = sin(2 * pi * seq(dim(dat_z)[2]) / 12)
        cc = rbind(cos_t,sin_t)
    }
    return(cc)
} #(needs work - set up for season and add climate data?)
cc <- ccgen() #covariates for the state processes (seasonal and climatic)
# cc <- "zero" #use zero if not including seasonality
QQ <- "identity"  # 'QQ' is identity (would usually be tuned if we were doing actual marss)

# observation equation params
nn <- length(names(obs_ts)) # number of obs time series
aa <- "zero" # 'aa' is the offset/scaling (zero allowed because we standardize the data)
#(we want zero because we're not interested in what the offsets actually are)
DD <- 'unconstrained'
# DD <- 'zero' #use zero if not including climate covs
dd <- covs_z
# dd <- 'zero' #use zero if not including climate covs
RR <- obs_err_var_struc # 'RR' is var-cov matrix for obs errors (tune this)
ZZgen <- function(){
    ZZ <- matrix(list(0),nn,mm)
    ZZ[,1] <- paste0("z",names(obs_ts),1)
    if(mm > 1){
        for(i in 2:mm){
            ZZ[i:nn,i] <- paste0("z",names(obs_ts)[-(1:(i-1))],i)
        }
    }
    return(ZZ)
}
ZZ <- ZZgen() # 'ZZ' is loadings matrix (some elements set to zero for identifiability)

#combine seasonality and covariate matrices
cov_and_seas <- rbind(cc,covs_z)

# 3.1 - add time interval effects and covariate:focal_month:interval interactions to dd ####

#look inside function for details
designer <- function(design, sections, focal_months){
    interval_fac = month_fac = season_fac = NULL

    #this function creates an appropriate covariate matrix (d) for multiple possible model
    #designs. The 'sections' and 'focal months' arguments will be ignored if not applicable.
    #design options are:

    #'just_effect', which produces a d matrix containing only the covariates, with no seasonal
    #effects. (note that fixed factor effects are the only approach that can be applied via this
    #function. seasonality could also be modeled as a sine wave (univariate vs. 12-variate), which
    #would allow the estimation of covariate effect sizes AND covariate effect sizes by month,
    #but with less accurate accounting for seasonal effects that are not of interest. Email me
    #if you're interested in this approach.

    if(design=='just_effect'){
        out <- cov_and_seas[-(1:12),, drop=F]
    }

    #'effect_and_seasonality_without_interaction', which will account for seasonal regularities
    #and then estimate average covariate effect size(s) (i.e. not per month, and not across time)

    if(design=='effect_and_seasonality_without_interaction'){
        out <- cov_and_seas
    }

    #'effect_byMonth', which will account for seasonal regularities and then estimate the effect
    #of the FIRST covariate listed in cov_choices (section 1) for specified months
    #(supplied via "focal months" argument; can be a subset
    #like c(5,6,7,8) or all months: 1:12). Any other covariates will be included as individual rows and will
    #serve to soak up additional variation. Average, non-time-or-month-specific effect sizes can be examined for
    #these. This design will NOT estimate the average effect size for the first covariate listed.

    if(design %in% c('effect_byMonth','effect_byMonth_noSeas','effect_byMonth_acrossTime')){
        month_fac <- rep(0, 12)
        month_fac[focal_months] <- month.abb[focal_months]
        month_fac <- factor(rep(month_fac, length.out=ncol(covs_z)))
    }

    if(design=='effect_byMonth'){
        interactions <- model.matrix( ~ t(covs_z)[,1]:month_fac - 1)
        # interactions <- interactions[,c(5,4,8,1,9,7,6,2,12,11,10,3)] #sort by month instead of alphabetically
        out <- rbind(cov_and_seas[-13,], t(interactions))
    }

    #if you don't want to include a fixed factor to soak up additional unknown seasonal variation,
    #specify 'effect_byMonth_noSeas'. This will still include covariate effects and will split
    #the effect of the first covariate listed across the focal_months.

    if(design=='effect_byMonth_noSeas'){
        interactions <- model.matrix( ~ t(covs_z)[,1]:month_fac - 1)
        # interactions <- interactions[,c(5,4,8,1,9,7,6,2,12,11,10,3)] #sort by month instead of alphabetically
        out <- rbind(cov_and_seas[-(1:13),], t(interactions))
    }

    #finally, there's 'effect_byMonth_acrossTime', which will capture change in a covariate
    #per month during specified time intervals (to do this for multiple covariates will
    #almost certainly be too expensive). the intervals are supplied via the 'sections'
    #argument. The number of sections you specify will determine how many even (as even as
    #possible given the length of your series) intervals the dataset will be chopped into.
    #for each interval, you'll get the effect size of the covariate for each focal month.
    #WARNING: this setup assumes your time series begins in January (it will also be safe
    #to ensure that it ends in december, though this may not be a requirement. untested.).
    #as with 'effect_byMonth', only the first covariate listed in cov_choices will be structured
    #for by-month-across-time examination.

    if(design=='effect_byMonth_acrossTime'){
        section_length <- ncol(covs_z) / sections

        if(!(is.integer(section_length))){
            values_to_add <- ncol(dat_z) - sections*floor(section_length)
            message(paste('NOTICE: fractional number of months in interval length. first',values_to_add,
                          'section(s) will be one observation longer'))
            vals <- rep(floor(section_length), sections)
            for(i in 1:values_to_add) {vals[i] <- vals[i]+1}
        }

        interval_fac <- factor(sort(rep(1:sections, times=vals)))

        interactions <- model.matrix( ~ t(covs_z)[,1]:month_fac:interval_fac - 1)
        out <- rbind(cov_and_seas[-13,], t(interactions))

    }

    return(out)
}

#regardless of design, I'm keeping the name "cov_and_seas" just so I don't break stuff
cov_and_seas <- designer(design=design, sections=sections, focal_months=focal_months)

# 3.2 - see if region 3 and 4 climate variables differ ####
#feel free to uncomment and experiment with this stuff just to verify that
#the regions barely differ in climate

#visualize
# par(mfrow=c(3,2))
# for(i in 1:ncol(covs3)){
#     plot(scale(covs3[,i]), scale(covs4[,i]), main=colnames(covs3)[i])
# }
# par(mfrow=c(3,2))
# for(i in 1:ncol(covs3)){
#     plot(1:nrow(covs3), scale(covs3[,i]), main=colnames(covs3)[i], col='red', type='l')
#     lines(1:nrow(covs3), scale(covs4[,i]), col='blue')
# }
# plot(1,1, type='n', axes=F, ann=F)
# legend('center', legend=c('3','4'), col=c('red','blue'), lty=1)
#
# #compare model outputs based on region 3 vs 4 covars (focusing on the most different from above)
# #precip
# dfa3 <- runDFA(obs=dat_z, NumStates=mm, ErrStruc='DE', EstCovar=TRUE, Covars=rbind(cc,covs3_z[3,]))
# dfa4 <- runDFA(obs=dat_z, NumStates=mm, ErrStruc='DE', EstCovar=TRUE, Covars=rbind(cc,covs4_z[3,]))
# process_plotter_TMB(dfa3, mm) #this function is defined in section 5.1
# process_plotter_TMB(dfa4, mm)
# #hydrologic drought
# dfa5 <- runDFA(obs=dat_z, NumStates=mm, ErrStruc='DE', EstCovar=TRUE, Covars=rbind(cc,covs3_z[4,]))
# dfa6 <- runDFA(obs=dat_z, NumStates=mm, ErrStruc='DE', EstCovar=TRUE, Covars=rbind(cc,covs4_z[4,]))
# process_plotter_TMB(dfa5, mm)
# process_plotter_TMB(dfa6, mm)

# 4 - run DFA ####
# library(MARSS)

#MARSS full specification
# dfa <- MARSS(y=dat_z, model=list(B=BB, U=uu, C='zero', c='zero', Q=QQ, Z=ZZ, A=aa, D='unconstrained', d=covs_z, R=RR),
#              inits=list(x0=matrix(rep(0,mm),mm,1)),
#              control=list(minit=200, maxit=20000, allow.degen=TRUE), silent=2)
# dfa <- MARSS(y=dat_z, model=list(B=BB, U=uu, C=DD, c=dd, Q=QQ, Z=ZZ, A=aa, D=CC, d=cc, R=RR),
#              inits=dfa$par,
#              control=list(minit=200, maxit=3000), method='BFGS') #can't use BFGS for equalvarcov

# #MARSS form=DFA
# dfa2 <- MARSS(y=dat_z, model=list(m=2, R='diagonal and equal', A='zero'),#, D='unconstrained'),
#               inits=list(x0='zero'), z.score=TRUE, #coef(dfa, type='matrix')$D
#               control=list(minit=1, maxit=100, allow.degen=TRUE), silent=2, form='dfa')#,
# # covariates=cov_and_seas)

param_counter <- function(){
    err_strucs = list('DE'=2, 'DUE'=nrow(dat_z), 'EVCV'=2, 'UNC'=nrow(dat_z)*2)
    R_params = unname(unlist(err_strucs[names(err_strucs)==obs_err_var_struc]))
    nparam = prod(dim(ZZ)) + (nrow(cov_and_seas)*nrow(ZZ)) + R_params
    ndata = length(unlist(which(!is.na(dat_z))))
    message(paste('number of data points (n):', ndata))
    message(paste('number of parameters (p):', nparam))
    message(paste('n/p ratio (ideally > 10):', ndata/nparam))
}
param_counter()

#prepare TMB script
if(obs_err_var_struc %in% c('DE', 'DUE', 'UNC')){
    source('../00_tmb_uncor_Rmat_DE_DUE_UNC.R')
} else if(obs_err_var_struc == 'EVCV'){
    source('../00_tmb_uncor_Rmat_EVCV.R')
} else message("invalid R error structure. must be 'DE', 'DUE', 'UNC' or 'EVCV'")

#run model with TMB
dfa <- runDFA(obs=dat_z, NumStates=mm, ErrStruc=obs_err_var_struc,
              EstCovar=TRUE, Covars=cov_and_seas)

# saveRDS(dfa, '../saved_structures/best_sussol_UNC2mFIXEDatpc19782015.rds')

# #get seasonal effects (I think this chunk is fully obsolete, even for MARSS, but
#i'll leave it here just in case)
# CC_out = coef(dfa, type="matrix")$C
# # The time series of net seasonal effects
# seas = CC_out %*% cc[,1:12]
# rownames(seas) = 1:mm
# colnames(seas) = month.abb
# seas

# 4.1 - save model object or global environment image ####

# save.image('../manuscript/figures/discharge_due_4m_atpc_byMo_allMos.rda')
# save.image('../manuscript/figures/discharge_due_4m_atpc_byMo_acrossTime_may-aug.rda')
# save.image('../manuscript/figures/discharge_due_4m_atpc_byMo_acrossTime_nov-feb.rda')
# save.image('../manuscript/figures/discharge_due_4m_atpc_byMo_acrossTime_MASO.rda')
#
# save.image('../manuscript/figures/temp_due_4m_at_byMo_allMos.rda')
# save.image('../manuscript/figures/temp_due_4m_at_byMo_acrossTime_may-aug.rda')
# save.image('../manuscript/figures/temp_due_4m_at_byMo_acrossTime_nov-feb.rda')
# save.image('../manuscript/figures/temp_due_4m_at_byMo_acrossTime_MASO.rda')

# save.image('../manuscript/figures/discharge_due_5m_atpcsn_byMo_allMos.rda')
# save.image('../manuscript/figures/temp_due_5m_atpcsn_byMo_allMos.rda')

# save.image('../single_trend_exploration/1trend.rda')
# save.image('../single_trend_exploration/1trendNoSeas.rda')
# save.image('../single_trend_exploration/1trendNoSeasNoSnow.rda')
# save.image('../single_trend_exploration/2trendNoSeas.rda')
# save.image('../single_trend_exploration/2trendNoSeasNoSnow.rda')

# 4.2 - or load desired model object ####

#load best temp model and all associated mumbo jumbo
# dfa <- readRDS('../round_11_newApproach_byMo_allMos/model_objects_temp/TEMP_DUE_4m_fixed_factors_at_1978-2015.rds')
# dfa <- readRDS('../round_11_newApproach_byMo_allMos/model_objects_discharge/DISCHARGE_DUE_4m_fixed_factors_atpc_1978-2015.rds')
#
# dfa <- readRDS('../round_12_byMoAcrossTime/may-aug/model_objects_temp/TEMP_DUE_4m_fixed_factors_at_1978-2015.rds')
# dfa <- readRDS('../round_12_byMoAcrossTime/nov-feb/model_objects_temp/TEMP_DUE_4m_fixed_factors_at_1978-2015.rds')
# dfa <- readRDS('../round_12_byMoAcrossTime/marAprSepOct/model_objects_temp/TEMP_DUE_4m_fixed_factors_at_1978-2015.rds')
# dfa <- readRDS('../round_12_byMoAcrossTime/may-aug/model_objects_discharge/DISCHARGE_DUE_4m_fixed_factors_atpc_1978-2015.rds')
# dfa <- readRDS('../round_12_byMoAcrossTime/nov-feb/model_objects_discharge/DISCHARGE_DUE_4m_fixed_factors_atpc_1978-2015.rds')
# dfa <- readRDS('../round_12_byMoAcrossTime/marAprSepOct/model_objects_discharge/DISCHARGE_DUE_4m_fixed_factors_atpc_1978-2015.rds')

# dfa <- readRDS('../round_13_byMo_allMos_scale_hiM/model_objects_discharge/DISCHARGE_DUE_5m_fixed_factors_atpcsn_1978-2015.rds')
# dfa <- readRDS('../round_13_byMo_allMos_scale_hiM/model_objects_temp/TEMP_DUE_5m_fixed_factors_atpcsn_1978-2015.rds')

# cov_and_seas <- readRDS('../saved_structures/fixed_at.rds')
# cc <- readRDS('../saved_structures/fixed.rds')
# trans <- readRDS('../saved_structures/temp_trans.rds')
# obs_ts <- readRDS('../saved_structures/temp.rds')
# covs <- readRDS('../saved_structures/at.rds')

# 5 - plot estimated state processes, loadings, and model fits (MARSS) ####

# varimax rotation to get Z loadings
# Z_est <- coef(dfa, type="matrix")$Z # get the estimated ZZ
# H_inv <- varimax(Z_est)$rotmat # get the inverse of the rotation matrix
# Z_rot = Z_est %*% H_inv # rotate factor loadings
# proc_rot = solve(H_inv) %*% dfa$states # rotate processes
#
# # plot hidden processes
# process_plotter <- function(){
#     par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(mm, 1))
#     xlbl <- int_dates
#     y_ts <- int_dates
#     ylm <- c(-1,1)*max(abs(proc_rot))
#     for(i in 1:mm) {
#         plot(y_ts,proc_rot[i,], type="n", bty="L",
#              ylim=ylm, xlab="", ylab="", xaxt="n")
#         abline(h=0, col="gray")
#         lines(y_ts,proc_rot[i,], lwd=2)
#         # lines(y_ts,proc_rot[i,] * seas[i,], lwd=2, col='green')
#         mtext(paste("Process",i), side=3, line=0.5)
#         axis(1, at=xlbl, labels=xlbl, cex.axis=0.8)
#     }
# }
# process_plotter()
#
# # plot loadings
# loading_plotter <- function(){
#     par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(mm, 1))
#     ylbl <- names(obs_ts)
#     clr <- viridis(nn) #colors may not line up with series plots in section 2
#     ylm <- c(-1,1)*max(abs(proc_rot))
#     minZ <- 0
#     ylm <- c(-1,1)*max(abs(Z_rot))
#     for(i in 1:mm) {
#         plot(c(1:nn)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
#              lwd=2, xlab="", ylab="", xaxt="n", ylim=ylm, xlim=c(0.5,nn+0.5), col=clr)
#         for(j in 1:nn) {
#             if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
#             if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
#             abline(h=0, lwd=1.5, col="gray")
#         }
#         mtext(paste("Factor loadings on process",i),side=3,line=0.5)
#     }
# }
# loading_plotter()
#
# # get model fits & CI's
# get_DFA_fits <- function(MLEobj,alpha=0.05) {
#     ## empty list for results
#     fits <- list()
#     ## extra stuff for var() calcs
#     Ey <- MARSS:::MARSShatyt(MLEobj)
#     ## model params
#     ZZ <- coef(MLEobj, type="matrix")$Z
#     ## number of obs ts
#     nn <- dim(Ey$ytT)[1]
#     ## number of time steps
#     TT <- dim(Ey$ytT)[2]
#     ## get the inverse of the rotation matrix
#     H_inv <- varimax(ZZ)$rotmat
#     ## model expectation
#     fits$ex <- ZZ %*% H_inv %*% MLEobj$states
#     ## Var in model fits
#     VtT <- MARSSkfss(MLEobj)$VtT
#     VV <- NULL
#     for(tt in 1:TT) {
#         RZVZ <- coef(MLEobj, type="matrix")$R - ZZ%*%VtT[,,tt]%*%t(ZZ)
#         SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop=FALSE] %*% t(MLEobj$states[,tt,drop=FALSE])
#         VV <- cbind(VV,diag(RZVZ + SS%*%t(ZZ) + ZZ%*%t(SS)))
#     }
#     SE <- sqrt(VV)
#     ## upper & lower (1-alpha)% CI
#     fits$up <- qnorm(1-alpha/2)*SE + fits$ex
#     fits$lo <- qnorm(alpha/2)*SE + fits$ex
#     return(fits)
# }
# mod_fit <- get_DFA_fits(dfa)
#
# # plot fits
# fits_plotter <- function(){
#     ylbl <- names(obs_ts)
#     xlbl = y_ts = 1:444
#     # par(mfrow=c(1,1), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
#     par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
#     ymin <- min(dat_z, na.rm=TRUE)
#     ymax <- max(dat_z, na.rm=TRUE)
#     for(i in 1:nn) {
#         lo <- mod_fit$lo[i,]
#         mn <- mod_fit$ex[i,]
#         up <- mod_fit$up[i,]
#         plot(y_ts,mn,xlab="",ylab=ylbl[i],xaxt="n",type="n", cex.lab=1.2,
#              ylim=c(ymin,ymax))
#         axis(1, at=xlbl, labels=xlbl, cex.axis=1)
#         points(y_ts,dat_z[i,], pch=16, col="darkblue")
#         lines(y_ts, up, col="darkgray")
#         lines(y_ts, mn, col="black", lwd=2)
#         lines(y_ts, lo, col="darkgray")
#     }
# }
# fits_plotter()

# 5.1 - plot processes, loadings, fits, residuals, ACF, CCF, and get R^2 (TMB) ####
library(viridis)

process_plotter_TMB <- function(dfa_obj, ntrends, chunk=NULL){
    #supply an integer to "chunk" in order to see the trend for the corresponding 4-year chunk of the time series
    #i.e. chunk=1 would limit the x axis range to just the first four years.
    xlm = c(1,int_dates[length(int_dates)])
    if(!is.null(chunk)) xlm=c(chunk*48-47,chunk*48)
    if(ntrends<=4){
        par(mai=c(0.5,0.5,0.5,0.1), oma=c(0,3,0,0), mfrow=c(ntrends, 1))
    } else par(mai=c(0,0,0,0), oma=c(0,3,0,0), mfrow=c(ntrends, 1))
    xlbl <- int_dates
    y_ts <- int_dates
    ylm <- c(-1,1)*max(abs(dfa_obj$Estimates$u))
    for(i in 1:ntrends){
        plot(y_ts,dfa_obj$Estimates$u[i,], type="n", bty="L",
             ylim=ylm, xlab='', xaxt="n", ylab='', xaxs='i', xlim=xlm)
        abline(h=0, col="gray")
        lines(y_ts,dfa_obj$Estimates$u[i,], lwd=2)
        axis(2); axis(1, at=seq(1,xlm[2],12)); abline(v=seq(1,xlm[2],12), col='gray', lty=2)
        # mtext(paste("Process",i), side=3, line=-2)
        mtext(paste("Process",i), side=3, line=-2)
        xlbl = xlbl*c(rep(0,11),1)
        # if(i==2){
        # axis(1, at=xlbl, labels=xlbl, cex.axis=0.8)
        # mtext('Month Index', side=1, line=2.5)
        # }
    }
    # mtext(expression(paste('Water ', degree, 'C (de-meaned)')),
    # side=2, line=-1.5, outer=TRUE)
}
# pdf('../manuscript/figures/diagnostic_plots/01_trends.pdf', width=7, height=5)
# png('../manuscript/figures/04_processes_and_loadings.png', width=7, height=6, units='in', res=96, type='cairo')
# process_plotter_TMB(dfa, mm, chunk=NULL)

# dev.off()

loading_plotter_TMB <- function(dfa_obj, ntrends){
    if(ntrends<=4){
        par(mai=c(0.5,0.5,0.5,0.1), oma=c(0,3,0,0), mfrow=c(ntrends, 1))
    } else par(mai=c(0,0,0,0), oma=c(0,3,0,0), mfrow=c(ntrends, 1))
    ylbl <- names(obs_ts)
    # clr <- viridis(nn) #colors may not line up with series plots in section 2
    ylm <- c(-1,1)*max(abs(dfa_obj$Estimates$u))
    minZ <- 0
    Z_rot <- dfa_obj$Estimates$Z
    ylm <- c(-1,1)*max(abs(Z_rot))
    for(i in 1:ntrends) {
        plot(c(1:nn)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
             lwd=2, xlab="", ylab="", xaxt="n", ylim=ylm, xlim=c(0.5,nn+0.5),
             bty='L', yaxt='n', col='black')
        axis(2, las=2)
        # mtext(paste("Process",i), side=3, line=-2)
        for(j in 1:nn) {
            # if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
            lab_size=1.2
            if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=lab_size)}
            if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=lab_size)}
            abline(h=0, lwd=1.5, col="gray")
        }
        mtext(paste("Factor loadings on process",i),side=3,line=-2, col='gray30')
        # if(i==2){
        # mtext('Site ID', side=1, line=2.5)
        # }
        # mtext("Factor loadings", side=2, line=-21.7, outer=TRUE)
    }
}
# pdf('../manuscript/figures/diagnostic_plots/02_loadings.pdf', width=7, height=4)
# loading_plotter_TMB(dfa, mm)
# dev.off()

# full_fit <- dfa$Estimates$Z %*% dfa$Estimates$u + dfa$Estimates$D %*% rbind(cc,covs_z)
# identical(full_fit, dfa$Fits)
# hiddenTrendOnly_fit <- dfa_obj$Estimates$Z %*% dfa_obj$Estimates$u

# pdf('../manuscript/figures/05_fits_and_residuals.pdf', width=7, height=6)
# png('../manuscript/figures/05_fits_and_residuals.png', width=7, height=6, units='in', res=96, type='cairo')
fits_plotter_TMB <- function(dfa_obj){
    hiddenTrendOnly_fit <- dfa_obj$Estimates$Z %*% dfa_obj$Estimates$u
    par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
    # par(mfrow=c(2,1), mar=c(.5,.5,.5,.5), oma=c(3,3,0,0))
    for(i in 1:ncol(obs_ts)){
        # for(i in 3){
        plot(dfa_obj$Fits[i,], type='n',
             ylim=c(min(dat_z[i,], na.rm=TRUE), max(dat_z[i,], na.rm=TRUE)),
             # col='gray40', xaxt='n', ylim=c(-10,10))
             ylab=rownames(dat_z)[i], xlab='month_index')
        points(dat_z[i,], col='blue', pch=20, cex=1)
        lines(dfa_obj$Fits[i,], lwd=1, col='gray40')
        lines(hiddenTrendOnly_fit[i,], col='black', lwd=1)
    }
    # mtext(expression(paste('Water ', degree, 'C (de-meaned)')), 2, 2)
    # mtext('Month Index', 1, 2, T)
}
# pdf('../manuscript/figures/diagnostic_plots/03_fits.pdf', width=7, height=7, onefile=TRUE)
# fits_plotter_TMB(dfa) #black is model fit, green is hidden-trend-only fit, blue is data
# dev.off()

residuals_plotter <- function(dfa_obj){
    par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
    # for(i in 3){
    for(i in 1:ncol(obs_ts)){
        plot(dat_z[i,] - dfa$Fits[i,],
             # ylab='', pch=20, col='blue')
             xlab='Month Index', ylab=rownames(dat_z)[i], pch=20, col='gray40')
        abline(h=0, lwd=2, lty=2, col='black')
    }
    # mtext('Residual Error', 2, 2)
}
# pdf('../manuscript/figures/diagnostic_plots/04_residuals.pdf', width=7, height=7, onefile=TRUE)
# residuals_plotter(dfa)
# dev.off()

ACF_plotter <- function(dfa_obj){
    par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
    # for(i in 3){
    for(i in 1:ncol(obs_ts)){
        acf(dat_z[i,] - dfa$Fits[i,], na.action=na.pass,
            # ylab='', pch=20, col='blue')
            ylab=rownames(dat_z)[i], pch=20)
    }
}
# pdf('../manuscript/figures/diagnostic_plots/05_acf.pdf', width=7, height=7, onefile=TRUE)
# ACF_plotter()
# dev.off()

#check out cross-correlations too, to see if there are lagged relationships between trends, effect sizes, and covariates.
# defpar = par(mfrow=c(2,1))
# ccf(cov_and_seas[1,], as.vector(dfa$Estimates$u[1,])) #13
# ccf(colSums(cov_and_seas[3:14,]), as.vector(dfa$Estimates$u[1,])) #15:26
# ccf(cov_and_seas[1,], as.vector(dfa$Estimates$u[2,])) #13
# ccf(colSums(cov_and_seas[3:14,]), as.vector(dfa$Estimates$u[2,])) #15:26
# par(defpar)

get_R2 <- function(dfa_obj){
    R2 <- rep(NA, nrow(dat_z))
    for(i in 1:nrow(dat_z)){
        mod <- lm(dat_z[i,] ~ dfa_obj$Fits[i,])
        R2[i] <- summary(mod)$r.squared
    }
    return(list(min=min(R2), median=median(R2), max=max(R2)))
}
# get_R2(dfa)

# 5.2 - landscape variable setup ####

#load and preprocess landscape data
land <- read.csv('watershed_data/watershed_data_simp.csv', stringsAsFactors=FALSE)
land$siteCode[land$siteCode == 'AA'] <- 'ZA' #rename sammamish @ bothell sitecode
land <- land[land$siteCode %in% names(obs_ts),] #remove sites not in analysis
land <- land[match(names(obs_ts), land$siteCode),] #sort landscape data by site order in model
rownames(land) <- 1:nrow(land)

#add some synthetic variables to the landscape dataset
WsAreaOver1000 <- read.csv('watershed_data/WsAreaOver1000.csv')
land <- merge(land, WsAreaOver1000, by='siteCode', all.x=TRUE)

WsSlope <- read.csv('watershed_data/slope/slope.csv')
land <- merge(land, WsSlope, by='siteCode', all.x=TRUE)

pcascores <- read.csv('watershed_data/pca_scores.csv')
land <- merge(land,pcascores, by.x='siteCode', by.y='X')

#choose landscape variables of interest
landvars <- c('BFIWs','ElevWs','PctImp2006WsRp100',
              'PctGlacLakeFineWs','PctAlluvCoastWs','PctIce2011Ws',
              'PctCrop2011Ws', 'PctUrbOp2011WsRp100','PctUrbLo2011WsRp100',
              'PctUrbMd2011WsRp100','PctUrbHi2011WsRp100',
              'RdDensWsRp100','RunoffWs','OmWs',
              'RckDepWs','WtDepWs','PermWs','PopDen2010Ws',
              'WsAreaSqKm','WsAreaOver1000','WsSlope','PC1','PC2')
# landvars <- c('ElevWs','PctImp2006WsRp100',
#               'PctGlacLakeFineWs','PctAlluvCoastWs','PctIce2011Ws',
#               'PctCrop2011Ws', 'PctUrbOp2011WsRp100','PctUrbLo2011WsRp100',
#               'PctUrbMd2011WsRp100','PctUrbHi2011WsRp100',
#               'OmWs',
#               'PermWs','PopDen2010Ws',
#               'WsAreaSqKm')

#get the indices of each of these variables in the main land dataframe
landcols <- rep(NA, length(landvars))
for(i in 1:length(landvars)){
    landcols[i] <- which(colnames(land) == landvars[i])
}

#this identifies the landscape vars that correlate best with the response (used below)
best_landvars <- function(response, top){
    #response = loadings or rescaled_effect_sizes; top = top n vars

    best_cols = best_names = best_cors = matrix(NA, nrow=top, ncol=ncol(response))
    for(i in 1:ncol(response)){
        all_cors <- as.vector(cor(response[,i], land[landcols]))
        rank <- order(abs(all_cors), decreasing=T)
        best_cols[,i] <- landcols[head(rank, top)]
        best_names[,i] <- colnames(land)[best_cols[,i]]
        best_cors[,i] <- all_cors[head(rank, top)]
    }
    return(list(cols=best_cols, names=best_names, cors=best_cors))
}

# nstream <- ncol(obs_ts)
# ncov <- ncol(covs)

# 5.3 - effect size regressions (not set up for boxcox/power transformations) ####
# look inside functions for details

#this converts the scaled effect sizes back to their original scales, based on the original
#SDs of the response and covariate(s). it cannot presently account for the effects of
#transformation (and i doubt it can in general). if the model can't converge on untransformed
#data, our only option is to log transform and report the effect sizes as such. see section
#1.3 for more.
#I tried to make this function account for all the different covariate designs (except ..._acrossTime), but
#there are tons of different contingencies to take into account, so I may have missed something.
#if your effect size regressions don't look like what you expected, start troubleshooting here
#(again, I'm happy to help). If you use design='effect_byMonth', this will
#assume you have seasonality method='fixed...'. it might also work with fourier. but you'll need
#seasonality to be included. If you used 'effect_byMonth_noSeas', you're fine.
eff_rescaler <- function(all_cov, seas, scaled=scale){
    #get covariate effect sizes (D coefficients) from model, isolated from seasonal effects
    if(design == 'effect_byMonth'){
        z_eff_1 <- NULL
        if(length(cov_choices) > 1){
            z_eff_1 <- as.matrix(dfa$Estimates$D[,(nrow(seas)+1):(nrow(seas)+length(cov_choices)-1)])
        }
        z_eff_2 <- as.matrix(rowMeans(dfa$Estimates$D[,(nrow(seas)+length(cov_choices)):
                              (nrow(seas)+length(cov_choices)-1+length(focal_months))]))
        z_effect_size <- cbind(z_eff_2, z_eff_1)
    } else {
        if(design == 'effect_byMonth_noSeas'){
            z_eff_1 <- NULL
            if(length(cov_choices) > 1){
                z_eff_1 <- as.matrix(dfa$Estimates$D[,1:(length(cov_choices)-1)])
            }
            z_eff_2 <- as.matrix(rowMeans(dfa$Estimates$D[,length(cov_choices):
                                                              (length(cov_choices)-1+length(focal_months))]))
            z_effect_size <- cbind(z_eff_2, z_eff_1)
        } else {
            if(design == 'effect_byMonth_acrossTime'){
                stop(writeLines(paste0("Not built to produce correct output for eff_regress_plotter\n",
                                "if design='effect_byMonth_acrossTime'.")))
            } else {
                if(nrow(all_cov) > 2){
                    # z_effect_size <- as.matrix(dfa$Estimates$D[,(nrow(seas)+1):ncol(dfa$Estimates$D)])
                    z_effect_size <- as.matrix(dfa$Estimates$D[,(nrow(seas)+1):(nrow(seas)+nrow(covs_z))])
                } else {
                    z_effect_size <- dfa$Estimates$D
                }
            }
        }
    }

    nstream <- nrow(z_effect_size)
    ncov <- ncol(z_effect_size)

    #convert effects back to original scale (units response/units covar)
    #x/sd(response) = 1/sd(covar); x is the effect size scale factor
    rescaled_effect_size <- matrix(NA, nrow=nstream, ncol=ncov)
    for(i in 1:nstream){
        sd_response <- trans$sds[i]
        for(j in 1:ncov){
            sd_covar <- sd(covs[,j], na.rm=TRUE)
            if(!scaled) sd_response = 1
            rescaled_effect_size[,j] <- z_effect_size[,j] * (sd_response/sd_covar)
        }
    }

    return(rescaled_effect_size)
}
rescaled_effect_size <- eff_rescaler(cov_and_seas, cc)

#grab top landscape vars that correlate with effect size, plot and get stats
# eff_best <- best_landvars(rescaled_effect_size, 6)

#look inside function for details (be sure to change the y axis label to 'D log(resp)/D cov'
#if you log transformed the response. (note that this may have been done by default in section 1.3
#if the response was anything other than OXYGEN, TEMP, PRESS, or PH.
#I have not yet built this function to handle multi-covariate models.
#also note that if you used a log transform above, the y-axis here should be
#'log(change_resp)/change_cov'
eff_regress_plotter <- function(mode, var=NA, col_scale='ElevWs'){
    #mode='exploration' is for use within the model fitting loop
    #automatically selects the best correlated landscape vars
    #mode='indiv' is for plotting against individual landscape vars once a model has been selected.
    #if using 'indiv', must select a var name
    #col_scale determines which variable to color the points by
    #green is high, black is low

    pal <- colorRampPalette(c('black', 'green'))
    ncovs = ncol(rescaled_effect_size)

    if(mode == 'exploration'){
        land_ind <- eff_best[[1]]
        top <- nrow(land_ind)
        par(mfrow=c(top/2, 2))
        for(j in 1:ncol(eff_best[[1]])){
            for(i in 1:top){
                cols <- pal(10)[as.numeric(cut(land[,landcols[landvars==col_scale]], breaks=10))]
                mod <- lm(rescaled_effect_size[,j] ~ land[,land_ind[i,j]])
                p <- round(summary(mod)$coefficients[2,4], 2)
                plot(land[,land_ind[i,j]], rescaled_effect_size[,j],
                     xlab=colnames(land)[land_ind[i,j]], ylab='D resp / D cov',
                     main=paste('covar =', cov_choices[j], '- slope p = ', p),
                     col=cols, pch=colnames(trans$trans))
                abline(mod)
            }
        }
    } else {
        if(mode == 'indiv'){
            par(mfrow=c(ncovs,1))
            cols <- pal(10)[as.numeric(cut(land[,landcols[landvars==col_scale]], breaks=10))]
            if(ncovs == 1){
                mod <- lm(rescaled_effect_size ~ land[,landcols[landvars==var]])
                p <- round(summary(mod)$coefficients[2,4], 2)
                plot(land[,landcols[landvars==var]], rescaled_effect_size,
                     xlab=var, ylab='D resp / D cov',
                     main=paste('covar =', cov_choices, '; slope p = ', p),
                     col=cols, pch=colnames(trans$trans))
                abline(mod)
            } else {
                for(i in 1:ncovs){
                    mod <- lm(rescaled_effect_size[,i] ~ land[,landcols[landvars==var]])
                    p <- round(summary(mod)$coefficients[2,4], 2)
                    plot(land[,landcols[landvars==var]], rescaled_effect_size[,i],
                         xlab=var, ylab='D resp / D cov',
                         main=paste('covar =', cov_choices[i], '; slope p = ', p),
                         col=cols, pch=colnames(trans$trans))
                    abline(mod)
                }
            }
        }
    }
}
eff_regress_plotter('indiv', 'PctIce2011Ws', 'ElevWs')
# eff_regress_plotter('exploration', , 'ElevWs')

# 5.4 - process loading regressions (look inside functions for details) ####

loadings <- dfa$Estimates$Z
load_best <- best_landvars(loadings, 6) #landscape vars that cor best with common trend loadings

load_regress_plotter <- function(mmm, mode, var=NA, col_scale='ElevWs'){
    #requires number of hidden trends as input (mmm); annoying, I know
    #mode='exploration' is for use within the model fitting loop
    #automatically selects the best correlated landscape vars
    #mode='indiv' is for plotting against individual landscape vars once a model has been selected
    #if using 'indiv', must select a var name
    #col_scale determines which variable to color the points by
    #green is high, black is low

    loadings <- dfa$Estimates$Z
    pal <- colorRampPalette(c('orange', 'blue'))

    if(mode == 'exploration'){
        land_ind <- load_best[[1]]
        top <- nrow(land_ind)
        par(mfrow=c(top/2, mmm))
        for(j in 1:mmm){
            for(i in 1:top){
                cols <- pal(10)[as.numeric(cut(land[,landcols[landvars==col_scale]], breaks=10))]
                mod <- lm(loadings[,j] ~ land[,land_ind[i,j]])
                p <- round(summary(mod)$coefficients[2,4], 2)
                plot(land[,land_ind[i,j]], loadings[,j],
                     xlab=colnames(land)[land_ind[i,j]], ylab='factor loading on hidden trend',
                     main=paste('hidden trend', j, 'slope p =', p),
                     col=cols, pch=colnames(trans$trans))
                abline(mod)
            }
        }
    } else {
        if(mode == 'indiv'){
            par(mfrow=c(1,1))
            cols <- pal(10)[as.numeric(cut(land[,landcols[landvars==col_scale]], breaks=10))]
            for(j in 1:mmm){
                mod <- lm(loadings[,j] ~ land[,landcols[landvars==var]])
                p <- round(summary(mod)$coefficients[2,4], 2)
                plot(land[,landcols[landvars==var]], loadings[,j],
                     xlab=var, ylab='factor loading on hidden trend',
                     main=paste('hidden trend', j, 'slope p =', p),
                     col=cols, pch=colnames(trans$trans))
                abline(mod)
            }
        }
    }
}
load_regress_plotter(mm, 'indiv', 'ElevWs', 'WsAreaOver1000')
load_regress_plotter(mm, 'exploration', , 'ElevWs')

# 6 - best TEMP model (abandoned, but there's some useful plotting stuff here) ####

# check out all covariate effect plots just to see if there's anything interesting
#that was missed during fitting
for(i in landvars){
    eff_regress_plotter('indiv', i, 'ElevWs') #covariate effect
}

#site K is a common leverage point and is unaffected by air temp. checking to see if it's
#tidally influenced. in the meantime, removing it from analysis, checking top cors.
#UPDATE: K has been removed from the analysis. See section 1.1
# K_ind <- which(land$siteCode=='K')
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis.
best <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, rescaled_effect_size))))))

#plot best cors along with fitted models
defpar <- par(mfrow=c(3,2))
# full_names <- c('mean water table depth', 'base flow index', '% ice 2011', 'mean elevation')
pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(10)[as.numeric(cut(land_sub$ElevWs, breaks=10))]
for(i in 1:6){
    plot(land_sub[,names(best)[1:6][i]], rescaled_effect_size,
         xlab=names(best)[i], ylab=expression(paste(Delta,'water temp /', Delta, 'air temp')),
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(rescaled_effect_size ~ land_sub[,names(best)[1:6][i]])
    abline(mod, col='gray', lty=2)
}
par(defpar)

#why is water table depth such a strong factor? what else is it correlated with?
rev(tail(sort(abs(apply(land[,43:ncol(land)], 2, function(x) cor(land$WtDepWs, x)))), 15))
#it's probably just elevation (note the returned correlations have been abs()'d)

#okay, so stream temp follows the regional air trend depending primarily on
#base flow, glaciation, and elevation

#check out all common trend plots to see if anything was missed during fitting
defpar <- par(mfrow=c(3,3))
for(i in landvars){
    load_regress_plotter(ncol(dfa$Estimates$Z), 'indiv', i, 'ElevWs') #common trend
}
par(defpar)

#get the best ones
best1 <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, dfa$Estimates$Z[,1]))))))
best2 <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, dfa$Estimates$Z[,2]))))))
# best3 <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, dfa$Estimates$Z[,3]))))))

#plot best cors with loadings along with fitted models

# png('../manuscript/figures/02c_loadings_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# pdf('../manuscript/figures/02b_loadings_reg.pdf', width=7, height=6)
pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(10)[as.numeric(cut(land_sub$ElevWs, breaks=10))]
# full_names <- c('Total runoff (mm)', 'Water table depth (cm)', '% Coastal alluvium',
#                 'Base Flow Index', '% Ice', 'Rock Depth (cm)')
# full_names <- c('% Ice', '% Coastal alluvium','Total runoff (mm)', 'Rock Depth (cm)',
#                 'Soil permeability (cm/hr)', '% Cropland')
defpar <- par(mfrow=c(3,2), mar=c(4,.5,.5,.5), oma=c(4,5,0,0))
for(i in 1:6){
    mod <- lm(dfa$Estimates$Z[,1] ~ land_sub[,names(best1)[i]])
    # mod <- lm(dfa$Estimates$Z[,2] ~ land_sub[,names(best2)[i]])
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.05, ' *', '')
    plot(land_sub[,names(best1)[i]], dfa$Estimates$Z[,1],
    # plot(land_sub[,names(best2)[i]], dfa$Estimates$Z[,2],
    xlab=paste0(names(best1)[i], sig),
    # xlab=paste0(full_names[i], sig),
    # main='blue=low elev, red=high elev',
    col=cols, pch=colnames(trans$trans), yaxt='n')
abline(mod, col='gray', lty=2, lwd=2)
if(i %in% c(1,3,5)){
    axis(2)
}
}
mtext('Loadings on common trend 1', 2, outer=TRUE, line=2.5)
# mtext('Loadings on common trend 2', 2, outer=TRUE, line=2.5)
mtext('Watershed variable', 1, outer=TRUE, line=1.5)
par(defpar)
# dev.off()

# full_names <- c('rock depth', 'riparian road density', 'soil permeability', '% ice 2011',
#                 'riparian open space development', 'coastal alluvium')
for(i in 1:6){
    plot(land_sub[,names(best2)[i]], dfa$Estimates$Z[,2],
         xlab=names(best2)[i], ylab='loading on common trend 2',
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(dfa$Estimates$Z[,2] ~ land_sub[,names(best2)[i]])
    abline(mod, col='gray', lty=2)
}

# for(i in 1:6){
#     plot(land_sub[,names(best2)[i]], dfa$Estimates$Z[,3],
#          xlab=names(best2)[i], ylab='loading on common trend 3',
#          main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
#     mod <- lm(dfa$Estimates$Z[,3] ~ land_sub[,names(best2)[i]])
#     abline(mod, col='gray', lty=2)
# }
par(defpar)

#use this code to compare monthly effects
# pdf("C:/Users/Mike/Desktop/with_sitenames.pdf", width=10)
rescaled_seas <- apply(dfa$Estimates$D[,1:12], 2, function(x) x * trans$sds)
defpar <- par(mfrow=c(3,2))
pal <- colorRampPalette(c('blue', 'green'))
cols <- pal(10)[as.numeric(cut(land$BFIWs, breaks=10))]
for(i in 1:12){
    mod <- lm(rescaled_seas[,i] ~ land$PctIce2011Ws, weights=land$watershedA)
    slope <- round(unname(mod$coefficients[2]), 2)
    plot(land$PctIce2011Ws, rescaled_seas[,i], main=paste(month.abb[i], 'slope =', slope),
         ylab=paste(month.abb[i], 'change in water temp'), xlab='% ice',
         ylim=c(min(rescaled_seas), max(rescaled_seas)), col=cols, cex=1,
         pch=colnames(trans$trans))
    abline(mod, col='gray', lty=2, lwd=2)
    abline(h=0, col='red')
}
par(defpar)
# dev.off()

# 6.1 - best SUSSOL model (abandoned, but there's some useful plotting stuff here) ####

# check out all covariate effect plots just to see if there's anything interesting
#that was missed during fitting
for(i in landvars){
    eff_regress_plotter('indiv', i, 'ElevWs') #covarite effect
}

#site K is a common leverage point and is unaffected by air temp. checking to see if it's
#tidally influenced. in the meantime, removing it from analysis, checking top cors.
#UPDATE: K has been removed from the analysis. See section 1.1
# K_ind <- which(land$siteCode=='K')
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis
best <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, rescaled_effect_size))))))

#plot best cors along with fitted models
defpar <- par(mfrow=c(2,2))
full_names <- c('riparian road density', 'rip. dev. open land', 'runoff', 'riparian imp. surfaces')
pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(10)[as.numeric(cut(land_sub$ElevWs, breaks=10))]
for(i in 1:4){
    plot(land_sub[,names(best)[1:4][i]], rescaled_effect_size,
         xlab=full_names[i], ylab=expression(paste('Log',Delta,'sus-sol. /', Delta, 'precip')),
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(rescaled_effect_size ~ land_sub[,names(best)[1:4][i]])
    abline(mod, col='gray', lty=2)
}
par(defpar)

#check out all common trend plots to see if anything was missed during fitting
defpar <- par(mfrow=c(3,3))
for(i in landvars){
    load_regress_plotter(ncol(dfa$Estimates$Z), 'indiv', i, 'ElevWs') #common trend
}
par(defpar)

#get the best ones
best1 <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, dfa$Estimates$Z[,1]))))))
best2 <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, dfa$Estimates$Z[,2]))))))

#plot best cors along with fitted models
defpar <- par(mfrow=c(3,2))

pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(10)[as.numeric(cut(land_sub$ElevWs, breaks=10))]
full_names <- c('WS area', 'Pop. density', 'Rip. urban (Hi)', 'rock depth',
                'glacial fine sed.', 'elevation')
for(i in 1:6){
    plot(land_sub[,names(best1)[i]], dfa$Estimates$Z[,1],
         xlab=full_names[i], ylab='loading on common trend 1',
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(dfa$Estimates$Z[,1] ~ land_sub[,names(best1)[i]])
    abline(mod, col='gray', lty=2)
}

full_names <- c('runoff', 'soil permeability', 'base flow index', 'Rip. road density',
                'water table depth', 'rip. dev. open land')
for(i in 1:6){
    plot(land_sub[,names(best2)[i]], dfa$Estimates$Z[,2],
         xlab=full_names[i], ylab='loading on common trend 2',
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(dfa$Estimates$Z[,2] ~ land_sub[,names(best2)[i]])
    abline(mod, col='gray', lty=2)
}

par(defpar)

#use this code to compare monthly effects
# pdf("C:/Users/Mike/Desktop/with_sitenames.pdf", width=10)
rescaled_seas <- apply(dfa$Estimates$D[,1:12], 2, function(x) x * trans$sds)
defpar <- par(mfrow=c(3,2))
pal <- colorRampPalette(c('blue', 'green'))
cols <- pal(10)[as.numeric(cut(land$BFIWs, breaks=10))]
for(i in 1:12){
    mod <- lm(rescaled_seas[,i] ~ land$PctIce2011Ws, weights=land$watershedA)
    slope <- round(unname(mod$coefficients[2]), 2)
    plot(land$PctIce2011Ws, rescaled_seas[,i], main=paste(month.abb[i], 'slope =', slope),
         ylab=paste(month.abb[i], 'change in sussol'), xlab='% ice',
         ylim=c(min(rescaled_seas), max(rescaled_seas)), col=cols, cex=1,
         pch=colnames(trans$trans))
    abline(mod, col='gray', lty=2, lwd=2)
    abline(h=0, col='red')
}
par(defpar)
# dev.off()

# 6.1 - best TURB model (abandoned, but there's some useful plotting stuff here) ####

# check out all covariate effect plots just to see if there's anything interesting
#that was missed during fitting
for(i in landvars){
    eff_regress_plotter('indiv', i, 'ElevWs') #covarite effect
}

#site K is a common leverage point and is unaffected by air temp. checking to see if it's
#tidally influenced. in the meantime, removing it from analysis, checking top cors
#UPDATE: K has been removed from the analysis. See section 1.1
# K_ind <- which(land$siteCode=='K')
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis
best <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, rescaled_effect_size))))))

#plot best cors along with fitted models
defpar <- par(mfrow=c(2,2))
full_names <- c('riparian road density', 'rip. dev. open land', 'Rip. urb. (Lo)', 'riparian imp. surfaces')
pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(10)[as.numeric(cut(land_sub$ElevWs, breaks=10))]
for(i in 1:4){
    plot(land_sub[,names(best)[1:4][i]], rescaled_effect_size,
         xlab=full_names[i], ylab=expression(paste('Log-',Delta,'turbidity /', Delta, 'precip')),
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(rescaled_effect_size ~ land_sub[,names(best)[1:4][i]])
    abline(mod, col='gray', lty=2)
}
par(defpar)

#check out all common trend plots to see if anything was missed during fitting
defpar <- par(mfrow=c(3,3))
for(i in landvars){
    load_regress_plotter(ncol(dfa$Estimates$Z), 'indiv', i, 'ElevWs') #common trend
}
par(defpar)

#get the best ones
best1 <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, dfa$Estimates$Z[,1]))))))
best2 <- rev(tail(sort(abs(apply(land_sub, 2, function(x) cor(x, dfa$Estimates$Z[,2]))))))

#plot best cors along with fitted models
defpar <- par(mfrow=c(3,2))

pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(10)[as.numeric(cut(land_sub$ElevWs, breaks=10))]
full_names <- c('runoff', 'soil permeability', 'Rip. road density',
                'water table depth', 'base flow index', 'rip. dev. open land')
for(i in 1:6){
    plot(land_sub[,names(best1)[i]], dfa$Estimates$Z[,1],
         xlab=full_names[i], ylab='loading on common trend 1',
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(dfa$Estimates$Z[,1] ~ land_sub[,names(best1)[i]])
    abline(mod, col='gray', lty=2)
}

full_names <- c('% ice', 'Pop. density', 'runoff', 'WS area',
                'glacial fine sed.', 'Rip. urban (Hi)')

for(i in 1:6){
    plot(land_sub[,names(best2)[i]], dfa$Estimates$Z[,2],
         xlab=full_names[i], ylab='loading on common trend 2',
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(dfa$Estimates$Z[,2] ~ land_sub[,names(best2)[i]])
    abline(mod, col='gray', lty=2)
}

par(defpar)

#use this code to compare monthly effects
# pdf("C:/Users/Mike/Desktop/with_sitenames.pdf", width=10)
rescaled_seas <- apply(dfa$Estimates$D[,1:12], 2, function(x) x * trans$sds)
defpar <- par(mfrow=c(3,2))
pal <- colorRampPalette(c('blue', 'green'))
cols <- pal(10)[as.numeric(cut(land$BFIWs, breaks=10))]
for(i in 1:12){i
    mod <- lm(rescaled_seas[,i] ~ land$RdDensWsRp100)
    slope <- round(unname(mod$coefficients[2]), 2)
    plot(land$RdDensWsRp100, rescaled_seas[,i], main=paste(month.abb[i], 'slope =', slope),
         ylab=paste(month.abb[i], 'change in turbidity'), xlab='Rip. Imp. Surfaces',
         ylim=c(min(rescaled_seas), max(rescaled_seas)), col=cols, cex=1,
         pch=colnames(trans$trans))
    abline(mod, col='gray', lty=2, lwd=2)
    abline(h=0, col='red')
}
par(defpar)
# dev.off()
