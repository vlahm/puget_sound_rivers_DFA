#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 8/10/2016

#NOTEs - should have at least 4 datapoints (observations x streams) for each parameter in model.
#collapse folds with ALT+O (windows, linux) or CMD+OPT+O (Mac); might have to do it twice.
#if R crashes when you try to use runDFA,
#use apply(dat_z, 2, function(x) sum(is.na(x))/length(x)) to see if you have any timepoints with
#no data or 1 data point. these timepoints must either be removed or imputed.
#to quickly access any of the function definitions, put the cursor on the function name and hit F2

rm(list=ls()); cat('\014') #clear env and console

# 0 - setup ####
# setwd('C:/Users/Mike/git/stream_nuts_DFA/data/')
# setwd('~/git/puget_sound_rivers_DFA/data')
# setwd('Z:/stream_nuts_DFA/data/')
setwd("C:/Users/vlahm/Desktop/stream_nuts_DFA/data")
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

# 1 - CHOICES ####

# response choices: COND FC NH3_N NO2_NO3 OP_DIS OXYGEN PH PRESS SUSSOL TEMP TP_P TURB
# also DISCHARGE (from USGS)
y_choice = 'TEMP'
# cov choices: meantemp meantemp_anom precip precip_anom hydroDrought hydroDrought_anom
# maxtemp maxtemp_anom hdd hdd_anom, snowmelt (snowmelt only available 1978-2015)
#also, Ihaven't actually used snowmelt in a model yet, so there could be bugs)
#specify all of the covs that will be used later in the fitting loop
cov_choices = c('meantemp', 'precip', 'snowmelt')
#region choices: '3' (lowland), '4' (upland), '3_4' (average of 3 and 4, or each separately)
region = '3_4'#code not set up to include snowmelt unless region='3_4' and average_regions=TRUE
#average regions 3 and 4? (if FALSE, sites from each region will be assigned their own climate covariates)
average_regions = TRUE #only used if region = '3_4'
#which years to include?
startyr = 1978
endyr = 2015
#model params (just leave these as they are for model fitting. their values wont be used
#but they still need to exist)
method = 'fixed_individual'
ntrends = 1
obs_err_var_struc = 'diagonal and equal'
#UPDATE: Mark Schueurell no longer scales his response data. scaling forces the variance of
#the D matrix to be small, thus artificially diminishing the impact of the covariates.
scale = FALSE
#exclude sites with >= this proportion of NA values.
na_thresh = 0.55
#transformations are 'log' and 'none' from here. can also explore 'power' and 'boxcox' in section 3.1
#run function transformables() to see whether your response needs to be transformed.
transform = 'none'

#be sure to visit section 3.1, where you can add time interval factors and interaction effects
#to the covariate matrix

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
if(!y_choice=='DISCHARGE') yy <- subset(yy, select=-K)

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
#         for(j in 2:ncol(yy)){
#         # for(j in names(screwy)){ #uncomment these and comment the ones with '#*' to jitter only the worst
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
#     # for(i in 1:137){ #uncomment these and comment the ones with '#*' to jitter only pre-june-1989
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

  #subset by year and covariate choices
  covs <- as.matrix(covs[substr(covs$date,1,4) >= startyr &
                           substr(covs$date,1,4) <= endyr,
                         colnames(covs) %in%
                           covdict[names(covdict) %in% cov_choices]])
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

  #it won't transform if the response is any of the four below. this was convenient once, but seems restrictive now.
  #if you want to transform these, go for it.
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

  out <- list(trans=scaled, sds=sds, lambdas=lambdas)
  return(out)
}
trans <- transformer(obs_ts, transform=transform, exp=NA, scale=scale, plot=T)

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

cov_and_seas <- rbind(cc,covs_z)

# 3.1 - add time interval effects and covariate:focal_month:interval interactions to dd ####

#look inside function for details
factorizer <- function(sections, focal_months){

  #create factor of even intervals across time series
  #argument 'sections' is the number of time intervals to divide the time series into
    if(sections!='none'){
        section_length <- ncol(covs_z) / sections

        if(!(is.integer(section_length))){
        values_to_add <- ncol(dat_z) - sections*floor(section_length)
        message(paste('NOTICE: fractional number of months in interval length. first',values_to_add,
                      'section(s) will be one observation longer'))
        vals <- rep(floor(section_length), sections)
        for(i in 1:values_to_add) {vals[i] <- vals[i]+1}
        }

        interval_fac <- factor(sort(rep(1:sections, times=vals)))
    } else interval_fac <- rep(NA, ncol(dat_z))

  #create model matrix for desired interactions (for which months would you like to see
  #the interaction of covariate:interval:month?
  #(i.e. change in covariate per month during each time interval)
  #I'm interested in hydrograph change in november-december and may-june, so I'm specifying their
  #corresponding numbers as the "focal_months"
  #WARNING: this function assumes your time series begins in January
  month_fac <- rep(0, 12)
  month_fac[focal_months] <- month.abb[focal_months]
  month_fac <- factor(rep(month_fac, length.out=ncol(covs_z)))

  return(data.frame(series_interval=interval_fac, focal_months=month_fac))
}
# facs <- factorizer(5, c(11,12,5,6))

#create model matrix for interval effect
# interval_effect <- model.matrix( ~ facs[,1] -1)

#create model matrix for covariate:(interval):focal_month interactions
#bind them to the covariates and seasonality components (keeping the same name to preserve code)
# if(sections != 'none'){
#     interactions <- model.matrix( ~ t(covs_z):facs[,1]:facs[,2] - 1)
#     cov_and_seas <- rbind(cov_and_seas, t(interval_effect), t(interactions))
# } else {
#     interactions <- model.matrix( ~ t(covs_z):facs[,2] - 1)
#     cov_and_seas <- rbind(cov_and_seas, t(interactions))
# }

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

# 5.1 - plot processes, loadings, fits, residuals, ACF, and get R^2 (TMB) ####
library(viridis)

#these functions are used later within the loop

process_plotter_TMB <- function(dfa_obj, ntrends){
  par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(ntrends, 1))
  xlbl <- int_dates
  y_ts <- int_dates
  ylm <- c(-1,1)*max(abs(dfa_obj$Estimates$u))
  for(i in 1:ntrends){
    plot(y_ts,dfa_obj$Estimates$u[i,], type="n", bty="L",
         ylim=ylm, xlab="", ylab="", xaxt="n")
    abline(h=0, col="gray")
    lines(y_ts,dfa_obj$Estimates$u[i,], lwd=2)
    mtext(paste("Process",i), side=3, line=0.5)
    axis(1, at=xlbl, labels=xlbl, cex.axis=0.8)
  }
}
# process_plotter_TMB(dfa, mm)

loading_plotter_TMB <- function(dfa_obj, ntrends){
  par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(ntrends, 1))
  ylbl <- names(obs_ts)
  clr <- viridis(nn) #colors may not line up with series plots in section 2
  ylm <- c(-1,1)*max(abs(dfa_obj$Estimates$u))
  minZ <- 0
  Z_rot <- dfa_obj$Estimates$Z
  ylm <- c(-1,1)*max(abs(Z_rot))
  for(i in 1:ntrends) {
    plot(c(1:nn)[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
         lwd=2, xlab="", ylab="", xaxt="n", ylim=ylm, xlim=c(0.5,nn+0.5), col=clr)
    for(j in 1:nn) {
      if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
      if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
      abline(h=0, lwd=1.5, col="gray")
    }
    mtext(paste("Factor loadings on process",i),side=3,line=0.5)
  }
}
# loading_plotter_TMB(dfa, mm)

# full_fit <- dfa$Estimates$Z %*% dfa$Estimates$u + dfa$Estimates$D %*% rbind(cc,covs_z)
# identical(full_fit, dfa$Fits)
# hiddenTrendOnly_fit <- dfa_obj$Estimates$Z %*% dfa_obj$Estimates$u

fits_plotter_TMB <- function(dfa_obj){
  hiddenTrendOnly_fit <- dfa_obj$Estimates$Z %*% dfa_obj$Estimates$u
  par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
  for(i in 1:ncol(obs_ts)){
    plot(dfa_obj$Fits[i,], type='l', lwd=2,
         ylim=c(min(dat_z[i,], na.rm=TRUE), max(dat_z[i,], na.rm=TRUE)),
         ylab=rownames(dat_z)[i], xlab='day_index')
    lines(hiddenTrendOnly_fit[i,], col='green', lwd=2)
    points(dat_z[i,], col='blue', pch=20, cex=1)
  }
}
# fits_plotter_TMB(dfa) #black is model fit, green is hidden-trend-only fit, blue is data

residuals_plotter <- function(dfa_obj){
  par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
  for(i in 1:ncol(obs_ts)){
    plot(dat_z[i,] - dfa$Fits[i,], ylab=rownames(dat_z)[i], xlab='day_index',
         pch=20)
    abline(h=0, col='gray', lwd=2)
  }
}
# residuals_plotter(dfa)

ACF_plotter <- function(dfa_obj){
  par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
  for(i in 1:ncol(obs_ts)){
    acf(dat_z[i,] - dfa$Fits[i,], na.action=na.pass,
        ylab=rownames(dat_z)[i], pch=20)
  }
}
# ACF_plotter()

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

#add watershed area over 1000m column
WsAreaOver1000 <- read.csv('watershed_data/WsAreaOver1000.csv')
land <- merge(land, WsAreaOver1000, by='siteCode', all.x=TRUE)

#choose landscape variables of interest
landvars <- c('BFIWs','ElevWs','PctImp2006WsRp100',
              'PctGlacLakeFineWs','PctAlluvCoastWs','PctIce2011Ws',
              'PctCrop2011Ws', 'PctUrbOp2011WsRp100','PctUrbLo2011WsRp100',
              'PctUrbMd2011WsRp100','PctUrbHi2011WsRp100',
              'RdDensWsRp100','RunoffWs','OmWs',
              'RckDepWs','WtDepWs','PermWs','PopDen2010Ws',
              'WsAreaSqKm','WsAreaOver1000')

#get the indices of each of these variables in the main land dataframe
landcols <- rep(NA, length(landvars))
for(i in 1:length(landvars)){
  landcols[i] <- which(colnames(land) == landvars[i])
}

#this identifies the landscape vars that correlate best with the response (used in the loop)
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
# these functions are used later within the loop
# look inside functions for details

#this converts the scaled effect sizes back to their original scales, based on the original
#SDs of the response and covariate(s). it cannot presently account for the effects of
#transformation (and i doubt it can in general). if the model can't converge on untransformed
#data, our only option is to log transform and report the effect sizes as such. see section
#1.3 for more.
eff_rescaler <- function(all_cov, seas, scaled=scale){
  #get covariate effect sizes (D coefficients) from model, isolated from seasonal effects
  if(nrow(all_cov) > 2){
    # z_effect_size <- as.matrix(dfa$Estimates$D[,(nrow(seas)+1):ncol(dfa$Estimates$D)])
    z_effect_size <- as.matrix(dfa$Estimates$D[,(nrow(seas)+1):(nrow(seas)+nrow(covs_z))])
  } else {
    z_effect_size <- dfa$Estimates$D
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
# rescaled_effect_size <- eff_rescaler(cov_and_seas, cc)

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

  if(mode == 'exploration'){
    land_ind <- eff_best[[1]]
    top <- nrow(land_ind)
    par(mfrow=c(top/2, 2))
    for(j in 1:ncol(eff_best[[1]])){
      for(i in 1:top){
        cols <- pal(10)[as.numeric(cut(land[,landcols[landvars==col_scale]], breaks=10))]
        plot(land[,land_ind[i,j]], rescaled_effect_size[,j],
             xlab=colnames(land)[land_ind[i,j]], ylab='D resp / D cov',
             main=paste('covar =', cov_choices[j]),
             col=cols, pch=colnames(trans$trans))
      }
    }
  } else {
    if(mode == 'indiv'){
      cols <- pal(10)[as.numeric(cut(land[,landcols[landvars==col_scale]], breaks=10))]
      plot(land[,landcols[landvars==var]], rescaled_effect_size,
           xlab=var, ylab='D resp / D cov',
           main=paste('covar =', cov_choices),
           col=cols, pch=colnames(trans$trans))
    }
  }
}
# eff_regress_plotter('indiv', 'WtDepWs')

# 5.4 - process loading regressions (look inside functions for details) ####
#these functions are used later within the loop

# loadings <- dfa$Estimates$Z
# load_best <- best_landvars(loadings, 6) #landscape vars that cor best with common trend loadings

load_regress_plotter <- function(mmm, mode, var=NA, col_scale='ElevWs'){
  #requires number of hidden trends as input (mmm); annoying, I know
  #mode='exploration' is for use within the model fitting loop
  #automatically selects the best correlated landscape vars
  #mode='indiv' is for plotting against individual landscape vars once a model has been selected
  #if using 'indiv', must select a var name
  #col_scale determines which variable to color the points by
  #green is high, black is low

  loadings <- dfa$Estimates$Z
  pal <- colorRampPalette(c('black', 'green'))

  if(mode == 'exploration'){
    land_ind <- load_best[[1]]
    top <- nrow(land_ind)
    par(mfrow=c(top/2, mmm))
    for(j in 1:mmm){
      for(i in 1:top){
        cols <- pal(10)[as.numeric(cut(land[,landcols[landvars==col_scale]], breaks=10))]
        plot(land[,land_ind[i,j]], loadings[,j],
             xlab=colnames(land)[land_ind[i,j]], ylab='factor loading on hidden trend',
             main=paste('hidden trend', j), col=cols, pch=colnames(trans$trans))
      }
    }
  } else {
    if(mode == 'indiv'){
      cols <- pal(10)[as.numeric(cut(land[,landcols[landvars==col_scale]], breaks=10))]
      for(j in 1:mmm){
        plot(land[,landcols[landvars==var]], loadings[,j],
             xlab=var, ylab='factor loading on hidden trend',
             main=paste('hidden trend', j),
             col=cols, pch=colnames(trans$trans))
      }
    }
  }
}
# load_regress_plotter(2, 'indiv', 'WtDepWs')

# 6 - model fitting/parameter tuning (MARSS - not parallelized) ####
# library(MARSS)
# R_strucs <- c('diagonal and equal', 'diagonal and unequal', 'equalvarcov')
# R_names <- c('DiagAndEq', 'DiagAndUneq', 'EqVarCov')
# ntrends <- 1:4
# model_out <- data.frame()
#
# # ntrends <- 5
# # model_out <- data.frame()
# # R_strucs <- c('diagonal and unequal')
# # R_strucs <- c('unconstrained')
# # RRR='diagonal and unequal'
# # RRR='unconstrained'
# # R_names <- c('DiagAndUneq')
# # R_names <- c('Unconst')
# # mmm=1
#
# for(RRR in R_strucs){
#     for(mmm in ntrends){
#
#         #create loadings matrix
#         mm <- mmm
#         ZZZ <- ZZgen()
#
#         #fit model with EM algorithm
#         print(paste(RRR,mmm))
#
#
#         if(mmm == 1){
#             print(paste('should be 1', mmm))
#
#             dfa <- try(MARSS(y=dat_z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
#                              inits=list(x0=0), silent=FALSE,
#                              control=list(maxit=20000, allow.degen=TRUE)))
#             if(isTRUE(class(dfa)=='try-error')) {next}
#         } else {
#             print(mmm)
#
#             dfa <- try(MARSS(y=dat_z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
#                              inits=list(x0=matrix(rep(0,mmm),mmm,1)), silent=FALSE,
#                              control=list(maxit=20000, allow.degen=TRUE)))
#             if(isTRUE(class(dfa)=='try-error')) {next}
#         }
#
#         #where possible, polish EM estimate with BFGS
#         if(RRR != 'equalvarcov'){
#             print(paste(RRR,mmm,'BFGS'))
#
#             dfa <- try(MARSS(y=dat_z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
#                              inits=dfa$par, silent=FALSE,
#                              # inits=coef(dfa, form='marss'), #alternate form? - see MARSSoptim examples
#                              control=list(dfa$control$maxit), method='BFGS'))
#             if(isTRUE(class(dfa)=='try-error')) {next}
#         }
#
#         #store params, etc. in dataframe
#         model_out <- rbind(model_out, data.frame(R=RRR, m=mmm, LogLik=dfa$logLik,
#                                                  K=dfa$num.params, AICc=dfa$AICc,
#                                                  AIC=dfa$AIC, converged=dfa$convergence,
#                                                  filter=dfa$fun.kf, method=dfa$method,
#                                                  iter=unname(dfa$numIter), n=dfa$samp.size,
#                                                  stringsAsFactors=FALSE))
#
#         #save model object
#         saveRDS(dfa, file=paste0("../stream_nuts_DFA/model_objects/",
#                                  R_names[R_strucs == RRR], '_', mmm, 'm_', startyr, y_choice, '.rds'))
#
#         if(mmm > 1){
#
#             #open plot device
#             pdf(file=paste0("../stream_nuts_DFA/model_outputs/",
#                             R_names[R_strucs == RRR], '_', mmm, 'm_', startyr, y_choice, '.pdf'), onefile=TRUE)
#
#             # varimax rotation to get Z loadings
#             Z_est <- coef(dfa, type="matrix")$Z # get the estimated ZZ
#             H_inv <- varimax(Z_est)$rotmat # get the inverse of the rotation matrix
#             Z_rot = Z_est %*% H_inv # rotate factor loadings
#             proc_rot = solve(H_inv) %*% dfa$states # rotate processes
#
#             #plot hidden processes, loadings, and model fits for each parameter combination
#             process_plotter()
#             loading_plotter()
#
#             mod_fit <- get_DFA_fits(dfa)
#             fits_plotter()
#
#             #close plot device
#             dev.off()
#
#         }
#     }
# }
# #save data frame
# write.csv(model_out, file=paste0("../stream_nuts_DFA/model_objects/",
#                                  'param_tuning_dataframe_', startyr, y_choice, '.csv'))

# 6.1 - model fitting/parameter tuning (TMB) ####
library(doParallel); library(foreach)

#specify system-dependent cluster type; ncores-1 to be used in parallel (no need to modify this)
if(.Platform$OS.type == "windows"){
  cl <- makeCluster(detectCores() - 1, type='PSOCK')
} else {
  cl <- makeCluster(detectCores() - 1, type='FORK')
}
registerDoParallel(cl)

# getDoParWorkers()
# getDoParRegistered() #make sure (this will not return NULL)
# getDoParName() #see what it is
# registerDoSEQ()
# stopCluster(cl)
# stopImplicitCluster()

unregister <- function() {
  env <- foreach:::.foreachGlobals
  rm(list=ls(name=env), pos=env)
}
# unregister() #this resets everything related to parallelization

#create vectors of parameter values to loop through
R_strucs <- c('DUE')#, 'DE') #'UNC'
# R_strucs <- c('EVCV')
#ISSUE: during parallelization, each new thread runs its own separate TMB script, but this original thread (the
#one from which you start the parallelizing loop) also has to run a TMB script, and I don't know how to make it run
#different ones depending on which R_strucs is assigned to it. So, EVCV should be run in its own separate
#parallel loop.
ntrends <- 3:4
seasonality <- list(fixed_factors=ccgen('fixed_individual'))#, #this is calling functions defined above
#fourier=ccgen('fourier'))#, no_seas=NULL)
#subset particular covariates from the full covariate matrix (use rbind to turn vectors into row vectors)
covs_z2 <- covs_z #for name preservation elsewhere
covariates <- list(pc=rbind(covs_z2[2,]),#, sn=rbind(covs_z2[3,]),
                   at=rbind(covs_z2[1,]),
                   # pcsn=covs_z2[2:3,])
                   atpc=covs_z2[1:2,])#, atpcsn=covs_z2[1:3,],
                   # atsn=covs_z2[c(1,3),])
#for the next two, I'm only interested in these values, but you could add loops for them
#just like the other loops below

#number of intervals to divide the time series into, for examining change over time
sections <- 'none' #can be 'none'
#the months to focus on for by-month effect size (1 is jan...)
focal_months <- c(5,6,7,8)
#switch this to TRUE if you're only using covs and seasonal effects
covSeas_only = FALSE

# troubleshooting
# sss = 1 #uncomment to fix seasonality (must also comment **s below)
# RRR = 'UNC' #uncomment to fix error structure (must also comment **R below)
#stock parameters for troubleshooting
# RRR='DE'; mmm=1; cov=1; sss=1
# rm(RRR, mmm, cov, sss)
# rm(cov_and_seas, dfa, all_cov, seas)
#sourcing the appropriate TMB script based on the chosen error structure
if('EVCV' %in% R_strucs) source('../00_tmb_uncor_Rmat_EVCV.R') else source('../00_tmb_uncor_Rmat_DE_DUE_UNC.R')

#this is the model fitting loop. it opens new R processes for each processor core
#and allocates model runs efficiently to each. the outputs are stored in stream_nuts_DFA/model_objects
#and stream_nuts_DFA/model_outputs. The variable "model_out" is a summary data frame of diagnostic stuff
#that will be dealt with next
model_out <-
    foreach(RRR=R_strucs, .combine=rbind) %:% # **R
    foreach(mmm=ntrends, .combine=rbind) %:%
    foreach(sss=1:length(seasonality), .combine=rbind) %:% # **s
    foreach(cov=1:length(covariates), .combine=rbind,
    .packages='viridis') %dopar% {

    #re-run TMB script for each new R process
    if('EVCV' %in% R_strucs) source('../00_tmb_uncor_Rmat_EVCV.R') else
      source('../00_tmb_uncor_Rmat_DE_DUE_UNC.R')

    # print(paste(RRR,mmm,names(seasonality)[sss],names(covariates)[cov]))

    #create covariate + seasonality + interaction + (interval/section) matrix
      if(!covSeas_only){
        #create model matrix for covariate:(interval):focal_month interactions
        #bind them to the covariates and seasonality components (keeping the same name to preserve code)
        facs <- factorizer(sections, focal_months)
        if(sections != 'none'){
          interval_effect <- model.matrix( ~ facs[,1] -1)
          interactions <- model.matrix( ~ t(covs_z):facs[,1]:facs[,2] - 1)
          cov_and_seas <- rbind(seasonality[[sss]], covs_z, t(interval_effect), t(interactions))
        } else {
          interactions <- model.matrix( ~ t(covs_z):facs[,2] - 1)
          cov_and_seas <- rbind(seasonality[[sss]], covs_z, t(interactions))
        }
      } else {
        cov_and_seas <- rbind(seasonality[[sss]],covs_z)
      }

    #fit model with TMB
    if (is.null(seasonality[[sss]]) == TRUE & is.null(covariates[[cov]]) == TRUE){
      dfa <- runDFA(obs=dat_z, NumStates=mmm, ErrStruc=RRR,
                    EstCovar=FALSE)
    } else {
      dfa <- runDFA(obs=dat_z, NumStates=mmm, ErrStruc=RRR,
                    EstCovar=TRUE, Covars=cov_and_seas)
    }

    #save model object
    saveRDS(dfa, file=paste0("../model_objects_", tolower(y_choice), "/",
                             y_choice, '_', RRR, '_', mmm, 'm_',
                             names(seasonality)[sss], '_', names(covariates)[cov], '_',
                             startyr, '-', endyr, '.rds'))

    #landscape variable correlations
    if(!is.null(covs_z)){
      rescaled_effect_size <- eff_rescaler(cov_and_seas, seasonality[[sss]]) #effect sizes on original scale
      eff_best <- best_landvars(rescaled_effect_size, 6) #vars best correlated with effect size
    }
    loadings <- dfa$Estimates$Z
    load_best <- best_landvars(loadings, 6) #vars best correlated with loadings

    #format influential landscape var names and cors for export
    load_names1=load_names2=load_names3=load_names4=load_cors1=load_cors2=load_cors3=
      load_cors4=eff_names1=eff_names2=eff_cors1=eff_cors2=NA #placeholders

    for(i in 1:ncol(load_best[[1]])){
      assign(paste0('load_names', i), paste(load_best[[2]][,i], collapse=','))
      assign(paste0('load_cors', i), paste(load_best[[3]][,i], collapse=','))
    }
    if(!is.null(covs_z)){
      for(i in 1:ncol(eff_best[[1]])){
        assign(paste0('eff_names', i), paste(eff_best[[2]][,i], collapse=','))
        assign(paste0('eff_cors', i), paste(eff_best[[3]][,i], collapse=','))
      }
    }

    #open plot device
    pdf(file=paste0("../model_outputs_", tolower(y_choice), "/",
                    y_choice, '_', RRR, '_', mmm, 'm_', names(seasonality)[sss], '_',
                    names(covariates)[cov], '_',
                    startyr, '-', endyr, '.pdf'),
        onefile=TRUE)

    #plot hidden processes, loadings, model fits, and landscape regressions
    process_plotter_TMB(dfa, mmm)
    loading_plotter_TMB(dfa, mmm)
    fits_plotter_TMB(dfa)
    residuals_plotter(dfa)
    ACF_plotter(dfa)
    if(!is.null(covs_z)) eff_regress_plotter(mode='exploration')
    load_regress_plotter(mmm, mode='exploration')

    #close plot device
    dev.off()

    #get min, median, max R^2
    R2 <- get_R2(dfa)

    #store params, etc. in dataframe
    data.frame(R=RRR, m=mmm, cov=names(covariates)[cov],
               seasonality=names(seasonality)[sss],
               nparams=length(dfa$Optimization$par),
               ndatapts=length(unlist(which(!is.na(dat_z)))),
               n_p_ratio=length(unlist(which(!is.na(dat_z)))) / length(dfa$Optimization$par),
               NLL=dfa$Optimization$value, AIC=dfa$AIC,
               min_R2=R2[[1]], median_R2=R2[[2]], max_R2=R2[[3]],
               counts_func=unname(dfa$Optimization$counts[1]),
               counts_gradient=unname(dfa$Optimization$counts[2]),
               convergence=dfa$Optimization$convergence,
               message=paste('0',dfa$Optimization$message),
               load_names1=load_names1,load_names2=load_names2,
               load_names3=load_names3,load_names4=load_names4,
               load_cors1=load_cors1,load_cors2=load_cors2,
               load_cors3=load_cors3,load_cors4=load_cors4,
               eff_names1=eff_names1,eff_names2=eff_names2,
               eff_cors1=eff_cors1,eff_cors2=eff_cors2,
               stringsAsFactors=FALSE)
  }

#save data frame mentioned above
write.csv(model_out, file=paste0("../model_objects_", tolower(y_choice), "/",
                                 'param_tuning_dataframe_',
                                 startyr, '-', endyr, '_', y_choice, '.csv'))

stopCluster(cl) #free parallelized cores for other uses
