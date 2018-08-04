#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 8/10/2016
#updated: 8/3/2018

#NOTES:

#collapse folds with ALT+O (windows, linux) or CMD+OPT+O (Mac); might have to
#do it twice

#you should have at least 4 datapoints (observations x streams) for each
#parameter in model.
#you can see the datapoints:params ratio in the "results" dataframe in section 6.2

#if R crashes when you try to use runDFA,
#use apply(dat_z, 2, function(x) sum(is.na(x))/length(x)) to see if you have
#any timepoints with
#no data or 1 data point. these timepoints must either be removed or
#imputed (see section 1.2).

rm(list=ls()); cat('\014') #clear env and console

# 0 - setup ####
setwd('C:/Users/vlahm/Desktop/git/puget_sound_rivers_DFA/data/')

load('chemPhys_data/yys_bymonth.rda')
DISCHARGE <- read.csv('discharge_data/discharge.csv', stringsAsFactors=FALSE, colClasses=c('date'='Date'))
snowmelt <- read.csv('climate_data/snow_data/snowmelt.csv')

#install packages that aren't already installed (see https://github.com/kaskr/adcomp for TMB package)
#imputeTS, RColorBrewer, cluster, fpc
package_list <- c('MARSS','viridis','vegan', 'e1071', 'imputeTS', 'stringr',
    'foreach', 'doParallel', 'caret', 'Matrix')
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages, repos="http://cran.rstudio.com/")

# 1 - OPTIONS ####

# response choices: COND FC NH3_N NO2_NO3 OP_DIS OXYGEN PH PRESS SUSSOL
#TEMP TP_P TURB
# also DISCHARGE (from USGS)
y_choice = 'NO2_NO3'
# cov choices: meantemp meantemp_anom precip precip_anom hydroDrought
#hydroDrought_anom
# maxtemp maxtemp_anom hdd hdd_anom, snowmelt
#(snowmelt only available 1978-2015. also I haven't actually used snowmelt
#in a model yet, so there could be bugs with it)
#specify all of the covs here that will be used later in the fitting loop
cov_choices = c('meantemp', 'precip', 'snowmelt')
#region choices: '3' (lowland), '4' (upland), '3_4'
#(average of 3 and 4, or each separately)
#regions 3 and 4 were discovered to be very similar early on,
#so most of this script will only work if
#you choose '3_4' here and 'average_regions=TRUE'.
region = '3_4' #code not set up to include snowmelt unless region='3_4' and
#average_regions=TRUE
#average regions 3 and 4? (if FALSE, sites from each region will be assigned
#their own climate covariates)
average_regions = TRUE #only used if region = '3_4'
#which years to include?
startyr = 1978
endyr = 2015
#UPDATE: Mark Schueurell no longer scales his response data. scaling forces the variance of
#the D matrix to be small, thus artificially diminishing the impact of the covariates.
#Scaling is still recommended by the community at large, and I didn't see a
#big difference either way.
scale = FALSE
#exclude sites with >= this proportion of NA values.
na_thresh = 0.55
#transformations are 'log' and 'none'. you can also explore 'power' and 'boxcox' in section 3.1
#run function transformables() in section 1.3 to see whether your response needs to be transformed.
transform = 'none'

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

# subset by region if requested
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

#i just added this part for interpolating NO2_NO3 on 8/3/18.
#There's a huge chunk of missing data for all sites in the 80s (god damned reagan)
if(y_choice == 'NO2_NO3'){
    yts <- ts(obs_ts, start=1, frequency=12)
    obs_ts_pre <- na.seasplit(yts, 'interpolation')
    obs_ts <- data.frame(obs_ts_pre)
}

#make sure it worked as expected (the purple parts are imputed)
defpar <- par(mfrow=c(3,1))
for(i in 1:3){
    plot(obs_ts_pre[,i], col='purple')
    lines(yts[,i], col='orange')
}
par(defpar)

# 1.3 - (transform), center, (scale) ####
library(caret); library(e1071)

#visualize which responses should be transformed to normal (looks like
#NO2_NO3 could use a transform)
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
transformables()

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
trans <- transformer(obs_ts, transform=transform, exp=NA, scale=scale, plot=F)

dat_z <- t(trans$trans)
mean(dat_z[1,], na.rm=T); sd(dat_z[1,], na.rm=T) #verify centering and scaling (if chosen)

#scale and center covariate data
covs_z <- t(scale(as.matrix(covs)))

# 3.0 - some leftover stuff from MARSS experimentation ####
#just some functions for generating c and Z matrices
ccgen <- function(meth=method){
    if(meth %in% c('fixed_collective', 'fixed_individual')){
        year_block <- diag(12)
        nyears <- endyr-startyr+1

        cc <- Reduce(function(x,y) {cbind(x,y)},
            eval(parse(text=paste0('list(',
                paste(rep('year_block', nyears), sep=', ', collapse=', '),
                ')'))))
        rownames(cc) <- month.abb
    } else { #fourier method
        cos_t = cos(2 * pi * seq(dim(dat_z)[2]) / 12)
        sin_t = sin(2 * pi * seq(dim(dat_z)[2]) / 12)
        cc = rbind(cos_t,sin_t)
    }
    return(cc)
} #(needs work - set up for season and add climate data?)
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
nn <- length(names(obs_ts))
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
# cov_and_seas <- designer(design=design, sections=sections, focal_months=focal_months)

# 5.1 - plot processes, loadings, fits, residuals, ACF, and get R^2 (TMB) ####
library(viridis)

#these functions are used later within the loop

process_plotter_TMB <- function(dfa_obj, ntrends){
    if(ntrends<=4){
        par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(ntrends, 1))
    } else par(mai=c(0,0,0,0), omi=c(0,0,0,0), mfrow=c(ntrends, 1))
    xlbl <- int_dates
    y_ts <- int_dates
    ylm <- c(-1,1)*max(abs(dfa_obj$Estimates$u))
    for(i in 1:ntrends){
        plot(y_ts,dfa_obj$Estimates$u[i,], type="n", bty="L",
            ylim=ylm, xlab="", ylab="", xaxt="n")
        abline(h=0, col="gray")
        lines(y_ts,dfa_obj$Estimates$u[i,], lwd=2)
        mtext(paste("Process",i), side=3, line=-2)
        axis(1, at=xlbl, labels=xlbl, cex.axis=0.8)
    }
}
# process_plotter_TMB(dfa, mm)

loading_plotter_TMB <- function(dfa_obj, ntrends){
    if(ntrends<=4){
        par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(ntrends, 1))
    } else par(mai=c(0,0,0,0), omi=c(0,0,0,0), mfrow=c(ntrends, 1))
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
        mtext(paste("Factor loadings on process",i),side=3,line=-2)
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
#I tried to make this function account for all the different covariate designs (except ..._acrossTime), but
#there are tons of different contingencies to take into account, so I may have missed something.
#if your effect size regressions don't look like what you expected, start troubleshooting here
#(again, I'm happy to help). If you use design='effect_byMonth', this will
#assume you have seasonality method='fixed...'. it might also work with fourier. but you'll need
#seasonality to be included. If you used 'effect_byMonth_noSeas', you're fine.
eff_rescaler <- function(all_cov, seas, ncov, scaled=scale){
    #get covariate effect sizes (D coefficients) from model, isolated from seasonal effects
    if(design == 'effect_byMonth'){
        z_eff_1 <- NULL
        if(ncov > 1){
            z_eff_1 <- as.matrix(dfa$Estimates$D[,(nrow(seas)+1):(nrow(seas)+ncov-1)])
        }
        z_eff_2 <- as.matrix(rowMeans(dfa$Estimates$D[,(nrow(seas)+ncov):
                (nrow(seas)+ncov-1+length(focal_months))]))
        z_effect_size <- cbind(z_eff_2, z_eff_1)
    } else {
        if(design == 'effect_byMonth_noSeas'){
            z_eff_1 <- NULL
            if(ncov > 1){
                z_eff_1 <- as.matrix(dfa$Estimates$D[,1:(ncov-1)])
            }
            z_eff_2 <- as.matrix(rowMeans(dfa$Estimates$D[,ncov:
                    (ncov-1+length(focal_months))]))
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

# 6.1 - model fitting/parameter tuning (TMB; parallelized) ####
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

#create vectors of parameter values to loop through:

#choose observation error structures to include.
#ISSUE: during parallelization, each new thread runs its own separate TMB script, but this original thread (the
#one from which you start the parallelizing loop) also has to run a TMB script, and I don't know how to make it run
#different ones depending on which R_strucs is assigned to it. So, EVCV should be run in its own separate
#parallel loop.
R_strucs <- c('DUE')#, 'DE') #'UNC'
# R_strucs <- c('EVCV')

#number of common trends to fit
ntrends <- 3:4

#should seasonality be modeled as fixed factors for each month ('fixed_individual'), as
#a Fourier process ('fourier'), or not at all (NULL)?
#CAUTION: if you choose a design of 'effect_byMonth_acrossTime' or 'effect_byMonth'below, you can
#only use 'fixed_individual' here.
seasonality <- list(fixed_factors=ccgen('fixed_individual'))
#fourier=ccgen('fourier'))#, no_seas=NULL)

#subset particular covariates from the full covariate matrix, depending on what you included
#in "OPTIONS" (section 1) above. Can be NULL.
#if you're looking at effects by month, or by month over time (see 'design' option below),
#then whichever covariate is listed first for each list element will be the one that can be examined by month/
#over time. any others will still be included, and their average effect sizes (not by month, not over time) can
#still be determined.
covs_z2 <- covs_z #for name preservation elsewhere
covariates <- list(pc=covs_z2[2,, drop=FALSE],#, sn=rbind(covs_z2[3,]),
    at=covs_z2[1,, drop=FALSE],
    # pcsn=covs_z2[2:3,])
    atpc=covs_z2[1:2,])#, atpcsn=covs_z2[1:3,],
# atsn=covs_z2[c(1,3),])

#the 3 options below will not vary among the models in the loop. only change these between runs.

#choose the covariate matrix design here. options are 'just_effect', 'effect_and_seasonality_without_interaction'
#'effect_byMonth', 'effect_byMonth_noSeas', and 'effect_byMonth_acrossTime'.
#(see "designer" function in section 3.1 for details)
design = 'effect_byMonth'
#sections = number of intervals to divide the time series into, if examining change over time
#(see "designer" function in section 3.1 for details).
#value must be an integer. will be ignored unless design='effect_byMonth_acrossTime'.
sections <- 5
#focal_months = the months to focus on for by-month effect size (1 is jan...).
#if looking at effect_byMonth_acrossTime, including all months will be too expensive
#(see "designer" function in section 3.1 for details).
#value should be a vector of integers between 1 and 12. will be ignored unless
#design='effect_byMonth' or 'effect_byMonth_acrossTime'. If the former, you can probably
#afford to include all months. Otherwise it might get too costly. You can track this by
#examining the parameter to data ratio in the output table.
focal_months <- 1:12

# troubleshooting:
# sss = 1 #uncomment to fix seasonality (must also comment **s below)
# RRR = 'UNC' #uncomment to fix error structure (must also comment **R below)
# RRR='DE'; mmm=1; cov=3; sss=1 #convenience variables for stepping through loop
# rm(RRR, mmm, cov, sss)
# rm(cov_and_seas, dfa, all_cov, seas)

#start timing
ptm <- proc.time()

#source the appropriate TMB script based on the chosen error structure
if('EVCV' %in% R_strucs) source('../00_tmb_uncor_Rmat_EVCV.R') else source('../00_tmb_uncor_Rmat_DE_DUE_UNC.R')

#this is the model fitting loop. it opens new R processes for each processor core
#and allocates model runs efficiently to each. the outputs are stored in
#stream_nuts_DFA/model_objects_<lower(y_choice)> and
#stream_nuts_DFA/model_outputs_<lower(y_choice)>. for example, if you're looking at DISCHARGE,
#you'll need to create directories called model_objects_discharge and model_outputs_discharge.
#The variable "model_out" is a summary data frame of diagnostic stuff, which will be stored with
#model objects.
model_out <-
    foreach(RRR=R_strucs, .combine=rbind) %:% # **R
    foreach(mmm=ntrends, .combine=rbind) %:%
    foreach(sss=1:length(seasonality), .combine=rbind) %:% # **s
    foreach(cov=1:length(covariates), .combine=rbind,
        .packages='viridis') %dopar% {

            #re-run TMB script for each new R process
            if('EVCV' %in% R_strucs) source('../00_tmb_uncor_Rmat_EVCV.R') else
                source('../00_tmb_uncor_Rmat_DE_DUE_UNC.R')

            #create covariate matrix (d) - see function definition in section 3.1 for details
            covs_z <- covariates[[cov]]
            cov_and_seas <- rbind(seasonality[[sss]], covs_z)
            #keeping name "cov_and_seas" here for compatibility elsewhere
            cov_and_seas <- designer(design=design, sections=sections, focal_months=focal_months)

            #fit model with TMB
            if (is.null(seasonality[[sss]]) == TRUE & is.null(covs_z) == TRUE){
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

runtime <- proc.time() - ptm

if(exists('model_out')){
    writeLines(paste0('Loop ran to completion.\nElapsed time: ',
        round(unname(runtime[3])/60, 2),
        ' minutes.\nRestart Rstudio before using the loop again.'))
} else writeLines(paste("Check out the 'troubleshooting' chunk just above the loop.\n",
    "You'll probably have to restart Rstudio before trying again."))

# 6.2 - model fitting/parameter tuning (TMB; one-off; created 8/3/18) ####

#error structures 'DE' (diagonal and equal), 'DUE'
#(diag and uneq), 'UNC' (unconstrained), 'EVCV' (equal var covar)
R_struct = 'DUE'

#number of common trends to fit
ntrends <- 2

#should seasonality be modeled as fixed factors for each month ('fixed_individual'), as
#a Fourier process ('fourier'), or not at all (NULL)?
#CAUTION: if you choose a design of 'effect_byMonth_acrossTime' or 'effect_byMonth'below, you can
#only use 'fixed_individual' here.
seasonality <- ccgen('fixed_individual')

#subset particular covariates from the full covariate matrix, depending on what you included
#in "OPTIONS" (section 1) above. Can be NULL.
#if you're looking at effects by month, or by month over time (see 'design' option below),
#then whichever covariate is listed first for each list element will be the one that can be examined by month/
#over time. any others will still be included, and their average effect sizes (not by month, not over time) can
#still be determined.
covariates = covs_z[1:2,] #air temperature and precipitation
# covariates = covs_z[1,, drop=FALSE] #air temperature
# covariates = covs_z[2,, drop=FALSE] #precipitation
# covariates = covs_z[3,, drop=FALSE] #snowmelt

#choose the covariate matrix design here. options are 'just_effect', 'effect_and_seasonality_without_interaction'
#'effect_byMonth', 'effect_byMonth_noSeas', and 'effect_byMonth_acrossTime'.
#(see "designer" function in section 3.1 for details)
design = 'effect_byMonth'

#sections = number of intervals to divide the time series into, if examining change over time
#(see "designer" function in section 3.1 for details).
#value must be an integer. will be ignored unless design='effect_byMonth_acrossTime'.
sections <- 5

#focal_months = the months to focus on for by-month effect size (1 is jan...).
#if looking at effect_byMonth_acrossTime, including all months will be too expensive
#(see "designer" function in section 3.1 for details).
#value should be a vector of integers between 1 and 12. will be ignored unless
#design='effect_byMonth' or 'effect_byMonth_acrossTime'. If the former, you can probably
#afford to include all months. Otherwise it might get too costly. You can track this by
#examining the parameter to data ratio in the output table.
focal_months <- 1:12

#source the appropriate TMB script based on the chosen error structure
if('EVCV' %in% R_struct){
    source('../00_tmb_uncor_Rmat_EVCV.R')
} else {
    source('../00_tmb_uncor_Rmat_DE_DUE_UNC.R')
}

#create covariate matrix (d) - see function definition in section 3.1 for details
cov_and_seas <- rbind(seasonality, covariates)
#designer modifies the existing cov_and_seas matrix
#(bad practice not to pass it as argument, i know)
cov_and_seas <- designer(design=design, sections=sections,
    focal_months=focal_months)

#fit model with TMB
dfa <- runDFA(obs=dat_z, NumStates=ntrends, ErrStruc=R_struct,
    EstCovar=TRUE, Covars=cov_and_seas)

#landscape variable correlations
rescaled_effect_size <- eff_rescaler(cov_and_seas, seasonality,
    ncov=nrow(covariates)) #effect sizes on original scale
eff_best <- best_landvars(rescaled_effect_size, 6) #vars best correlated with effect size

loadings <- dfa$Estimates$Z
load_best <- best_landvars(loadings, 6) #vars best correlated with loadings

#plot hidden processes, loadings, model fits, and landscape regressions
process_plotter_TMB(dfa, ntrends)
loading_plotter_TMB(dfa, ntrends)
fits_plotter_TMB(dfa)
residuals_plotter(dfa)
ACF_plotter(dfa)
eff_regress_plotter(mode='exploration')
load_regress_plotter(ntrends, mode='exploration')

#get min, median, max R^2
R2 <- get_R2(dfa)

#store params, etc. in dataframe
results = data.frame(R=R_struct, m=ntrends,
    nparams=length(dfa$Optimization$par),
    ndatapts=length(unlist(which(!is.na(dat_z)))),
    n_p_ratio=length(unlist(which(!is.na(dat_z)))) / length(dfa$Optimization$par),
    NLL=dfa$Optimization$value, AIC=dfa$AIC,
    min_R2=R2[[1]], median_R2=R2[[2]], max_R2=R2[[3]],
    counts_func=unname(dfa$Optimization$counts[1]),
    counts_gradient=unname(dfa$Optimization$counts[2]),
    convergence=dfa$Optimization$convergence,
    message=paste('0',dfa$Optimization$message),
    stringsAsFactors=FALSE)
