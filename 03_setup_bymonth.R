#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 8/10/2016

#NOTE - should have 4 datapoints (observations x streams) for each parameter in model.

rm(list=ls()); cat('\014')

# 0 - setup ####
setwd('C:/Users/Mike/git/stream_nuts_DFA/data/')
setwd('~/git/puget_sound_rivers_DFA/data')
setwd('Z:/stream_nuts_DFA/data/')
load('chemPhys_data/yys_bymonth.rda')
source('../00_tmb_uncor_Rmat.R')

#install packages that aren't already installed
package_list <- c('MARSS','viridis','imputeTS','vegan','cluster','fpc',
                  'RColorBrewer', 'foreach', 'doParallel')
new_packages <- package_list[!(package_list %in% installed.packages()[,"Package"])]
if(length(new_packages)) install.packages(new_packages)
if (!require("manipulateR")) {
    if (!require("devtools")) install.packages('devtools')
    library(devtools)
    install_github('vlahm/manipulateR') #this one is from github
}
for(i in c(package_list, 'manipulateR')) library(i, character.only=TRUE) #and load them all

# if (is.null(dev.list()) == TRUE){windows(record=TRUE)} #open new plot window unless already open

# 1 - CHOICES ####

# response choices: COND FC NH3_N NO2_NO3 OP_DIS OXYGEN PH PRESS SUSSOL TEMP TP_P TURB
y_choice = 'TEMP'
# cov choices: meantemp meantemp_anom precip precip_anom hydroDrought hydroDrought_anom
# maxtemp maxtemp_anom hdd hdd_anom
#this selects the full cov matrix during testing, but selects covs
#to be entered individually or in pairs during model fitting
# cov_choices = c('meantemp', 'precip', 'maxtemp', 'hydroDrought', 'hdd')
cov_choices = c('meantemp')
#region choices: '3' (lowland), '4' (upland), '3_4' (average of 3 and 4)
region = '3_4'
#average regions 3 and 4? (if FALSE, sites from each region will be assigned their own climate covariates)
average_regions = TRUE #only used if region = '3_4'
#seasonality model choices: 'fixed_collective', 'fixed_individual', 'fourier'
method = 'fixed_individual'
#which years to include?
startyr = 1978
endyr = 2015
#model params
ntrends = 2
obs_err_var_struc = 'diagonal and equal'

# 1.1 - subset datasets according to choices ####

#chem/phys data manipulations
# subset by year and exclude columns with >= na_thresh proportion of NAs
yy <- eval(parse(text=y_choice))
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
yy <- subsetter(yy, start=startyr, end=endyr, na_thresh=0.75) #beware reducing threshold, may add sites and break stuff

# subset by region
if(region == '3'){
    yy <- yy[,colnames(yy) %in% c('date','K','J','A','B','ZA','Q','H','T','G','F','E','C')]
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
                'hdd_anom'='hd_anom_1900.99')

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
    covs <- as.matrix(covs[substr(covs$date,1,4) >= startyr &
                               substr(covs$date,1,4) <= endyr,
                           colnames(covs) %in%
                               covdict[names(covdict) %in% cov_choices]])
}

# 2 - plot response variable time series (unintelligible by month) ####
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

# 3 - Set up input matrices to MARSS function call ####
# scale and center y and cov data
dat_z <- t(scale(as.matrix(obs_ts)))
mean(dat_z[1,], na.rm=T); sd(dat_z[1,], na.rm=T)

covs_z <- t(scale(as.matrix(covs)))

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

# 4 - run DFA (testing) ####

#MARSS full specification
dfa <- MARSS(y=dat_z, model=list(B=BB, U=uu, C='zero', c='zero', Q=QQ, Z=ZZ, A=aa, D='unconstrained', d=covs_z, R=RR),
             inits=list(x0=matrix(rep(0,mm),mm,1)),
             control=list(minit=200, maxit=20000, allow.degen=TRUE), silent=2)
dfa <- MARSS(y=dat_z, model=list(B=BB, U=uu, C=DD, c=dd, Q=QQ, Z=ZZ, A=aa, D=CC, d=cc, R=RR),
             inits=dfa$par,
             control=list(minit=200, maxit=3000), method='BFGS') #can't use BFGS for equalvarcov
#MARSS form=DFA
dfa2 <- MARSS(y=dat_z, model=list(m=2, R='diagonal and equal', A='zero'),#, D='unconstrained'),
             inits=list(x0='zero'), z.score=TRUE, #coef(dfa, type='matrix')$D
             control=list(minit=1, maxit=100, allow.degen=TRUE), silent=2, form='dfa')#,
             # covariates=cov_and_seas)
#TMB
dfa <- runDFA(obs=dat_z, NumStates=mm, ErrStruc='UNC', EstCovar=T, Covars=cov_and_seas)


# #get seasonal effects
# CC_out = coef(dfa, type="matrix")$C
# # The time series of net seasonal effects
# seas = CC_out %*% cc[,1:12]
# rownames(seas) = 1:mm
# colnames(seas) = month.abb
# seas

# 5 - plot estimated state processes, loadings, and model fits (MARSS-testing) ####

# varimax rotation to get Z loadings
Z_est <- coef(dfa, type="matrix")$Z # get the estimated ZZ
H_inv <- varimax(Z_est)$rotmat # get the inverse of the rotation matrix
Z_rot = Z_est %*% H_inv # rotate factor loadings
proc_rot = solve(H_inv) %*% dfa$states # rotate processes

# plot hidden processes
process_plotter <- function(){
    par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(mm, 1))
    xlbl <- int_dates
    y_ts <- int_dates
    ylm <- c(-1,1)*max(abs(proc_rot))
    for(i in 1:mm) {
        plot(y_ts,proc_rot[i,], type="n", bty="L",
             ylim=ylm, xlab="", ylab="", xaxt="n")
        abline(h=0, col="gray")
        lines(y_ts,proc_rot[i,], lwd=2)
        # lines(y_ts,proc_rot[i,] * seas[i,], lwd=2, col='green')
        mtext(paste("Process",i), side=3, line=0.5)
        axis(1, at=xlbl, labels=xlbl, cex.axis=0.8)
    }
}
process_plotter()

# plot loadings
loading_plotter <- function(){
    par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(mm, 1))
    ylbl <- names(obs_ts)
    clr <- viridis(nn) #colors may not line up with series plots in section 2
    ylm <- c(-1,1)*max(abs(proc_rot))
    minZ <- 0
    ylm <- c(-1,1)*max(abs(Z_rot))
    for(i in 1:mm) {
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
loading_plotter()

# get model fits & CI's
get_DFA_fits <- function(MLEobj,alpha=0.05) {
    ## empty list for results
    fits <- list()
    ## extra stuff for var() calcs
    Ey <- MARSS:::MARSShatyt(MLEobj)
    ## model params
    ZZ <- coef(MLEobj, type="matrix")$Z
    ## number of obs ts
    nn <- dim(Ey$ytT)[1]
    ## number of time steps
    TT <- dim(Ey$ytT)[2]
    ## get the inverse of the rotation matrix
    H_inv <- varimax(ZZ)$rotmat
    ## model expectation
    fits$ex <- ZZ %*% H_inv %*% MLEobj$states
    ## Var in model fits
    VtT <- MARSSkfss(MLEobj)$VtT
    VV <- NULL
    for(tt in 1:TT) {
        RZVZ <- coef(MLEobj, type="matrix")$R - ZZ%*%VtT[,,tt]%*%t(ZZ)
        SS <- Ey$yxtT[,,tt] - Ey$ytT[,tt,drop=FALSE] %*% t(MLEobj$states[,tt,drop=FALSE])
        VV <- cbind(VV,diag(RZVZ + SS%*%t(ZZ) + ZZ%*%t(SS)))
    }
    SE <- sqrt(VV)
    ## upper & lower (1-alpha)% CI
    fits$up <- qnorm(1-alpha/2)*SE + fits$ex
    fits$lo <- qnorm(alpha/2)*SE + fits$ex
    return(fits)
}
mod_fit <- get_DFA_fits(dfa)

# plot fits
fits_plotter <- function(){
    ylbl <- names(obs_ts)
    xlbl = y_ts = 1:444
    # par(mfrow=c(1,1), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
    par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
    ymin <- min(dat_z, na.rm=TRUE)
    ymax <- max(dat_z, na.rm=TRUE)
    for(i in 1:nn) {
        lo <- mod_fit$lo[i,]
        mn <- mod_fit$ex[i,]
        up <- mod_fit$up[i,]
        plot(y_ts,mn,xlab="",ylab=ylbl[i],xaxt="n",type="n", cex.lab=1.2,
             ylim=c(ymin,ymax))
        axis(1, at=xlbl, labels=xlbl, cex.axis=1)
        points(y_ts,dat_z[i,], pch=16, col="darkblue")
        lines(y_ts, up, col="darkgray")
        lines(y_ts, mn, col="black", lwd=2)
        lines(y_ts, lo, col="darkgray")
    }
}
fits_plotter()

# 5.1 - plot processes, loadings, fits, and get R^2 (TMB-testing) ####

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
process_plotter_TMB(dfa, mm)

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
loading_plotter_TMB(dfa, mm)

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
        points(dat_z[i,], col='blue', pch=20, cex=1.5)
    }
}
fits_plotter_TMB(dfa) #black is model fit, green is hidden-trend-only fit, blue is data

get_R2 <- function(){
    R2 <- rep(NA, nrow(dat_z))
    for(i in 1:nrow(dat_z)){
        mod <- lm(dat_z[i,] ~ dfa$Fits[i,])
        R2[i] <- summary(mod)$r.squared
    }
    return(list(min(R2), median(R2), max(R2)))
}
get_R2()

# 6 - model fitting/parameter tuning (MARSS - not parallelized) ####

R_strucs <- c('diagonal and equal', 'diagonal and unequal', 'equalvarcov')
R_names <- c('DiagAndEq', 'DiagAndUneq', 'EqVarCov')
ntrends <- 1:4
model_out <- data.frame()

# ntrends <- 5
# model_out <- data.frame()
# R_strucs <- c('diagonal and unequal')
# R_strucs <- c('unconstrained')
# RRR='diagonal and unequal'
# RRR='unconstrained'
# R_names <- c('DiagAndUneq')
# R_names <- c('Unconst')
# mmm=1

for(RRR in R_strucs){
    for(mmm in ntrends){

        #create loadings matrix
        mm <- mmm
        ZZZ <- ZZgen()

        #fit model with EM algorithm
        print(paste(RRR,mmm))


        if(mmm == 1){
            print(paste('should be 1', mmm))

            dfa <- try(MARSS(y=dat_z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
                             inits=list(x0=0), silent=FALSE,
                             control=list(maxit=20000, allow.degen=TRUE)))
            if(isTRUE(class(dfa)=='try-error')) {next}
        } else {
            print(mmm)

            dfa <- try(MARSS(y=dat_z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
                             inits=list(x0=matrix(rep(0,mmm),mmm,1)), silent=FALSE,
                             control=list(maxit=20000, allow.degen=TRUE)))
            if(isTRUE(class(dfa)=='try-error')) {next}
        }

        #where possible, polish EM estimate with BFGS
        if(RRR != 'equalvarcov'){
            print(paste(RRR,mmm,'BFGS'))

            dfa <- try(MARSS(y=dat_z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
                             inits=dfa$par, silent=FALSE,
                             # inits=coef(dfa, form='marss'), #alternate form? - see MARSSoptim examples
                             control=list(dfa$control$maxit), method='BFGS'))
            if(isTRUE(class(dfa)=='try-error')) {next}
        }

        #store params, etc. in dataframe
        model_out <- rbind(model_out, data.frame(R=RRR, m=mmm, LogLik=dfa$logLik,
                                                 K=dfa$num.params, AICc=dfa$AICc,
                                                 AIC=dfa$AIC, converged=dfa$convergence,
                                                 filter=dfa$fun.kf, method=dfa$method,
                                                 iter=unname(dfa$numIter), n=dfa$samp.size,
                                                 stringsAsFactors=FALSE))

        #save model object
        saveRDS(dfa, file=paste0("../stream_nuts_DFA/model_objects/",
                                 R_names[R_strucs == RRR], '_', mmm, 'm_', startyr, y_choice, '.rds'))

        if(mmm > 1){

            #open plot device
            pdf(file=paste0("../stream_nuts_DFA/model_outputs/",
                            R_names[R_strucs == RRR], '_', mmm, 'm_', startyr, y_choice, '.pdf'), onefile=TRUE)

            # varimax rotation to get Z loadings
            Z_est <- coef(dfa, type="matrix")$Z # get the estimated ZZ
            H_inv <- varimax(Z_est)$rotmat # get the inverse of the rotation matrix
            Z_rot = Z_est %*% H_inv # rotate factor loadings
            proc_rot = solve(H_inv) %*% dfa$states # rotate processes

            #plot hidden processes, loadings, and model fits for each parameter combination
            process_plotter()
            loading_plotter()

            mod_fit <- get_DFA_fits(dfa)
            fits_plotter()

            #close plot device
            dev.off()

        }
    }
}
#save data frame
write.csv(model_out, file=paste0("../stream_nuts_DFA/model_objects/",
                                 'param_tuning_dataframe_', startyr, y_choice, '.csv'))

# 6.1 - model fitting/parameter tuning (TMB) ####
#round 2 determined that individual fixed factors provide the best seasonality absorption

#specify system-dependent cluster type; ncores-1 to be used in parallel
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
# unregister()

R_strucs <- c('DE','DUE','UNC')
ntrends <- 1:4
seasonality <- list('fixed_factors'=ccgen('fixed_individual'),
                    'fourier'=ccgen('fourier'), 'no_seas'=NULL)

sss <- 1

model_out <-
    foreach(RRR=R_strucs, .combine=rbind) %:%
        foreach(mmm=ntrends, .combine=rbind) %:%
            # foreach(sss=1:length(seasonality), .combine=rbind,
            foreach(
                    .packages='viridis') %dopar% {

                source('../00_tmb_uncor_Rmat.R')

                print(paste(RRR,mmm,names(seasonality)[sss]))

                #fit model with TMB
                cov_and_seas <- rbind(seasonality[[sss]],covs_z)
                dfa <- runDFA(obs=dat_z, NumStates=mmm, ErrStruc=RRR,
                              EstCovar=TRUE, Covars=cov_and_seas)

                #save model object
                saveRDS(dfa, file=paste0("../model_objects/",
                                         RRR, '_', mmm, 'm_',
                                         names(seasonality)[sss], '_',
                                         startyr, '-', endyr, '_', y_choice, '.rds'))

                #open plot device
                pdf(file=paste0("../model_outputs/",
                                RRR, '_', mmm, 'm_', names(seasonality)[sss], '_',
                                startyr, '-', endyr, '_', y_choice, '.pdf'),
                    onefile=TRUE)

                #plot hidden processes, loadings, and model fits for each parameter combination
                process_plotter_TMB(dfa, mmm)
                loading_plotter_TMB(dfa, mmm)
                fits_plotter_TMB(dfa)

                #close plot device
                dev.off()

                #get min, median, max R^2
                R2 <- get_R2

                #store params, etc. in dataframe
                data.frame(R=RRR, m=mmm,
                           seasonality=names(seasonality)[sss],
                           nparams=length(dfa$Optimization$par),
                           LogLik=dfa$Optimization$value, AIC=dfa$AIC,
                           min_R2=R2[1], median_R2=R2[2], max_R2=R2[3],
                           counts_func=unname(dfa$Optimization$counts[1]),
                           counts_gradient=unname(dfa$Optimization$counts[2]),
                           convergence=dfa$Optimization$convergence,
                           message=paste('0',dfa$Optimization$message),
                           stringsAsFactors=FALSE)
                }

#save data frame
write.csv(model_out, file=paste0("../model_objects/",
                                 'param_tuning_dataframe_', startyr, '-',
                                 endyr, '_', y_choice, '.csv'))

stopCluster(cl) #free parallelized cores for other uses

# 6.2 - load desired model output ####
# mod_out <- readRDS("../round_2_tmb_covs/model_objects/UNC_2m_fixed_factors_1978-2014_TEMP.rds")
dfa <- readRDS("../round_2_tmb_covs/model_objects/UNC_2m_fixed_factors_1978-2014_TEMP.rds")

# 7 - landscape factor regressions (setup) ####

#regress factor loadings (Z) on hidden trends against landscape vars
land <- read.csv('watershed_data/watershed_data_simp.csv', stringsAsFactors=FALSE)
land$siteCode[land$siteCode == 'AA'] <- 'ZA' #rename sammamish @ bothell sitecode
land <- land[land$siteCode %in% names(obs_ts),] #remove sites not in analysis
land <- land[match(names(obs_ts), land$siteCode),] #sort landscape data by site order in model

#plot trend loadings against all landscape variables of interest
landvars <- c('BFIWs','ElevWs','PctImp2006WsRp100',
              'PctGlacLakeFineWs','PctAlluvCoastWs','PctIce2011Ws',
              'PctCrop2011Ws', 'PctUrbOp2011WsRp100','PctUrbLo2011WsRp100',
              'PctUrbMd2011WsRp100','PctUrbHi2011WsRp100',
              'RdDensWsRp100','RunoffWs','OmWs',
              'RckDepWs','WtDepWs','PermWs','PopDen2010Ws')
landcols <- which(colnames(land) %in% landvars)

#multiple covariates?
# for(i in 1:ncol(obs_ts)){
#     sd_response <- sd(obs_ts[,i], na.rm=TRUE)
#     sd_covar <- sd(covs[,i], na.rm=TRUE)
#     if(exists('cc')){
#         dfa$Estimates$D[,]
#     }
#     effect_size <-  ( / )
# }

#effect size regression (skip if no covars) ####
nstream <- ncol(obs_ts)
ncov <- ncol(covs)

#get covariate effect sizes (D coefficients) from model, isolated from seasonal effects
if(nrow(cov_and_seas) > 2){
    z_effect_size <- matrix(dfa$Estimates$D[,(nrow(cc)+1):ncol(dfa$Estimates$D)])
} else {
    z_effect_size <- dfa$Estimates$D
}

#convert effects back to original scale (units response/units covar)
#x/sd(response) = 1/sd(covar); x is the coefficient scale factor
for(i in 1:nstream){
    sd_response <- sd(obs_ts[,i], na.rm=TRUE)
    rescaled_effect_size <- matrix(NA, nrow=nstream, ncol=ncov)
    for(j in 1:ncov){
        sd_covar <- sd(covs[,j], na.rm=TRUE)
        rescaled_effect_size[,j] <- z_effect_size[,j] * (sd_response/sd_covar)
    }
}

#grab top 8 landscape vars that correlate with effect size, plot and get stats ####
eff_cors <- unname(apply(land[landcols], 2, function(i) cor(rescaled_effect_size, i)))
eff_rank <- order(abs(eff_cors), decreasing=TRUE)
best_landvars <- landcols[head(eff_rank, top)]

top <- 8 #top x columns to plot/test
for(i in 1:){
    plot(mod_out$Estimates$D[??], land[,best_landvars])
}
#left off here. check dim(D), plot land (x) vs rescaled effect size (y) - dunno
#how this will work for multiple covariates...


for(i in landcols){
    lm(mod_out$Estimates$Z[,1], land[,i])
}
for(i in landcols){
    plot(mod_out$Estimates$Z[,1], land[,i])
}

plot(covs_z, obs_ts$A)
summary(lm(obs_ts$A ~ t(covs_z)))

#regress effect sizes of climate covariates against landscape vars


#Z loadings
for(i in 1:length(landcols)){
    plot(mod_out$Estimates$Z[,1], land[,landcols[]])
}
