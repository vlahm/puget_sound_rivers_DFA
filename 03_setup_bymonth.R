#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 8/10/2016
#last edit: 8/22/2016

###
#change the names of DDgen and ddgen if it turns out that seasonal effects go in the process eqn
###

#haven't set up fourier CC matrix yet. adding monthly effects in cc but not CC doesn't do what i need.
rm(list=ls()); cat('\014')

# 0 - setup ####
setwd('C:/Users/Mike/git/stream_nuts_DFA/data/')
setwd('Z:/stream_nuts_DFA/data/')
load('chemPhys_data/yys_bymonth_mean.rda')

library(MARSS)
library(viridis)
library(imputeTS)
library(vegan)
library(cluster)
library(fpc)
library(RColorBrewer)
if (!require("manipulateR")) {
    if (!require("devtools")) install.packages('devtools')
    library(devtools)
    install_github('vlahm/manipulateR')
}
library(manipulateR)

# if (is.null(dev.list()) == TRUE){windows(record=TRUE)} #open new plot window unless already open

# 1 - CHOICES ####

# response choices: COND FC NH3_N NO2_NO3 OP_DIS OXYGEN PH PRESS SUSSOL TEMP TP_P TURB
y_choice = 'TEMP'
# cov choices: meantemp meantemp_anom precip precip_anom hydroDrought hydroDrought_anom
# maxtemp maxtemp_anom hdd hdd_anom
cov_choices = c('meantemp', 'precip')
#region choices: '3' (downstream), '4' (upstream), '3_4' (average of 3 and 4)
region = '3_4'
#average regions 3 and 4? (if FALSE, sites from each region will be assigned their own climate covariates)
average_regions = TRUE #only used if region = '3_4'
#seasonality model choices: 'fixed_collective', 'fixed_individual', 'fourier'
method = 'fixed_individual'
#which years to include?
startyr = 1978
endyr = 2014
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
if(region == 3){
    yy <- yy[,colnames(yy) %in% c('date','K','J','A','B','ZA','Q','H','T','G','F','E','C')]
} else {
    if(region == 4){
        yy <- yy[,colnames(yy) %in% c('date','Z','I','L','M','N','O','P','S','U','X','V','W')]
    }
}

#convert dates to integers
months <- as.numeric(format(yy[,1], '%m'))
years <- as.numeric(format(yy[,1], '%Y'))
int_dates <- (years - startyr) * 12 + months

#create objects for response variables and covariates
obs.ts <- yy[,-1]

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

# 1.2 - identify seasons (cluster time points into 4 groups by water temp) ####

#check for outliers and excessive noise
fivenum(obs.ts, na.rm=TRUE)
sum(obs.ts > 20, na.rm=TRUE)
sd(as.matrix(obs.ts), na.rm=TRUE)
#no outliers detected; reasonable noise; proceeding with k-means clustering

#remove columns with lots of NAs before imputing values
matfilter(obs.ts, fun=is.na, na.rm=FALSE, filter=FALSE)
obs_sub <- obs.ts[,-c(2,8,16,18,21,24,25)]

#impute values
for(i in 1:ncol(obs_sub)){
    y <- obs_sub[,i]
    yts <- ts(y, int_dates[1], int_dates[length(y)], frequency=12)
    yimp <- na.seasplit(yts, 'interpolation')
    # plot(int_dates, yimp[1:length(y)], col='green', type='l', lwd=2)
    # lines(int_dates, y, lwd=2)
    obs_sub[,i] <- round(yimp[1:length(y)], 2)
}

#k-means clustering (based on euclidean distance)
euc <- vegdist(obs_sub, method='euclidean')
kmns3 <- kmeans(obs_sub, centers=3, iter.max=10000, nstart=25)

#plot cluster output (jul-sep at one end, dec-feb at the other, mar-jun continuous, oct-nov continuous)
cols <- brewer.pal(12, 'Paired')
pchvec <- as.numeric(as.vector(factor(kmns3$cluster, labels=c(15,16,17))))
par(bg='black')
plotcluster(obs_sub, kmns3$cluster, col=cols, pch=pchvec)
legend('topleft', legend=month.abb, fill=cols, bg='white')
# kmns4 <- kmeans(obs_sub, centers=4, iter.max=10000, nstart=25)
# plotcluster(obs_sub, kmns4$cluster)

#aggregate data by season
seas_bymo <- rep(c(4,4,1,1,1,1,2,2,2,3,3,4), length.out=444) #spr-win = 1-4

obs_seas <- aggregate(obs.ts, list(seas_bymo, years), mean, na.rm=TRUE)
obs_seas[is.na(obs_seas)] <- NA #replace NaNs with NAs
colnames(obs_seas)[1] <- 'season'

cov_seas <- aggregate(covs, list(seas_bymo, years), mean)
colnames(cov_seas)[1:2] <- c('season', 'year')

seasons <- cov_seas$season
years_seas <- cov_seas$year

# 2 - plot response variable time series (unintelligible by month)####
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
dat.z <- t(scale(as.matrix(obs.ts)))

Sigma = sqrt(apply(as.matrix(obs.ts), 1, var, na.rm=TRUE))
y.bar = apply(as.matrix(obs.ts), 1, mean, na.rm=TRUE)
dat.z2 = (as.matrix(obs.ts) - y.bar) * (1/Sigma)
identical(t(dat.z), dat.z2); mean(dat.z[1,], na.rm=T); mean(dat.z2[1,], na.rm=T)

covs <- t(scale(as.matrix(covs)))

# state equation params
mm <- ntrends #number of hidden processes (trends)
BB <- "identity" #'BB' is identity: 1's along the diagonal & 0's elsewhere (not part of DFA)
uu <- "zero" # 'uu' is a column vector of 0's (not part of DFA)
CCgen <- function(){
    if(method %in% c('fixed_collective', 'fourier_collective')){
        out <- matrix(month.abb, mm, 12, byrow=TRUE)
        rownames(out) <- NULL
    } else {
        if(method %in% c('fixed_individual', 'fourier_individual')){
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
ccgen <- function(){
    if(method %in% c('fixed_collective', 'fixed_individual')){
        year_block <- diag(12)
        nyears <- endyr-startyr+1

        cc <- Reduce(function(x,y) {cbind(x,y)},
                     eval(parse(text=paste0('list(',
                                            paste(rep('year_block', nyears,), sep=', ', collapse=', '),
                                            ')'))))
        rownames(cc) <- month.abb
    } else { #fourier method
        cos_t = cos(2 * pi * seq(dim(dat.z)[2]) / 12)
        sin_t = sin(2 * pi * seq(dim(dat.z)[2]) / 12)
        cc = rbind(cos_t,sin_t)
    }
    return(cc)
}
cc <- ccgen() #covariates for the state processes (seasonal and climatic)
# cc <- "zero" #use zero if not including seasonality
QQ <- "identity"  # 'QQ' is identity (would usually be tuned if we were doing actual marss)

# observation equation params
nn <- length(names(obs.ts)) # number of obs time series
aa <- "zero" # 'aa' is the offset/scaling (zero allowed because we standardize the data)
#(we want zero because we're not interested in what the offsets actually are)
DD <- 'unconstrained'
# DD <- 'zero' #use zero if not including climate covs
dd <- covs
# dd <- 'zero' #use zero if not including climate covs
RR <- obs_err_var_struc # 'RR' is var-cov matrix for obs errors (tune this)
ZZgen <- function(){
    ZZ <- matrix(list(0),nn,mm)
    ZZ[,1] <- paste0("z",names(obs.ts),1)
    if(mm > 1){
        for(i in 2:mm){
            ZZ[i:nn,i] <- paste0("z",names(obs.ts)[-(1:(i-1))],i)
        }
    }
    return(ZZ)
}
ZZ <- ZZgen() # 'ZZ' is loadings matrix (some elements set to zero for identifiability)

# 4 - run DFA ####
dfa <- MARSS(y=dat.z, model=list(B=BB, U=uu, C='zero', c='zero', Q=QQ, Z=ZZ, A=aa, D='unconstrained', d=covs, R=RR),
             inits=list(x0=matrix(rep(0,mm),mm,1)),
             control=list(minit=200, maxit=20000, allow.degen=TRUE), silent=2)
dfa <- MARSS(y=dat.z, model=list(B=BB, U=uu, C=DD, c=dd, Q=QQ, Z=ZZ, A=aa, D=CC, d=cc, R=RR),
             inits=dfa$par,
             control=list(minit=200, maxit=3000), method='BFGS') #can't use BFGS for equalvarcov

dfa <- MARSS(y=dat.z, model=list(m=2, R='diagonal and equal', A='zero'),
             inits=list(x0=matrix(rep(0,mm),mm,1)), z.score=TRUE, #coef(dfa, type='matrix')$D
             control=list(minit=200, maxit=500, allow.degen=TRUE), silent=2, form='dfa', covariates=rbind(cc,covs))

dfa <- MARSS(y=dat.z, model=list(m=2, R='diagonal and equal', A='zero'),
             inits=coef(dfa, type='matrix'), method='BFGS', z.score=TRUE,
             control=list(maxit=20000), silent=2, form='dfa', covariates=rbind(cc,covs))

par(mfrow=c(5,3))
D_out <- coef(dfa, type='matrix')$D
rownames(D_out) <- colnames(yy)[-1]
for(i in 1:12){
    barplot(D_out[,i], main=month.name[i])
}
for(i in length(cov_choices):1){
    barplot(D_out[,ncol(D_out)-i], main=rev(cov_choices)[i])
}

#get seasonal effects
CC_out = coef(dfa, type="matrix")$C
# The time series of net seasonal effects
seas = CC_out %*% cc[,1:12]
rownames(seas) = 1:mm
colnames(seas) = month.abb
seas

#covariate shiz
dfa$par

# 5 - plot estimated state processes, loadings, and model fits####

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
    ylbl <- names(obs.ts)
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
    ylbl <- names(obs.ts)
    xlbl <- years
    y_ts <- years
    par(mfrow=c(5,2), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
    ymin <- min(dat.z, na.rm=TRUE)
    ymax <- max(dat.z, na.rm=TRUE)
    for(i in 1:nn) {
        lo <- mod_fit$lo[i,]
        mn <- mod_fit$ex[i,]
        up <- mod_fit$up[i,]
        plot(y_ts,mn,xlab="",ylab=ylbl[i],xaxt="n",type="n", cex.lab=1.2,
             ylim=c(ymin,ymax))
        axis(1, at=xlbl, labels=xlbl, cex.axis=1)
        points(y_ts,dat.z[i,], pch=16, col="darkblue")
        lines(y_ts, up, col="darkgray")
        lines(y_ts, mn, col="black", lwd=2)
        lines(y_ts, lo, col="darkgray")
    }
}
fits_plotter()

# model fitting - 6 - tune parameters without covs included ####
#should repeat this in entirety with covs included

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

            dfa <- try(MARSS(y=dat.z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
                             inits=list(x0=0), silent=FALSE,
                             control=list(maxit=20000, allow.degen=TRUE)))
            if(isTRUE(class(dfa)=='try-error')) {next}
        } else {
            print(mmm)

            dfa <- try(MARSS(y=dat.z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
                             inits=list(x0=matrix(rep(0,mmm),mmm,1)), silent=FALSE,
                             control=list(maxit=20000, allow.degen=TRUE)))
            if(isTRUE(class(dfa)=='try-error')) {next}
        }

        #where possible, polish EM estimate with BFGS
        if(RRR != 'equalvarcov'){
            print(paste(RRR,mmm,'BFGS'))

            dfa <- try(MARSS(y=dat.z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZZ, A=aa, D=DD, d=dd, R=RRR),
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
