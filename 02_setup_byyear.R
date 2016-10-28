#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 8/10/2016
#last edit: 8/11/2016

# setup - 0 - CHOOSE yy maxes or means ####
load("C:/Users/Mike/git/stream_nuts_DFA/data/chemPhys_data/yys_byyear_mean.rda")
meantemp <- read.csv("C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_year/ meantemp3 .csv")
precip <- read.csv("C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_year/ precip3 .csv")
library(MARSS)
library(viridis)
# if (is.null(dev.list()) == TRUE){windows(record=TRUE)} #open new plot window unless already open

# pdf("C:/Users/Mike/Desktop/Grad/Projects/Thesis/stream_nuts_DFA/time_series_maxes.pdf", onefile=T)

# 1 - CHOOSE and subset response and covs####
# response choices: COND FC NH3_N NO2_NO3 OP_DIS OXYGEN PH PRESS SUSSOL TEMP TP_P TURB
y_choice <- 'TEMP'
# cov choices: meantemp meantemp_anom precip precip_anom hydroDrought hydroDrought_anom meteoDrought
# meteoDrought_anom ZDrought ZDrought_anom
cov_choices <- c('meantemp')

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
startyr <- 1978
endyr <- 2014
yy <- subsetter(yy, start=startyr, end=endyr, na_thresh=0.7)

names(yy)
years <- as.numeric(format(yy[,1], '%Y')) #could use as dates if necessary, as integer years for now
obs.ts <- yy[,-1]

covs <- climate[climate$Date >= startyr & climate$Date <= endyr, colnames(climate) %in% cov_choices]

# 2 - plot response variable time series ####
par(mar=c(4,4,4,4))
series_plotter <- function(){
    colors1 <- viridis(ncol(yy)-1, end=1)
    ymin <- min(yy[,-1], na.rm=TRUE)
    ymax <- max(yy[,-1], na.rm=TRUE)

    plot(years, yy[,2], type='l', col=colors1[1], ylim=c(ymin,ymax), xlab='Time',
         # main=paste(y_choice),
         ylab='Turbidity (NTU)', xaxt='n', xaxs='i', yaxs='i')
    for(i in 3:ncol(yy)){
        lines(years, yy[,i], type='l', col=colors1[i-1])
    }
    axis(side=1, at=years[c(T,F)], labels=years[c(T,F)])
}
series_plotter()

# fit Lowess line
xxx <- as.numeric(format(rep(yy[,1], ncol(yy)-1), '%Y'))
yyy <- as.vector(as.matrix(yy[,-1]))
xxx <- xxx[which(!is.na(yyy))]
yyy <- yyy[!is.na(yyy)]
agg <- aggregate(yyy, by=list(xxx), FUN=mean)
loess_out <- loess(agg$x ~ agg$Group.1, span=.2, degree=2)
lines(agg$Group.1, loess_out$fit, col='blue', lwd=3)
#overlay air temp trend
# lines(agg$Group.1, meantemp$at[which(as.numeric(substr(meantemp$date,1,4)) %in% agg$Group.1)],
#       col='red', lwd=3)
# legend(x='bottomleft', legend=c('Water Temp', 'Air temp'), lwd=3, col=c('blue','red'))
#overlay precip trend
par(new=TRUE)
plot(agg$Group.1, precip$pc[which(as.numeric(substr(precip$date,1,4)) %in% agg$Group.1)],
      col='red', lwd=3, type='l', axes=F, xlab='', ylab='')
legend(x='topleft', legend=c('Mean Turbidity', 'Precip'), lwd=3, col=c('blue','red'))
axis(side=4)
mtext(text='Precipitation (cm)', side=4, line=2.5)


# 3 - Set up input matrics to MARSS function call ####
# scale and center y and cov data
dat.z <- t(scale(as.matrix(obs.ts)))
covs <- t(scale(as.matrix(covs)))

# state equation params
mm <- 2 #number of hidden processes (trends)
BB <- "identity"  #'BB' is identity: 1's along the diagonal & 0's elsewhere (not part of DFA)
uu <- "zero"  # 'uu' is a column vector of 0's (not part of DFA)
CC <- "zero"  # 'CC' and 'cc' are for covariates in the state process equation, which we do not have
cc <- "zero"  # i.e. matrix(0,1,wk_last)
QQ <- "identity"  # 'QQ' is identity (would usually be tuned if we were doing actual marss)

# observation equation params
nn <- length(names(obs.ts)) # number of obs time series
aa <- "zero" # 'aa' is the offset/scaling (zero allowed because we standardize the data)
             #(we want zero because we're not interested in what the offsets actually are)
DD <- "unconstrained"  # 'DD' and 'd' are for covariates
dd <- covs
RR <- "equalvarcov" # 'RR' is var-cov matrix for obs errors (tune this)
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
dfa <- MARSS(y=dat.z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZ, A=aa, D=DD, d=dd, R=RR),
             inits=list(x0=matrix(rep(0,mm),mm,1)),
             control=list(maxit=2000, allow.degen=TRUE), silent=FALSE)
# dfa <- MARSS(y=dat.z, model=list(B=BB, U=uu, C=CC, c=cc, Q=QQ, Z=ZZ, A=aa, D=DD, d=dd, R=RR),
#               inits=dfa$par,
#               control=list(maxit=3000), method='BFGS') #can't use BFGS for equalvarcov
dfa <- MARSS(y=dat.z, model=list(m=2, R='diagonal and equal', A='zero'),
             inits=list(x0='zero'), z.score=TRUE,
             control=list(minit=200, maxit=500, allow.degen=TRUE), silent=2, form='dfa')

names(dfa)

# 5 - plot estimated state processes, loadings, and model fits####

# varimax rotation to get Z loadings
Z_est <- coef(dfa, type="matrix")$Z # get the estimated ZZ
H_inv <- varimax(Z_est)$rotmat # get the inverse of the rotation matrix
Z_rot = Z_est %*% H_inv # rotate factor loadings
proc_rot = solve(H_inv) %*% dfa$states # rotate processes

# plot hidden processes
process_plotter <- function(){
    par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(mm, 1))
    xlbl <- years
    y_ts <- years
    ylm <- c(-1,1)*max(abs(proc_rot))
    for(i in 1:mm) {
      plot(y_ts,proc_rot[i,], type="n", bty="L",
           ylim=ylm, xlab="", ylab="", xaxt="n")
      abline(h=0, col="gray")
      lines(y_ts,proc_rot[i,], lwd=2)
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
        saveRDS(dfa, file=paste0("C:/Users/vlahm/Desktop/stream_nuts_DFA/stream_nuts_DFA/model_objects/",
                                 R_names[R_strucs == RRR], '_', mmm, 'm_', startyr, y_choice, '.rds'))

        if(mmm > 1){

            #open plot device
            pdf(file=paste0("C:/Users/vlahm/Desktop/stream_nuts_DFA/stream_nuts_DFA/model_outputs/",
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
write.csv(model_out, file=paste0("C:/Users/vlahm/Desktop/stream_nuts_DFA/stream_nuts_DFA/model_objects/",
                         'param_tuning_dataframe_', startyr, y_choice, '.csv'))
