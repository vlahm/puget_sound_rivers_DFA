# library(MARSS)
library(viridis)

mm=2
nts=30
x = 1:120
# y = z = rep(c(1,1.4,3,1.6),per)
# z[seq(1,120,8)] <- 1.5
# y = c(1,2,3,4,5,5,4,3,3,1,1,2,4,5,6,4,3,2,1,2,2,3,5,5,4,3,2,1)
# z = c(1,1,3,4,4,5,4,3,2,1,1,2,3,5,5,4,2,2,1,2,3,4,4,5,4,3,1,1)
# w = c(1,2,3,4,4,5,4,3,2,1,2,2,3,4,5,5,3,2,1,1,3,4,5,5,4,3,2,2)
# x = 1:length(y)
y = z = w = rep(NA,120)
for(i in 1:30){
    w[(i*4-3):(i*4)] <- c(rnorm(1,0,0.1), rnorm(1,0.4,0.1), rnorm(1,1,0.1), rnorm(1,0.6,0.1))
    y[(i*4-3):(i*4)] <- c(rnorm(1,0,0.1), rnorm(1,0.4,0.1), rnorm(1,1,0.1), rnorm(1,0.6,0.1))
    z[(i*4-3):(i*4)] <- c(rnorm(1,0,0.1), rnorm(1,0.4,0.1), rnorm(1,1,0.1), rnorm(1,0.6,0.1))
}
for(i in 1:120){
    w[i] <- w[i]+0.01*i
    y[i] <- y[i]+0.01*i
    z[i] <- z[i]-0.01*i
}
# for(i in 1:120){
#     w[i] <- w[i]-0.02*i
#     y[i] <- y[i]-0.02*i
#     z[i] <- z[i]-0.02*i
# }

data <- t(scale(cbind(y,z,w)))
rownames(data) <- c('y','z','w')

cov <- t(scale(seq(1.2,-1.2,length.out=length(y))))

plot(x,data[1,], type='l', col='red')
lines(x,data[2,],col='blue')
lines(x,data[3,],col='green')

ccgen <- function(per, nyr, fourier){
    if(fourier==F){
        year_block <- diag(per)
        cc <- Reduce(function(x,y) {cbind(x,y)},
                     eval(parse(text=paste0('list(',
                                            paste(rep('year_block', nyr,), sep=', ', collapse=', '),
                                            ')'))))
    } else{
        cos_t = cos(2 * pi * seq(dim(data)[2]) / per)
        sin_t = sin(2 * pi * seq(dim(data)[2]) / per)
        cc = rbind(cos_t,sin_t)
    }
    return(cc)
}
cc <- ccgen(4, nts, F)



# dfa <- MARSS(y=data, model=list(m=mm, R='diagonal and equal', A='zero'),
#              inits=list(x0='zero'), z.score=T,
#              control=list(minit=1, maxit=500, allow.degen=T, safe=T, trace=1),
#              # silent=2, form='dfa')
#              silent=2, form='dfa',
#              # covariates=cc)
#              covariates=cov)
#              # covariates=rbind(cc,cov))

dfa2 <- runDFA(obs=data, NumStates=2, ErrStruc='DUE', EstCovar=T, Covars=cov)
dfa3 <- runDFA(obs=data, NumStates=1, ErrStruc='DUE', EstCovar=T, Covars=cov)

plot(dfa2$Estimates$u[1,], type='l')
lines(dfa2$Estimates$u[2,])
lines(dfa3$Estimates$u[1,], type='l', col='red')
lines(dfa3$Estimates$u[2,], col='red')

dfa2$Fits
#plot data series 1 and then fit on top, e.g.

dfa2$Estimates$Z #loadings



#plot
Z_est <- coef(dfa, type="matrix")$Z
H_inv <- varimax(Z_est)$rotmat
Z_rot = Z_est %*% H_inv
proc_rot = solve(H_inv) %*% dfa$states

process_plotter <- function(){
    defpar <- par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(mm, 1))
    xlbl = y_ts =  1:length(y)
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
    par(defpar)
}
process_plotter()

loading_plotter <- function(){
    par(mai=c(0.5,0.5,0.5,0.1), omi=c(0,0,0,0), mfrow=c(mm, 1))
    ylbl <- colnames(t(data))
    clr <- viridis(nrow(data)) #colors may not line up with series plots in section 2
    ylm <- c(-1,1)*max(abs(proc_rot))
    minZ <- 0
    ylm <- c(-1,1)*max(abs(Z_rot))
    for(i in 1:mm) {
        plot(c(1:nrow(data))[abs(Z_rot[,i])>minZ], as.vector(Z_rot[abs(Z_rot[,i])>minZ,i]), type="h",
             lwd=2, xlab="", ylab="", xaxt="n", ylim=ylm, xlim=c(0.5,nrow(data)+0.5), col=clr)
        for(j in 1:nrow(data)) {
            if(Z_rot[j,i] > minZ) {text(j, -0.03, ylbl[j], srt=90, adj=1, cex=1.2, col=clr[j])}
            if(Z_rot[j,i] < -minZ) {text(j, 0.03, ylbl[j], srt=90, adj=0, cex=1.2, col=clr[j])}
            abline(h=0, lwd=1.5, col="gray")
        }
        mtext(paste("Factor loadings on process",i),side=3,line=0.5)
    }
}
loading_plotter()

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
fits_plotter <- function(){
    ylbl <- colnames(t(data))
    xlbl = y_ts =  1:length(x)
    # par(mfrow=c(1,1), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
    par(mfrow=c(nrow(data),1), mai=c(0.6,0.7,0.1,0.1), omi=c(0,0,0,0))
    ymin <- min(data, na.rm=TRUE)
    ymax <- max(data, na.rm=TRUE)
    for(i in 1:nrow(data)) {
        lo <- mod_fit$lo[i,]
        mn <- mod_fit$ex[i,]
        up <- mod_fit$up[i,]
        plot(y_ts,mn,xlab="",ylab=ylbl[i],xaxt="n",type="n", cex.lab=1.2,
             ylim=c(ymin,ymax))
        axis(1, at=xlbl, labels=xlbl, cex.axis=1)
        points(y_ts,data[i,], pch=16, col="darkblue")
        lines(y_ts, up, col="darkgray")
        lines(y_ts, mn, col="black", lwd=2)
        lines(y_ts, lo, col="darkgray")
    }
}
fits_plotter()

D_out <- coef(dfa, type='matrix')$D
for(i in 1:4){
    barplot(D_out[,i], main=c('spr','sum','aut','win')[i])
}
for(i in length(cov_choices):1){
    barplot(D_out[,ncol(D_out)-i], main=rev(cov_choices)[i])
}
