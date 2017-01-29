#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 1/13/17
#plots of best temperature model: temp_due_2m_fixed_at_1978-2015

rm(list=ls()); cat('\014') #clear env and console

# 0 - setup part 1 ####
# load("C:/Users/Mike/git/stream_nuts_DFA/data/chemPhys_data/yys_byyear_mean.rda")
# meantemp <- read.csv("C:/Users/Mike/git/stream_nuts_DFA/data/climate_data/by_year/obsolete/meantemp3 .csv")
#
# library(viridis)
#
# #open plot window if not already open
# if (is.null(dev.list()) == TRUE){
#     if(.Platform$OS.type == "windows"){
#         windows(record=TRUE, width=16, height=9)
#     } else {
#         x11(width=16, height=9)
#     }
# }
#
# #prepare by-year temp data
# y_choice <- 'TEMP'
# cov_choices <- c('meantemp')
# startyr <- 1978
# endyr <- 2014
#
# yy <- eval(parse(text=y_choice))
#
# # subset by year and exclude columns with >= na_thresh proportion of NAs
# subsetter <- function(yy, start, end, na_thresh=1){
#
#     #extract year subset
#     years <- as.numeric(format(yy$date, '%Y'))
#     yy <- yy[years >= start & years <= end, ]
#
#     #extract subset of columns with <= na_thresh proportion of NAs
#     col_na_prop <- apply(yy, 2, function(x) sum(is.na(x)) / length(x))
#     yy <- yy[,col_na_prop <= na_thresh]
#
#     rownames(yy) <- 1:nrow(yy)
#
#     return(yy)
# }
# yy <- subsetter(yy, start=startyr, end=endyr, na_thresh=0.55)
#
# years <- as.numeric(format(yy[,1], '%Y')) #could use as dates if necessary, as integer years for now
# obs.ts <- yy[,-1]
#
# covs <- climate[climate$Date >= startyr & climate$Date <= endyr, colnames(climate) %in% cov_choices]
#
# # 2 - plot response variable time series
# par(mar=c(4,4,4,4))
# colors1 <- viridis(ncol(yy)-1, end=1)
# ymin <- min(yy[,-1], na.rm=TRUE)
# ymax <- max(yy[,-1], na.rm=TRUE)
#
# plot(years, yy[,2], type='l', col=colors1[1], ylim=c(ymin,ymax), xlab='Time',
#      # main=paste(y_choice),
#      ylab='TEMP', xaxt='n', xaxs='i', yaxs='i')
# for(i in 3:ncol(yy)){
#     lines(years, yy[,i], type='l', col=colors1[i-1])
# }
# axis(side=1, at=years[c(T,F)], labels=years[c(T,F)])
#
# #overlay time series average
# xxx <- as.numeric(format(rep(yy[,1], ncol(yy)-1), '%Y'))
# yyy <- as.vector(as.matrix(yy[,-1]))
# xxx <- xxx[which(!is.na(yyy))]
# yyy <- yyy[!is.na(yyy)]
# agg <- aggregate(yyy, by=list(xxx), FUN=mean)
# # loess_out <- loess(agg$x ~ agg$Group.1, span=.2, degree=2)
# # lines(agg$Group.1, loess_out$fit, col='blue', lwd=3)
#
# #overlay air temp trend
# lines(agg$Group.1, meantemp$at[which(as.numeric(substr(meantemp$date,1,4)) %in% agg$Group.1)],
#       col='red', lwd=3)
# legend(x='bottomleft', legend=c('Air temp'), lwd=3, col=c('red'))
# # legend(x='bottomleft', legend=c('Water Temp', 'Air temp'), lwd=3, col=c('blue','red'))

# 0 - setup part 2 ####

setwd('C:/Users/Mike/git/stream_nuts_DFA/manuscript/figures')
setwd('~/git/puget_sound_rivers_DFA/manuscript/figures')

library(viridis)
library(plotrix)

print.letter <- function(label="(a)",xy=c(0.1,0.925),...) {
    tmp <- par("usr")
    text.x <- tmp[1]+xy[1]*diff(tmp[1:2])   #x position, diff=difference
    text.y <- tmp[3]+xy[2]*diff(tmp[3:4])   #y position
    text(x=text.x, y=text.y, labels=label, ...)
}

#load all objects generated during the creation and postprocessing of the best model
load('Renv_image3.rda')
#this is the "best temp" dfa from 02_testing_or_evaluation as of 1/28/17
dfa <- readRDS('../manuscript/figures/other_dfa.rds')
# dfa <- readRDS('other_dfa.rds')

#add percent watershed ice cover data from 2006, average with those from 2011
ice2006 <- read.csv('PctIce2006Ws.csv', stringsAsFactors=FALSE)
land <- merge(land,ice2006[,2:3], by.x='siteCode', by.y='site')
land$Ice06_11 <- rowMeans(cbind(land$PctIce2006Ws, land$PctIce2011Ws))

#open plot window if not already open
if (is.null(dev.list()) == TRUE){
    if(.Platform$OS.type == "windows"){
        windows(record=TRUE, width=16, height=9)
    } else {
        x11(width=16, height=9)
    }
}

# 1 - series plot ####

par(mar=c(4,4,4,4))

yearsbymo <- substr(yy$date,1,4)
yy <- as.data.frame(apply(yy[,-1], 2, function(i) tapply(i, yearsbymo, mean, na.rm=TRUE)))
covs <- tapply(covs, yearsbymo, mean)
library(data.table)
yy <- as.data.frame(setDT(yy, keep.rownames = TRUE)[])
yy[is.na(yy)] <- NA #turn NaNs into NAs

ymin <- min(yy[,-1], na.rm=TRUE) - .2
ymax <- max(yy[,-1], na.rm=TRUE)

# colors1 <- viridis(ncol(yy)-1, end=1)
palette <- colorRampPalette(colors=c("gray90", "black"))
colors1 <- palette(24)
col_ind <- rep_len(1:24, length.out=24)

plot(yy[,1], yy[,2], type='l', col=colors1[col_ind[1]], ylim=c(ymin,ymax), xlab='Time',
     ylab='TEMP', xaxt='n', xaxs='i', yaxs='i')
for(i in 3:(ncol(yy)-1)){
    lines(yy[,1], yy[,i], type='l', col=colors1[col_ind[i-1]])
    # chili <- readline('>')

    #left off here. plot the climate trend first, then step through each river and note the ones that
    #dont follow suit. designate these with color in the plot

    #also maybe choose something other than grayscale for the background
}
axis(side=1, at=yy[,1][c(T,F)], labels=yy[,1][c(T,F)])

lines(yy[,1], covs, col='red', lwd=3)
legend(x='bottomleft', legend=c('Air temp'), lwd=3, col=c('red'))


# 2 - effect size regression ####

land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
#% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)
# pdf('01_effect_size_reg.pdf', width=7, height=6)
defpar <- par(mar=c(5,5,2,5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
plot(land$Ice06_11, rescaled_effect_size, type='n',
     xlab='Watershed % ice cover', main='',
     ylab=expression(paste(Delta, ' stream', degree, 'C ',
                           Delta, ' air', degree, 'C'^-1)),
     cex.lab=1.3, cex.axis=1, font=2)
mod <- lm(rescaled_effect_size ~ land$Ice06_11)
abline(mod, col='gray', lty=2, lwd=3)
points(land$Ice06_11, rescaled_effect_size,
     bg=cols, col='black',
     pch=21, cex=1.5, lwd=2)
color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
             rect.col=colorRampPalette(c('brown', 'white'))(10),
             align='r', gradient='y')
text(3.39, 0.56, labels='Mean watershed')
text(3.52, 0.54, labels='elevation (m)')
par(defpar)
# dev.off()

# 3 - loading regression ####

pdf('02b_loadings_reg.pdf', width=7, height=6)
defpar <- par(mar=c(5,5,2,5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
plot(land$Ice06_11, dfa$Estimates$Z[,1],
     xlab='Watershed % ice cover', ylab='Trend 1 loadings', type='n',
     main='', cex.lab=1.3, cex.axis=1, font=2)
abline(mod, col='gray', lty=2, lwd=3)
points(land$Ice06_11, dfa$Estimates$Z[,1], col='black', bg=cols,
       pch=21, cex=1.5, lwd=2)
mod <- lm(dfa$Estimates$Z[,1] ~ land$Ice06_11)
color.legend(xl=4,xr=4.4,yb=1.32, yt=1.82, legend=c('147', '1349'),
             rect.col=colorRampPalette(c('brown', 'white'))(10),
             align='r', gradient='y')
text(3.39, 1.61, labels='Mean watershed')
text(3.52, 1.53, labels='elevation (m)')
par(defpar)
dev.off()

# 4 - effect size by month ####

pdf('03_eff_size_bymonth.pdf', width=8, height=8)
# seas <- apply(dfa$Estimates$D[,1:12], 2, function(x) x * trans$sds)
seas <- dfa$Estimates$D[,1:12]
seas <- seas+matrix(rep(trans$means,ncol(seas)), ncol=ncol(seas))
defpar <- par(mfrow=c(4,3), oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$BFIWs, breaks=10))]
overall_mean <- mean(seas)
for(i in 1:12){
    mod <- lm(seas[,i] ~ land$Ice06_11)
    slope <- round(unname(mod$coefficients[2]), 2)
    plot(land$Ice06_11, seas[,i], main='', yaxt='n', xaxt='n',
         ylab=paste(month.abb[i]), xlab='', type='n', bty='l',
         ylim=c(min(seas), max(seas)))
    abline(h=overall_mean, col='royalblue', lwd=2, lty=1)
    abline(mod, col='gray40', lty=2, lwd=2.5)
    points(land$Ice06_11, seas[,i], col='black', pch=21,
           cex=1.5, cex.lab=1.3, cex.axis=1, font=2, bg=cols)
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.05, '*', '')
    print.letter(label=substitute(paste(x, '. ', italic('m'), ' = ', y, z),
                                  list(x=month.abb[i], y=sprintf('%+1.2f', slope), z=sig)),
                 xy=c(0.5,0.9), cex=1.2, font=1, col="black", pos=4)
    if(i %in% 10:12) axis(1)
    if(i %in% c(1,4,7,10)) axis(2, las=2)
}
mtext('Watershed % ice cover', side=1, outer=TRUE, line=3, font=2)
mtext(expression(paste(bold('Mean monthly water temp (')~bold(degree)~bold('C)'))),
                     # ~bold(Delta)~bold('air')~bold(degree)~bold('C')^-1, sep='')),
      side=2, outer=TRUE, line=3)
par(defpar)
dev.off()

# 5 - slope comparison ####

# plot the slopes for each monthly series, see if they're different. if so, change river color on the map
# and have it correspond to the line color in the plot. consider a "mean line+sd polygon" approach, wither
# with linear models or with the time series themselves

