#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 1/13/17
#plots of best temperature model: temp_due_2m_fixed_at_1978-2015

#before making figures from a particular model run, gotta run through the
#02_testing_or_evaluation.R script and save the image as an .rda file. Instructions are there.

#some of these plots were designed before I made big changes to some of the model designs.
#the ... plots should work unless otherwise noted.

#the "effect size by month" plots are set up to work with all 12 months. If you specified fewer than
#all 12, you'll have to make modifications.

#the "effect size by month over time" plots are set up to work with 4 focal months.
#if you specified a different number, you'll have to make changes.


rm(list=ls()); cat('\014') #clear env and console

# 0 - setup (choose temp/discharge here) ####

setwd('C:/Users/Mike/git/stream_nuts_DFA/manuscript/figures')
setwd('~/git/puget_sound_rivers_DFA/manuscript/figures')

library(viridis)
library(plotrix)
library(stringr)
library(RColorBrewer)

print.letter <- function(label="(a)",xy=c(0.1,0.925),...) {
    tmp <- par("usr")
    text.x <- tmp[1]+xy[1]*diff(tmp[1:2])   #x position, diff=difference
    text.y <- tmp[3]+xy[2]*diff(tmp[3:4])   #y position
    text(x=text.x, y=text.y, labels=label, ...)
}

#load all objects generated during the creation and postprocessing of the best model
# load('temp_due_4m_at_byMo_allMos.rda')
# load('temp_due_4m_at_byMo_acrossTime_may-aug.rda')
# load('temp_due_4m_at_byMo_acrossTime_nov-feb.rda')
# load('temp_due_4m_at_byMo_acrossTime_MASO.rda')
# load('discharge_due_4m_atpc_byMo_allMos.rda')
# load('discharge_due_4m_atpc_byMo_acrossTime_may-aug.rda')
# load('discharge_due_4m_atpc_byMo_acrossTime_nov-feb.rda')
load('discharge_due_4m_atpc_byMo_acrossTime_MASO.rda')

#add percent watershed ice cover data from 2006, average with those from 2011.
#NOTE: WsAreaOver1000 is probably a better metric
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

# 1 - TEMP series plot ####

defpar <- par(mar=c(4,4,4,4))

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
par(defpar)

# 2 - TEMP effect size regression ####

land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
#% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)
# pdf('01_effect_size_reg.pdf', width=7, height=6)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
defpar <- par(mar=c(5,5,2,5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
plot(land$WsAreaOver1000*100, rescaled_effect_size, type='n',
     # plot(land$Ice06_11, rescaled_effect_size, type='n',
     xlab='Watershed area over 1000m (%)', main='',
     # xlab='Watershed % ice cover', main='',
     ylab=expression(paste(Delta, ' water', degree, 'C ',
                           Delta, ' air', degree, 'C'^-1)),
     cex.lab=1.3, cex.axis=1, font=2,
     xaxt='n', xlim=c(0,80))
mod <- lm(rescaled_effect_size ~ I(land$WsAreaOver1000*100))
# mod <- lm(rescaled_effect_size ~ land$Ice06_11)
abline(mod, col='gray', lty=2, lwd=3)
points(land$WsAreaOver1000*100, rescaled_effect_size,
       # points(land$Ice06_11, rescaled_effect_size,
       bg=cols, col='black',
       pch=21, cex=1.5, lwd=2)
# color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
#              rect.col=colorRampPalette(c('brown', 'white'))(10),
#              align='r', gradient='y')
# text(3.39, 0.56, labels='Mean watershed')
# text(3.52, 0.54, labels='elevation (m)')
axis(1, at=seq(0,80,20))
par(defpar)
# dev.off()

# 3 - TEMP loading regression ####

# png('02_loadings_reg.png', width=7, height=6, units='in', res=96, type='cairo')
pdf('02_loadings_reg.pdf', width=7, height=6)
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

# 4 - TEMP water temp by month ####

# png('04_temp_bymonth.png', width=8, height=8, units='in', res=96, type='cairo')
# pdf('04_temp_bymonth.pdf', width=8, height=8)
# seas <- apply(dfa$Estimates$D[,1:12], 2, function(x) x * trans$sds)
seas <- dfa$Estimates$D[,1:12]
seas <- seas+matrix(rep(trans$means,ncol(seas)), ncol=ncol(seas))
defpar <- par(mfrow=c(4,3), oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
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
      side=2, outer=TRUE, line=3)
par(defpar)
# dev.off()

# 5 - TEMP effect size by month ####

# png('03_eff_size_bymonth.png', width=8, height=8, units='in', res=96, type='cairo')
pdf('03a_temp_eff_size_bymonth_vsPctIce.pdf', width=8, height=8)
pdf('03b_temp_eff_size_bymonth_vsWsArea.pdf', width=8, height=8)
# seas <- apply(dfa$Estimates$D[,1:12], 2, function(x) x * trans$sds)

# landcov = land$WsAreaOver1000*100
landcov = land$Ice06_11

csi <- dfa$Estimates$D
csi_names <- rownames(cov_and_seas)
ints <- which(grepl('t\\(covs_z\\).*', csi_names))
mos <- intersect(month.abb, unique(substr(csi_names[ints],
                                          nchar(csi_names[ints][1])-2, nchar(csi_names[ints][1]))))
moInts <- matrix(NA, nrow=nrow(csi), ncol=length(mos), dimnames=list(NULL,mos))
for(i in 1:length(mos)){
    moInts[,i] <- csi[,grepl(paste0('.+', mos[i]), csi_names)]
    # moInts[,i] <- rowSums(csi[,grepl(paste0('.+', mos[i]), csi_names)], na.rm=TRUE)# &
    # grepl(paste0('.+', 'at'), csi_names)], na.rm=TRUE)
}
moInts <- moInts/sd(covs) #restore to original scale
# colnames(moInts)[length(mos)] <- 'All other months'
defpar <- par(mfcol=c(length(mos)/2,2), oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
overall_mean <- mean(moInts)

for(i in 1:length(mos)){
    mod <- lm(moInts[,i] ~ landcov)
    slope <- round(unname(mod$coefficients[2]), 4)
    plot(landcov, moInts[,i], main='', yaxt='n', xaxt='n',
         ylab=paste(month.abb[i]), xlab='', type='n', bty='l',
         ylim=c(min(moInts), max(moInts)))
    polygon(x=c(-5000,5000,5000,-5000), y=c(0,0,overall_mean,overall_mean),
            col='gray80', border=NA)
    abline(mod, col='darkred', lty=2, lwd=2.5)
    points(landcov, moInts[,i], col='black', pch=21,
           cex=1.5, cex.lab=1.3, cex.axis=1, font=2, bg=cols)
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    print.letter(label=substitute(paste(x, ' ', italic('m'), ' = ', y, z),
                                  list(x=colnames(moInts)[i], y=sprintf('%+1.2f', slope), z=sig)),
                 xy=c(0.5,0.9), cex=1.2, font=1, col="black", pos=4)
    if(i %in% c(6,12)) axis(1)
    if(i %in% 1:6) axis(2, las=2)
}
# mtext('Watershed area over 1000m (%)', side=1, outer=TRUE, line=3, font=2)
mtext('Watershed % ice cover', side=1, outer=TRUE, line=3, font=2)
mtext(expression(paste(bold(Delta)~bold('water')~bold(degree)~bold('C')
                       ~bold(Delta)~bold('air')~bold(degree)~bold('C')^-1, sep='')),
      side=2, outer=TRUE, line=3)
par(defpar)
dev.off()

# 6 - TEMP effect size by month over time ####

# pdf('05a_temp_effSize_byMonth_acrossTime_may-aug.pdf', width=14, height=9)
# pdf('05b_temp_effSize_byMonth_acrossTime_nov-feb.pdf', width=14, height=9)
pdf('05c_temp_effSize_byMonth_acrossTime_MASO.pdf', width=14, height=9)

# landcov = land$Ice06_11
landcov = land$WsAreaOver1000*100

csi <- dfa$Estimates$D
csi_names <- rownames(cov_and_seas)
ints <- which(grepl('t\\(covs_z\\).*', csi_names))
mos <- c(intersect(month.abb, unique(str_match(csi_names[ints], '.*?month_fac(.+):.*')[,2])), '0')
sects <- unique(str_match(csi_names[ints], '.*?interval_fac(\\d)')[,2])
moSectInts <- array(NA, dim=c(nrow(csi), length(mos), length(sects)),
                    dimnames=list(rownames(dat_z),mos,sects))
for(i in 1:length(mos)){
    for(j in 1:length(sects)){
        moSectInts[,i,j] <- csi[,which(grepl(paste0('.+', mos[i]), csi_names) &
                                           grepl(paste0(sects[j], '\\z'), csi_names, perl=TRUE))]
    }
}
moSectInts <- moSectInts/sd(covs) #restore to original scale
dimnames(moSectInts)[[2]][length(mos)] <- 'All other months'
defpar <- par(mfcol=c(9,2), oma=c(5,5,1,5), mar=c(0,0,0,0))
snowpal <- colorRampPalette(c('red', 'blue'))
snowcols <- snowpal(10)[as.numeric(cut(landcov, breaks=10))]
# mopalsum <- brewer.pal(9, 'Reds')[c(4,8)]
# mopalwin <- brewer.pal(9, 'Blues')[c(4,8)]
# mopal <- c(mopalsum, mopalwin, 'gray40')
mopal <- brewer.pal(8, 'Set1')[1:5]

overall_mean <- mean(moSectInts)

covOrder <- order(landcov)
covOrder <- c(covOrder[1:9], 99, covOrder[10:17])
cn=0
sct <- dim(moSectInts)[3]
for(i in covOrder){ #rivers
    if(i==99){
        plot(1,1,ann=FALSE,axes=FALSE,type='n')
    } else {
        cn = cn + 1
        plot(1:sct, moSectInts[i,1,], type='n',
             ylim=c(min(moSectInts[,,]), max(moSectInts[,,])),
             bty='n', xaxt='n', yaxt='n')
        polygon(x=c(1,1,sct,sct), y=c(0,overall_mean,overall_mean,0),
                col='gray85', border=NA)
        if(cn %in% c(9,17)) axis(1)
        # if(cn==1) axis(4, line=-1.5, labels=FALSE)
        if(cn==10) axis(2, at=c(0,.4,.8, 1.2), line=-.9, las=2)
        if(cn %in% 1:9) mtext(land$siteCode[i], 2, las=2, col=snowcols[i], font=2)
        if(cn %in% 10:17) mtext(land$siteCode[i], 4, las=2, col=snowcols[i], font=2)
        for(j in 1:dim(moSectInts)[2]){ #months
            lines(1:dim(moSectInts)[3], moSectInts[i,j,], col=mopal[j], lwd=2)
        }
    }
}

# legend(3, 11, legend=c(mos[1:(length(mos)-1)], 'All other months'), lty=1, lwd=2, #may-aug
# legend(3, 20, legend=c(mos[1:(length(mos)-1)], 'All other months'), lty=1, lwd=2, #nov-feb
legend(3, 32, legend=c(mos[1:(length(mos)-1)], 'All other months'), lty=1, lwd=2, #MASO
       col=mopal, horiz=TRUE, xpd=NA, xjust=0.5)
mtext(expression(paste(bold(Delta)~bold('water')~bold(degree)~bold('C')
                       ~bold(Delta)~bold('air')~bold(degree)~bold('C')^-1, sep='')),
      side=2, outer=TRUE, line=3)
mtext('Time series quintile (1978-2015)', side=1, outer=TRUE, line=3, font=2)

par(defpar)
dev.off()
# 7 - TEMP slope comparison ####

# plot the slopes for each monthly series, see if they're different. if so, change river color on the map
# and have it correspond to the line color in the plot. consider a "mean line+sd polygon" approach, wither
# with linear models or with the time series themselves
pdf('09_temp_byMo_overTime.pdf', width=10, height=10, onefile=TRUE)
landcov = land$WsAreaOver1000*100
# landcov = land$Ice06_11
snowpal <- colorRampPalette(c('red', 'blue'))
snowcols <- snowpal(10)[as.numeric(cut(landcov, breaks=10))]
covOrder <- order(landcov)
defpar <- par(mfcol=c(12,2), oma=c(5,5,5,5), mar=c(0,2,0,0))
for(j in 1:12){
    cn=0
    # j=1
    for(i in covOrder){
        cn = cn + 1
        mo <- rep(FALSE,12)
        mo[j] <- TRUE
        mo_ind <- rep(mo, length.out=456)
        plot(1:(456/12), obs_ts[,i][mo_ind], type='l',
             ylim=c(min(obs_ts[mo_ind,],na.rm=T),max(obs_ts[mo_ind,],na.rm=T)),
             xaxt='n', yaxt='n', col='gray60', bty='n')
        if(cn %in% c(12,24)) axis(1, at=c(2,12,22,32), labels=seq(1980,2010,10))
        if(cn == 13) axis(2, las=2)
        abline(h=mean(as.matrix(obs_ts)[mo_ind,],na.rm=T), lty=2, lwd=1, col='gray40')
        mod <- lm(obs_ts[,i][mo_ind] ~ I(1:(456/12)))
        abline(mod, lty=1, col=snowcols[i], lwd=2)
        if(cn %in% 1:12) mtext(land$siteCode[i], 2, las=2, col=snowcols[i], font=2, line=1)
        if(cn %in% 13:24) mtext(land$siteCode[i], 4, las=2, col=snowcols[i], font=2, line=1)
    }
    mtext(paste(month.abb[j], 'stream temp 1978-2015: bluer watersheds have more area over 1000m'),
          3, outer=TRUE, line=2, cex=1, font=2)
    mtext(expression(paste(bold('Water')~bold(degree)~bold('C'), sep='')),
          side=2, outer=TRUE, line=2)
    mtext('Year', side=1, outer=TRUE, line=3, font=2)
}
dev.off()


# 8 - DISCHARGE series plot (needs work) ####

defpar <- par(mar=c(4,4,4,4))

yearsbymo <- substr(yy$date,1,4)
yy <- as.data.frame(apply(yy[,-1], 2, function(i) tapply(i, yearsbymo, mean, na.rm=TRUE)))
pc <- tapply(covs[,1], yearsbymo, mean)
sn <- tapply(covs[,2], yearsbymo, mean)
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
     ylab='Discharge', xaxt='n', xaxs='i', yaxs='i')
for(i in 3:(ncol(yy)-1)){
    lines(yy[,1], yy[,i], type='l', col=colors1[col_ind[i-1]])
    # chili <- readline('>')

    #left off here. plot the climate trend first, then step through each river and note the ones that
    #dont follow suit. designate these with color in the plot

    #also maybe choose something other than grayscale for the background
}
axis(side=1, at=yy[,1][c(T,F)], labels=yy[,1][c(T,F)])
# par(new=T)
# plot(yy[,1], pc, col='darkgreen', lwd=3, type='l', ylim=c(ymin,ymax))
# par(new=T)
# plot(yy[,1], sn, col='blue', lwd=3)
legend(x='topleft', legend=c('Precip'), lwd=3, col=c('red'))
par(defpar)

# 9 - DISCHARGE effect size regression (untouched) ####

land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
#% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)
# pdf('01_effect_size_reg.pdf', width=7, height=6)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
defpar <- par(mar=c(5,5,2,5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
plot(land$WsAreaOver1000*100, rescaled_effect_size, type='n',
     # plot(land$Ice06_11, rescaled_effect_size, type='n',
     xlab='Watershed area over 1000m (%)', main='',
     # xlab='Watershed % ice cover', main='',
     ylab=expression(paste(Delta, ' water', degree, 'C ',
                           Delta, ' air', degree, 'C'^-1)),
     cex.lab=1.3, cex.axis=1, font=2,
     xaxt='n', xlim=c(0,80))
mod <- lm(rescaled_effect_size ~ I(land$WsAreaOver1000*100))
# mod <- lm(rescaled_effect_size ~ land$Ice06_11)
abline(mod, col='gray', lty=2, lwd=3)
points(land$WsAreaOver1000*100, rescaled_effect_size,
       # points(land$Ice06_11, rescaled_effect_size,
       bg=cols, col='black',
       pch=21, cex=1.5, lwd=2)
# color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
#              rect.col=colorRampPalette(c('brown', 'white'))(10),
#              align='r', gradient='y')
# text(3.39, 0.56, labels='Mean watershed')
# text(3.52, 0.54, labels='elevation (m)')
axis(1, at=seq(0,80,20))
par(defpar)
# dev.off()

# 10 - DISCHARGE loading regression (untouched) ####

# png('02_loadings_reg.png', width=7, height=6, units='in', res=96, type='cairo')
pdf('02_loadings_reg.pdf', width=7, height=6)
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

# 11 - DISCHARGE water temp by month (untouched) ####

# png('04_temp_bymonth.png', width=8, height=8, units='in', res=96, type='cairo')
# pdf('04_temp_bymonth.pdf', width=8, height=8)
# seas <- apply(dfa$Estimates$D[,1:12], 2, function(x) x * trans$sds)
seas <- dfa$Estimates$D[,1:12]
seas <- seas+matrix(rep(trans$means,ncol(seas)), ncol=ncol(seas))
defpar <- par(mfrow=c(4,3), oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
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
      side=2, outer=TRUE, line=3)
par(defpar)
# dev.off()

# 12 - DISCHARGE effect size by month ####

# png('03_eff_size_bymonth.png', width=8, height=8, units='in', res=96, type='cairo')
pdf('06a_discharge_eff_size_bymonth_vsPctIce.pdf', width=8, height=8)
# pdf('06b_discharge_eff_size_bymonth_vsWsArea.pdf', width=8, height=8)

# landcov = land$WsAreaOver1000*100
landcov = land$Ice06_11

csi <- dfa$Estimates$D
csi_names <- rownames(cov_and_seas)
ints <- which(grepl('t\\(covs_z\\).*', csi_names))
# ints <- which(grepl('t\\(covs_z\\)p?c?:.*', csi_names))
# mos <- c(intersect(month.abb, unique(str_match(csi_names[ints], '.*?\\](?:.*])?(.+)')[,2])), '0')
mos <- intersect(month.abb, unique(substr(csi_names[ints],
                                          nchar(csi_names[ints][1])-2, nchar(csi_names[ints][1]))))
moInts <- matrix(NA, nrow=nrow(csi), ncol=length(mos), dimnames=list(NULL,mos))
# for(i in 1:length(mos)){
#     moInts[,i] <- rowSums(csi[,grepl(paste0('.+', mos[i]), csi_names)], na.rm=TRUE)
# }
for(i in 1:length(mos)){
    moInts[,i] <- csi[,grepl(paste0('.+', mos[i]), csi_names)]
}
moInts <- moInts/sd(covs[,1]) #restore to original scale (but don't backtransform from log space)
# colnames(moInts)[length(mos)] <- 'All other months'
defpar <- par(mfcol=c(length(mos)/2,2), oma=c(5,5,1,1), mar=c(0.5,0.5,0.5,0.5))
pal <- colorRampPalette(c('brown', 'white'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
overall_mean <- mean(moInts)

for(i in 1:length(mos)){
    mod <- lm(moInts[,i] ~ landcov)
    slope <- round(unname(mod$coefficients[2]), 2)
    plot(landcov, moInts[,i], main='', yaxt='n', xaxt='n',
         ylab=paste(month.abb[i]), xlab='', type='n', bty='l',
         ylim=c(min(moInts), max(moInts)))
    polygon(x=c(-5000,5000,5000,-5000), y=c(0,0,overall_mean,overall_mean),
            col='gray80', border=NA)
    # abline(h=overall_mean, col='royalblue', lwd=2, lty=1)
    abline(mod, col='darkred', lty=2, lwd=2.5)
    points(landcov, moInts[,i], col='black', pch=21,
           cex=1.5, cex.lab=1.3, cex.axis=1, font=2, bg=cols)
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    print.letter(label=substitute(paste(x, '. ', italic('m'), ' = ', y, z),
                                  list(x=colnames(moInts)[i], y=sprintf('%+1.2f', slope), z=sig)),
                 xy=c(0.5,0.9), cex=1.2, font=1, col="black", pos=4)
    if(i %in% c(6,12)) axis(1)
    if(i %in% 1:6) axis(2, las=2)
}
# mtext('Watershed area over 1000m (%)', side=1, outer=TRUE, line=3, font=2)
mtext('Watershed % ice cover', side=1, outer=TRUE, line=3, font=2)
mtext(expression(paste(~italic('ln')~bold(Delta)~bold('Q (cfs) / ')
                       ~bold(Delta)~bold('air')~bold(degree)~bold('C'), sep='')),
      # ~bold(Delta)~bold('precip. (cm)'), sep='')),
      side=2, outer=TRUE, line=3)
par(defpar)
dev.off()

# 13 - DISCHARGE effect size by month over time ####
library(viridis)

# pdf('08a_discharge_effSize_byMonth_acrossTime_may-aug.pdf', width=14, height=9)
# pdf('08b_discharge_effSize_byMonth_acrossTime_nov-feb.pdf', width=14, height=9)
pdf('08c_discharge_effSize_byMonth_acrossTime_MASO.pdf', width=14, height=9)

# landcov = land$Ice06_11
landcov = land$WsAreaOver1000*100

csi <- dfa$Estimates$D
csi_names <- rownames(cov_and_seas)
ints <- which(grepl('t\\(covs_z\\).*', csi_names))
mos <- c(intersect(month.abb, unique(str_match(csi_names[ints], '.*?month_fac(.+):.*')[,2])), '0')
sects <- unique(str_match(csi_names[ints], '.*?interval_fac(\\d)')[,2])
moSectInts <- array(NA, dim=c(nrow(csi), length(mos), length(sects)),
                    dimnames=list(rownames(dat_z),mos,sects))
for(i in 1:length(mos)){
    for(j in 1:length(sects)){
        moSectInts[,i,j] <- csi[,which(grepl(paste0('.+', mos[i]), csi_names) &
                                           grepl(paste0(sects[j], '\\z'), csi_names, perl=TRUE))]
    }
}
moSectInts <- moSectInts/sd(covs[,1]) #restore to original scale
dimnames(moSectInts)[[2]][length(mos)] <- 'All other months'
defpar <- par(mfcol=c(9,2), oma=c(5,5,1,5), mar=c(0,0,0,0))
snowpal <- colorRampPalette(c('red', 'blue'))
snowcols <- snowpal(10)[as.numeric(cut(landcov, breaks=10))]
# mopalsum <- brewer.pal(9, 'Reds')[c(4,8)]
# mopalwin <- brewer.pal(9, 'Blues')[c(4,8)]
# mopal <- c(mopalsum, mopalwin, 'gray40')
mopal <- brewer.pal(8, 'Set1')[1:5]
overall_mean <- mean(moSectInts)
covOrder <- order(landcov)
covOrder <- c(covOrder[1:9], 99, covOrder[10:17])
cn=0
sct <- dim(moSectInts)[3]
for(i in covOrder){ #rivers

    if(i==99){
        plot(1,1,ann=FALSE,axes=FALSE,type='n')
    } else {
        cn = cn + 1
        plot(1:sct, moSectInts[i,1,], type='n',
             ylim=c(min(moSectInts[,,]), max(moSectInts[,,])),
             bty='n', xaxt='n', yaxt='n')
        if(cn %in% 1:9){
            polygon(x=c(.95,.95,sct+2,sct+2), y=c(0,overall_mean,overall_mean,0),
                    col='gray75', border=NA)
        } else {
            polygon(x=c(.8,.8,sct+.05,sct+.05), y=c(0,overall_mean,overall_mean,0),
                    col='gray75', border=NA)
        }
        if(cn %in% c(9,17)) axis(1)
        # if(cn==1) axis(4, line=-1.5, labels=FALSE)
        if(cn==10) axis(2, at=c(-.1,0,.1,.2), line=-1.1, las=2)
        if(cn %in% 1:9) mtext(land$siteCode[i], 2, las=2, col=snowcols[i], font=2)
        if(cn %in% 10:17) mtext(land$siteCode[i], 4, las=2, col=snowcols[i], font=2)
        for(j in 1:dim(moSectInts)[2]){ #months
            lines(1:dim(moSectInts)[3], moSectInts[i,j,], col=mopal[j], lwd=2)
        }
    }

    # cn = cn + 1
    # plot(1:dim(moSectInts)[3], moSectInts[i,1,], type='n',
    #      ylim=c(min(moSectInts[,,]), (max(moSectInts[,,]))),
    #      bty='n', xaxt='n', yaxt='n')
    # if(cn %in% c(12,19)) axis(1)
    # # if(cn==1) axis(4, line=-1.5, labels=FALSE)
    # if(cn==13) axis(2, line=-.7, las=2)
    # if(cn %in% 1:12) mtext(land$siteCode[i], 2, las=2, col=snowcols[i], font=2)
    # if(cn %in% 13:20) mtext(land$siteCode[i], 4, las=2, col=snowcols[i], font=2)
    # for(j in 1:dim(moSectInts)[2]){ #months
    #     lines(1:dim(moSectInts)[3], moSectInts[i,j,], col=mopal[j], lwd=2)
    # }
}

# legend(3, 3.2, legend=c(mos[1:(length(mos)-1)], 'All other months'), lty=1, lwd=2, #may-aug
# legend(3, 2.7, legend=c(mos[1:(length(mos)-1)], 'All other months'), lty=1, lwd=2, #nov-dec
legend(3, 6.2, legend=c(mos[1:(length(mos)-1)], 'All other months'), lty=1, lwd=2, #MASO
       col=mopal, horiz=TRUE, xpd=NA, xjust=0.5)
mtext(expression(paste(italic('ln')~bold(Delta)~bold('Q (cfs) / ')
                       ~bold(Delta)~bold('precip. (cm)'), sep='')),
      side=2, outer=TRUE, line=3)
mtext('Time series quintile (1978-2015)', side=1, outer=TRUE, line=3, font=2)

par(defpar)
dev.off()
# 14 - DISCHARGE slope comparison ####

pdf('10_discharge_byMo_overTime.pdf', width=10, height=10, onefile=TRUE)
landcov = land$WsAreaOver1000*100
# landcov = land$Ice06_11
snowpal <- colorRampPalette(c('red', 'blue'))
snowcols <- snowpal(10)[as.numeric(cut(landcov, breaks=10))]
covOrder <- order(landcov)
defpar <- par(mfcol=c(10,2), oma=c(5,5,5,5), mar=c(0,2,0,0))
for(j in 1:12){
    cn=0
    # j=1
    for(i in covOrder){
        cn = cn + 1
        mo <- rep(FALSE,12)
        mo[j] <- TRUE
        mo_ind <- rep(mo, length.out=456)
        plot(1:(456/12), log(obs_ts[,i][mo_ind]), type='l',
             ylim=c(min(log(obs_ts[mo_ind,]),na.rm=T),max(log(obs_ts[mo_ind,]),na.rm=T)),
             xaxt='n', yaxt='n', col='gray60', bty='n')
        if(cn %in% c(10,19)) axis(1, at=c(2,12,22,32), labels=seq(1980,2010,10))
        if(cn == 11) axis(2, las=2)
        abline(h=mean(log(as.matrix(obs_ts)[mo_ind,]),na.rm=T), lty=2, lwd=1, col='gray40')
        mod <- lm(log(obs_ts[,i][mo_ind]) ~ I(1:(456/12)))
        abline(mod, lty=1, col=snowcols[i], lwd=2)
        if(cn %in% 1:10) mtext(land$siteCode[i], 2, las=2, col=snowcols[i], font=2, line=1)
        if(cn %in% 11:19) mtext(land$siteCode[i], 4, las=2, col=snowcols[i], font=2, line=1)
    }
    plot(1,1,type='n',ann=F,axes=F)
    mtext(paste(month.abb[j], 'discharge 1978-2015: bluer watersheds have more area over 1000m'),
          side=3, outer=TRUE, line=2, cex=1, font=2)
    mtext('log Q (CFS)', font=2,
          side=2, outer=TRUE, line=2)
    mtext('Year', side=1, outer=TRUE, line=3, font=2)
}
dev.off()


# 15 - save stuff for later ####

#save effect sizes for discharge vs. temp plot
# saveRDS(moInts, '../../saved_structures/moInts_discharge_due_4m_at.rds')
# saveRDS(moInts, '../../saved_structures/moInts_temp_due_4m_at.rds')

#save stuff for air water discharge plot
saveRDS(list(covs, obs_ts), '../../saved_structures/air_water_discharge.rds') # from temp_due_4m_at_byMo_allMos.rda


# 16 - DISCHARGE vs. TEMP ####

#load discharge_due_4m_atpc_byMo_allMos.rda in the setup section, or else names
#will be screwed up.

disch_moInts <- readRDS('../../saved_structures/moInts_discharge_due_4m_at.rds')
temp_moInts <- readRDS('../../saved_structures/moInts_temp_due_4m_at.rds')
# "A"  "B"  "C"  "E"  "F"  "G"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "S"  "T"  "U"
# [20] "V"  "W"  "X"  "Z"  "ZA"
# "A"  "B"  "C"  "E"  "F"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "U"  "W"  "Z"  "ZA"

pdf('11a_disch_vs_air_byWsAreaOver1000.pdf', width=8, height=8)
# pdf('11b_disch_vs_air_byWsArea.pdf', width=8, height=8)
pal <- colorRampPalette(c('white', 'black'))
# cols <- pal(20)[as.numeric(cut(land$WsAreaSqKm, breaks=20))]
cols <- pal(10)[as.numeric(cut(land$WsAreaOver1000*100[-c(6,17,18,20,22)], breaks=10))]
defpar <- par(mfcol=c(6,2), oma=c(5,8,1,3), mar=c(0.5,0.5,0.5,0.5))
for(i in 1:12){
    mod <- lm(temp_moInts[,i][-c(6,17,18,20,22)] ~ disch_moInts[,i])
    plot(disch_moInts[,i], temp_moInts[,i][-c(6,17,18,20,22)], cex=1.5,
         # pch=names(obs_ts),
         yaxt='n', xaxt='n', bty='l', col='black', bg=cols, pch=21,
         xlim=c(min(disch_moInts), max(disch_moInts)),
         ylim=c(min(temp_moInts[-c(6,17,18,20,22)]), max(temp_moInts[-c(6,17,18,20,22)])))
    abline(mod, col='darkred', lty=2, lwd=2)
    if(i %in% 1:6) {axis(2, las=2); mtext(month.abb[i], 2, line=3, las=2)}
    if(i %in% 7:12) mtext(month.abb[i], 4, line=1, las=2)
    if(i %in% c(6,12)) axis(1)
}
mtext(expression(paste(~italic('ln')~bold(Delta)~bold('Q (cfs) / ')
                       ~bold(Delta)~bold('air')~bold(degree)~bold('C'), sep='')),
      # ~bold(Delta)~bold('precip. (cm)'), sep='')),
      side=1, outer=TRUE, line=3)
mtext(expression(paste(bold(Delta)~bold('water')~bold(degree)~bold('C')
                       ~bold(Delta)~bold('air')~bold(degree)~bold('C')^-1, sep='')),
      side=2, outer=TRUE, line=6)
par(defpar)
dev.off()

# 17 - air water discharge

#load temp_due_4m_at_byMo_allMos.rda for this one
awd <- readRDS('../../saved_structures/air_water_discharge.rds')

pdf('12_air_water_discharge.pdf', width=8, height=6)
defpar <- par(mar=c(5,4,1,6))
airtemps = watertemps = discharge = numeric()
for(i in 1:12){
    airtemps[i] <- mean(awd[[1]][seq(from=i, by=12, to=nrow(awd[[1]]))])
    watertemps[i] <- mean(as.matrix(awd[[2]][seq(from=i, by=12, to=nrow(awd[[1]])),]), na.rm=TRUE)
    discharge[i] <- mean(as.matrix(obs_ts[seq(from=i, by=12, to=nrow(awd[[1]])),]), na.rm=TRUE)
}
plot(1:12, airtemps, xaxt='n', xlab=expression(paste(~bold('Month'))),
     ylab=expression(paste(~bold('Mean temperature')~bold(degree)~bold('C'))),
     las=2)
points(1:12, watertemps, pch=20)
axis(1, at=1:12, labels=month.abb)
legend('left', legend=c('Air', 'Water', 'Discharge'), pch=c(1,20,17), bty='n',
       col=c('black','black','gray60'))
par(new=T)
plot(1:12, discharge, type='b', pch=17, col='gray60', xaxt='n', yaxt='n', xlab='', ylab='')
axis(4, las=1)
mtext('Mean discharge (CFS)', 4, font=2, las=3, line=3)
par(defpar)
dev.off()
