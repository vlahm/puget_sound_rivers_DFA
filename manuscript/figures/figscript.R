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
# load('discharge_due_4m_atpc_byMo_acrossTime_MASO.rda')
# load('discharge_due_5m_atpcsn_byMo_allMos.rda')
load('temp_due_5m_atpcsn_byMo_allMos.rda')
# load('../../single_trend_exploration/2trendNoSeasNoSnow.rda')

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

# 2 - TEMP effect size regression (single, old) ####

# land_sub <- land[,landcols] #subset landscape variables by those used in the analysis
#
# #% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
# #% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)
# # pdf('01_effect_size_reg.pdf', width=7, height=6)
# # png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# defpar <- par(mar=c(5,5,2,5))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
# plot(land$WsAreaOver1000*100, rescaled_effect_size, type='n',
#      # plot(land$Ice06_11, rescaled_effect_size, type='n',
#      xlab='Watershed area over 1000m (%)', main='',
#      # xlab='Watershed % ice cover', main='',
#      ylab=expression(paste(Delta, ' water', degree, 'C ',
#                            Delta, ' air', degree, 'C'^-1)),
#      cex.lab=1.3, cex.axis=1, font=2,
#      xaxt='n', xlim=c(0,80))
# mod <- lm(rescaled_effect_size ~ I(land$WsAreaOver1000*100))
# # mod <- lm(rescaled_effect_size ~ land$Ice06_11)
# abline(mod, col='gray', lty=2, lwd=3)
# points(land$WsAreaOver1000*100, rescaled_effect_size,
#        # points(land$Ice06_11, rescaled_effect_size,
#        bg=cols, col='black',
#        pch=21, cex=1.5, lwd=2)
# color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
#              rect.col=colorRampPalette(c('brown', 'white'))(10),
#              align='r', gradient='y')
# text(3.39, 0.56, labels='Mean watershed')
# text(3.52, 0.54, labels='elevation (m)')
# axis(1, at=seq(0,80,20))
# par(defpar)
# dev.off()
# 1.1 - TEMP effect size regression (linked with loading regression) ####

# system('taskkill /f /im AcroRd32.exe')
# pdf('16_temp_all_reg.pdf', width=7.5, height=7.5)
layout(matrix(c(1:6,9,7,8),nrow=3,byrow=TRUE))
defpar = par(oma=c(0,0,1,1), mar=c(4,3.5,0,0))
# landvar=land$WsAreaOver1000/land$WsAreaSqKm #??
landvar=land$ElevWs/100
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
#% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# defpar <- par(oma=c(2,4,2,1), mar=c(4,2,2,1))
# pal <- colorRampPalette(c('sienna2', 'dodgerblue4'))
# brks = as.numeric(cut(land$Ice06_11, breaks=10))
# cols = pal(length(unique(brks)))
# cols = as.character(factor(brks, labels=cols))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
for(covr in 1:3){
    res = rescaled_effect_size[,covr]
    if(covr==3){res = res * 2.54}
    plot(landvar, res, type='n', yaxt='n',
         # xlab='Watershed area over 1000m (%)', main='',
         xlab='', main='', ylab='', cex.lab=1.3, cex.axis=1, font=2, xaxt='n',bty='l')
    mod <- lm(res ~ landvar)
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    sig = NULL #turning off sig plotting
    if(covr == 1){
        mtext(bquote(textstyle(bold(Delta)*bold(T[water])~bold(Delta)*bold(T[air]^-1))~
                         scriptstyle(plain('(')*plain(degree)*plain('C')~plain(degree)*plain(C^-1)*plain(')'))),
              2, cex=.8, line=1.5)
        abline(mod, col='steelblue', lty=2, lwd=3)
        axis(2, padj=.9, tck=-.02)
    }
    if(covr == 2){
        mtext(bquote(textstyle(bold(Delta)*bold(T[water])~bold(Delta)*bold(precip^-1))~
                         textstyle(plain('(')*plain(degree)*plain('C')~plain(cm^-1)*plain(')'))),
              2, cex=.8, line=1.5)
        axis(2, at=c(-.02,0,.02,.04), labels=c(-0.02,0.00,0.02,0.04), padj=.9, tck=-.02)
        # mtext(bquote(textstyle(bold('Perennial watershed ice/snow coverage'))~textstyle(plain('(%)'))), side=1, cex=.8, line=1.7)
        mtext(bquote(textstyle(bold('Mean watershed elevation'))~textstyle(plain('(100 m)'))), side=1, cex=.8, line=1.7)
        abline(mod, col='steelblue', lty=2, lwd=3)
    }
    if(covr == 3){
        mtext(bquote(textstyle(bold(Delta)*bold(T[water])~bold(Delta)*bold(snowmelt^-1))~
                         textstyle(plain('(')*plain(degree)*plain('C')~plain(cm^-1)*plain(')'))), 2, cex=.8, line=1.5)
        # abline(mod, col='steelblue', lty=2, lwd=3)
        axis(2, padj=.9, tck=-.02)
    }
    print.letter(paste(letters[covr],sig), c(.9,.9), cex=1.8, font=2, col='steelblue')
    points(landvar, res,
           bg=cols, col='black', cex=rescale(log(land$WsAreaSqKm),c(1.2,3.2)),
           pch=21,
           # pch=land$siteCode,
           lwd=1)
    # color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
    #              rect.col=colorRampPalette(c('brown', 'white'))(10),
    #              align='r', gradient='y')
    # text(3.39, 0.56, labels='Mean watershed')
    # text(3.52, 0.54, labels='elevation (m)')
    axis(1, padj=-.9, tck=-.02)
    # if(covr %in% c(1,3)){axis(2)}else{axis(2, )}
}
lines(x=c(-19.9,13.8), y=c(-.342,-.342), xpd=NA, lwd=2, col='gray70')
# mtext(expression(paste(Delta,'Q ', Delta, theta^-1)), side=2, outer=TRUE, cex=1.3, line=1)
# par(defpar)
# dev.off()

#for results (gotta create res above):
summary(mod)
mean(res[elev_lo]); sd(res[elev_lo])
mean(res[elev_med]); sd(res[elev_med])
mean(res[elev_hi]); sd(res[elev_hi])

# 1.2 - TEMP loading regression ####

# png('02_loadings_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# pdf('02_loadings_reg.pdf', width=7, height=6)
# defpar <- par(mar=c(5,5,2,5))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
# plot(land$Ice06_11, dfa$Estimates$Z[,1],
#      xlab='Watershed % ice cover', ylab='Trend 1 loadings', type='n',
#      main='', cex.lab=1.3, cex.axis=1, font=2)
# abline(mod, col='gray', lty=2, lwd=3)
# points(land$Ice06_11, dfa$Estimates$Z[,1], col='black', bg=cols,
#        pch=21, cex=1.5, lwd=2)
# mod <- lm(dfa$Estimates$Z[,1] ~ land$Ice06_11)
# color.legend(xl=4,xr=4.4,yb=1.32, yt=1.82, legend=c('147', '1349'),
#              rect.col=colorRampPalette(c('brown', 'white'))(10),
#              align='r', gradient='y')
# text(3.39, 1.61, labels='Mean watershed')
# text(3.52, 1.53, labels='elevation (m)')
# par(defpar)
# dev.off()
plot(1,1,type='n',ann=FALSE,axes=FALSE)
legend(x=1,y=1.505, legend=c('Rain-dominated','Rain-and-snow','Snow-dominated'),
       xpd=NA, pt.bg=c('black','gray75','white'), pch=21, xjust=0.5,
       col='black', cex=1.3, horiz=FALSE, bty='o', box.col='gray70', box.lwd=2)

pal <- colorRampPalette(c('sienna2', 'dodgerblue4'))
brks = as.numeric(cut(land$Ice06_11, breaks=10))
cols = pal(length(unique(brks)))
cols = as.character(factor(brks, labels=cols))

# pdf('02b_discharge_loadings_reg.pdf', width=7, height=7)
landvar = list(land$WtDepWs, land$Ice06_11, NULL, land$BFIWs, land$WsSlope)
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
par(mar=c(3.5,3.5,0,0))
# defpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,1,1))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
for(trnd in c(1,2,4,5)){
    plot(landvar[[trnd]], dfa$Estimates$Z[,trnd], type='n', yaxt='n',
         xlab='', main='', ylab='', cex.lab=1, font=2, xaxt='n',bty='l')
    mod <- lm(dfa$Estimates$Z[,trnd] ~ landvar[[trnd]])
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    sig <- NULL #turning off sig plotting
    if(trnd == 1){
        mtext(bquote(textstyle(bold('Mean water table depth '))*textstyle(plain('(cm)'))), 1, cex=.8, font=2, line=2)
        print.letter(paste('d',sig),c(.9,.9), cex=1.8, font=2, col='springgreen4')
        mtext(bquote(bold('Shared trend loadings')~
                         textstyle(plain('(')*plain(Delta)*plain(T[water])~plain(Delta)*plain('?'^-1)*plain(')'))),
              side=2, cex=.8, line=1.8, adj=4)
    }
    if(trnd == 2){
        mtext(bquote(bold('Perennial ice/snow cover ')*textstyle(plain('(%)'))), 1, cex=.8, font=2, line=2)
        # mtext(bquote(bold('Total runoff ')*scriptstyle(plain('(mm ')*plain(mo^-1)*plain(')'))), 1, cex=.8, font=2, line=2)
        print.letter(paste('e',sig),c(.9,.9), cex=1.8, font=2, col='springgreen4')
    }
    if(trnd == 4){
        mtext(bquote(bold('BFI ')*textstyle(plain('(%)'))), 1, cex=.8, font=2, line=2)
        print.letter(paste('f',sig),c(.9,.9), cex=1.8, font=2, col='springgreen4')
    }
    if(trnd == 5){
        mtext(bquote(bold('Mean slope ')*textstyle(plain('(% rise)'))), 1, cex=.8, font=2, line=2)
        print.letter(paste('g',sig),c(.1,.9), cex=1.8, font=2, col='springgreen4')
    }
    # if(trnd == 3){
    #     mtext(bquote(bold('Mean slope (% rise)')~.(sig)), 1, cex=1, font=2, line=2.7)
    #     # mtext(bquote(plain(theta) == bold('Precip') ~.(sig)), 3, cex=1.1, font=2, line=1.2)
    # }
    abline(mod, col='springgreen4', lty=2, lwd=3)
    points(landvar[[trnd]], dfa$Estimates$Z[,trnd],
           bg=cols, col='black', cex=rescale(log(land$WsAreaSqKm),c(1.2,3.2)),
           # pch=21,
           pch=land$siteCode,
           lwd=1)
    # color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
    #              rect.col=colorRampPalette(c('brown', 'white'))(10),
    #              align='r', gradient='y')
    # text(3.39, 0.56, labels='Mean watershed')
    # text(3.52, 0.54, labels='elevation (m)')
    axis(1, padj=-.9, tck=-.02); axis(2, padj=.9, tck=-.02)
}

par(defpar)
# dev.off()
# shell('C:\\Users\\Mike\\git\\stream_nuts_DFA\\manuscript\\figures\\16_temp_all_reg.pdf')

#for results (loadings)
summary(mod)

# 1.3 - TEMP eff size and loadings with PCs ####

#load the full atpcsn model

land = readRDS('../../saved_structures/2trendNoSeasNoSnow_land.rds')
# dam_col = rep('black',nrow(land))
# dam_col[land$dam_upstream != 0] <- 'chocolate2'
dam_pch = rep(FALSE,nrow(land))
dam_pch[land$dam_upstream != 0] <- TRUE

dfa = readRDS('../../saved_structures/full_dfaOut.rds')
land = readRDS('../../saved_structures/full_land.rds')

# system('taskkill /f /im AcroRd32.exe')
# pdf('19_temp_simp_reg.pdf', width=7.5, height=5)
pdf('20_temp_simp_reg_NAMES.pdf', width=7.5, height=5)
# layout(matrix(c(1:6,9,7,8),nrow=3,byrow=TRUE))
defpar = par(mfrow=c(2,3), oma=c(0,0,1,1), mar=c(4,3.5,0,0))
landvar=land$ElevWs
# landvar=land$ElevWs/100
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis
#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
#% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# defpar <- par(oma=c(2,4,2,1), mar=c(4,2,2,1))
# pal <- colorRampPalette(c('sienna2', 'dodgerblue4'))
# brks = as.numeric(cut(land$Ice06_11, breaks=10))
# cols = pal(length(unique(brks)))
# cols = as.character(factor(brks, labels=cols))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
for(covr in 1:3){
    res = rescaled_effect_size[,covr]
    # if(covr==3){res = res * 2.54}
    plot(landvar, res, type='n', yaxt='n',
         # xlab='Watershed area over 1000m (%)', main='',
         xlab='', main='', ylab='', cex.lab=1.3, cex.axis=1, font=2, xaxt='n',bty='l')
    mod <- lm(res ~ landvar)
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    sig = NULL #turning off sig plotting
    if(covr == 1){
        mtext(bquote(textstyle(bold(Delta)*bold(T[water])~bold(Delta)*bold(T[air]^-1))~
                         scriptstyle(plain('(')*plain(degree)*plain('C')~plain(degree)*plain(C^-1)*plain(')'))),
              2, cex=.8, line=1.5)
        abline(mod, col='steelblue', lty=2, lwd=3)
        axis(2, padj=.9, tck=-.02)
    }
    if(covr == 2){
        mtext(bquote(textstyle(bold(Delta)*bold(T[water])~bold(Delta)*bold(precip^-1))~
                         textstyle(plain('(')*plain(degree)*plain('C')~plain(cm^-1)*plain(')'))),
              2, cex=.8, line=1.5)
        axis(2, at=c(-.02,0,.02,.04), labels=c(-0.02,0.00,0.02,0.04), padj=.9, tck=-.02)
        # mtext(bquote(textstyle(bold('Perennial watershed ice/snow coverage'))~textstyle(plain('(%)'))), side=1, cex=.8, line=1.7)
        mtext(bquote(textstyle(bold('Mean watershed elevation'))~textstyle(plain('(100 m)'))),
              side=1, cex=.8, line=2)
        abline(mod, col='steelblue', lty=2, lwd=3)
    }
    if(covr == 3){
        mtext(bquote(textstyle(bold(Delta)*bold(T[water])~bold(Delta)*bold(snowmelt^-1))~
                         textstyle(plain('(')*plain(degree)*plain('C')~plain(cm^-1)*plain(')'))), 2, cex=.8, line=1.5)
        # abline(mod, col='steelblue', lty=2, lwd=3)
        axis(2, padj=.9, tck=-.02)
    }
    print.letter(paste(letters[covr],sig), c(.9,.9), cex=1.8, font=2, col='steelblue')
    points(landvar, res,
           # bg=cols, col='black', cex=2,
           bg=cols, col='black', cex=rescale(log(land$WsAreaSqKm),c(1.2,3.2)),
           # pch=21,
           pch=land$siteCode,
           lwd=1)
    points(landvar[dam_pch], res[dam_pch], col='chocolate2', cex=2, lwd=1, pch=124)
    # color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
    #              rect.col=colorRampPalette(c('brown', 'white'))(10),
    #              align='r', gradient='y')
    # text(3.39, 0.56, labels='Mean watershed')
    # text(3.52, 0.54, labels='elevation (m)')
    axis(1, padj=-.9, tck=-.02)
    # if(covr %in% c(1,3)){axis(2)}else{axis(2, )}
}
lines(x=c(-19.9,13.8), y=c(-.342,-.342), xpd=NA, lwd=2, col='gray70')
# mtext(expression(paste(Delta,'Q ', Delta, theta^-1)), side=2, outer=TRUE, cex=1.3, line=1)
# par(defpar)
# dev.off()

#for results (gotta create res above):
# summary(mod)
# mean(res[elev_lo]); sd(res[elev_lo])
# mean(res[elev_med]); sd(res[elev_med])
# mean(res[elev_hi]); sd(res[elev_hi])

#TEMP loading regression

# png('02_loadings_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# pdf('02_loadings_reg.pdf', width=7, height=6)
# defpar <- par(mar=c(5,5,2,5))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
# plot(land$Ice06_11, dfa$Estimates$Z[,1],
#      xlab='Watershed % ice cover', ylab='Trend 1 loadings', type='n',
#      main='', cex.lab=1.3, cex.axis=1, font=2)
# abline(mod, col='gray', lty=2, lwd=3)
# points(land$Ice06_11, dfa$Estimates$Z[,1], col='black', bg=cols,
#        pch=21, cex=1.5, lwd=2)
# mod <- lm(dfa$Estimates$Z[,1] ~ land$Ice06_11)
# color.legend(xl=4,xr=4.4,yb=1.32, yt=1.82, legend=c('147', '1349'),
#              rect.col=colorRampPalette(c('brown', 'white'))(10),
#              align='r', gradient='y')
# text(3.39, 1.61, labels='Mean watershed')
# text(3.52, 1.53, labels='elevation (m)')
# par(defpar)
# dev.off()
plot(1,1,type='n',ann=FALSE,axes=FALSE)
legend(x=.95, y=1.2, legend=c('Rain-dominated','Rain-and-snow','Snow-dominated','Dam-influenced'),
       xpd=NA, pt.bg=c('black','gray75','white','white'), pch=c(21,21,21,124), xjust=0.5,
       col=c('black','black','black','chocolate2'), cex=1.3, horiz=FALSE, bty='n')#, box.col='gray70', box.lwd=2)

# pal <- colorRampPalette(c('sienna2', 'dodgerblue4'))
# brks = as.numeric(cut(land$Ice06_11, breaks=10))
# cols = pal(length(unique(brks)))
# cols = as.character(factor(brks, labels=cols))

# pdf('02b_discharge_loadings_reg.pdf', width=7, height=7)
land2 = land
dfa = readRDS('../../saved_structures/2trendNoSeasNoSnow_dfaOut.rds')
land = readRDS('../../saved_structures/2trendNoSeasNoSnow_land.rds')

landvar = list(land$PCo2, land$PCo1)
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# par(mar=c(3.5,3.5,0,0))
# defpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,1,1))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]


# elev_hi = which(land$Ice06_11 >= 0.7)
# elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
# elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
# cols = rep('black', nrow(land))
# cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
for(trnd in c(2,1)){
    plot(landvar[[trnd]], dfa$Estimates$Z[,trnd], type='n', yaxt='n',
         xlab='', main='', ylab='', cex.lab=1, font=2, xaxt='n',bty='l')
    mod <- lm(dfa$Estimates$Z[,trnd] ~ landvar[[trnd]])
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    sig <- NULL #turning off sig plotting
    if(trnd == 2){
        mtext(bquote(textstyle(bold('Principal axis 1 '))*textstyle(plain('(50.8%)'))), 1, cex=.8, font=2, line=2)
        print.letter(paste('d',sig),c(.1,.9), cex=1.8, font=2, col='springgreen4')
        mtext(bquote(bold('Shared trend loadings')~
                         textstyle(plain('(')*plain(Delta)*plain(T[water])~plain(Delta)*plain('?'^-1)*plain(')'))),
              side=2, cex=.8, line=1.8)#, adj=4)
    }
    if(trnd == 1){
        mtext(bquote(bold('Principal axis 2 ')*textstyle(plain('(28.6%)'))), 1, cex=.8, font=2, line=2)
        # mtext(bquote(bold('Total runoff ')*scriptstyle(plain('(mm ')*plain(mo^-1)*plain(')'))), 1, cex=.8, font=2, line=2)
        print.letter(paste('e',sig),c(.9,.9), cex=1.8, font=2, col='springgreen4')
    }
    abline(mod, col='springgreen4', lty=2, lwd=3)
    points(landvar[[trnd]], dfa$Estimates$Z[,trnd],
           # bg=cols, col='black', cex=2,
           bg=cols, col='black', cex=rescale(log(land2$WsAreaSqKm),c(1.2,3.2)),
           pch=21,
           # pch=land$siteCode,
           lwd=1)
    points(landvar[[trnd]][dam_pch], dfa$Estimates$Z[,trnd][dam_pch], col='chocolate2',
           cex=2, lwd=1, pch=124)
    # color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
    #              rect.col=colorRampPalette(c('brown', 'white'))(10),
    #              align='r', gradient='y')
    # text(3.39, 0.56, labels='Mean watershed')
    # text(3.52, 0.54, labels='elevation (m)')
    axis(1, padj=-.9, tck=-.02); axis(2, padj=.9, tck=-.02)
}

par(defpar)
dev.off()
shell('C:\\Users\\Mike\\git\\stream_nuts_DFA\\manuscript\\figures\\20_temp_simp_reg_NAMES.pdf')
# shell('C:\\Users\\Mike\\git\\stream_nuts_DFA\\manuscript\\figures\\19_temp_simp_reg.pdf')

#for results (loadings)
# summary(mod)

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

# 9.1 - DISCHARGE effect size regression (linked with loading regression) ####

# system('taskkill /f /im AcroRd32.exe')
pdf('14_discharge_all_reg.pdf', width=7.5, height=7.5)
layout(matrix(c(1:6,9,7,8),nrow=3,byrow=TRUE))
defpar = par(oma=c(0,0,1,1), mar=c(4,3.5,0,0))
landvar=land$ElevWs/100
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
#% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# defpar <- par(oma=c(2,4,2,1), mar=c(4,2,2,1))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
for(covr in 1:3){
res = rescaled_effect_size[,covr]
if(covr==3){res = res * 2.54}
plot(landvar, res, type='n', yaxt='n',
     # xlab='Watershed area over 1000m (%)', main='',
     xlab='', main='', ylab='', cex.lab=1.3, cex.axis=1, font=2, xaxt='n',bty='l')
mod <- lm(res ~ I(landvar))
sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
if(covr == 1){
    mtext(bquote(textstyle(bold(Delta)*bold('Q ')*bold(Delta)*bold(T[air]^-1))~
              scriptstyle(plain('(CFS')~plain(degree)*plain(C^-1)*plain(')'))),
          2, cex=.8, line=1.5)
    axis(2, at=c(.07,.09,.11), padj=.9, tck=-.02)
}
if(covr == 2){
    mtext(bquote(textstyle(bold(Delta)*bold('Q ')*bold(Delta)*bold(precip^-1))~
                     scriptstyle(plain('(CFS')~plain(cm^-1)*plain(')'))),
          2, cex=.8, line=1.5)
    abline(mod, col='steelblue', lty=2, lwd=3)
    axis(2, at=c(-.005,0,.005,.01), labels=c(expression(plain(-5)*plain(e^-3)),
        0,expression(plain(5)*plain(e^-3)),
        expression(plain(1)*plain(e^-2))), padj=.8, tck=-.02)
    mtext(bquote(textstyle(bold('Mean watershed elevation'))~scriptstyle(plain('(100 m)'))), side=1, cex=.8, line=1.7)
}
if(covr == 3){
    mtext(bquote(textstyle(bold(Delta)*bold('Q ')*bold(Delta)*bold(snowmelt^-1))~
              scriptstyle(plain('(CFS')~plain(cm^-1)*plain(')'))), 2, cex=.8, line=1.5)
    abline(mod, col='steelblue', lty=2, lwd=3)
    axis(2, at=c(.1,.3,.5,.7), padj=.9, tck=-.02)
}
print.letter(paste(letters[covr],sig), c(.1,.9), cex=1.2, font=2, col='steelblue')
points(landvar, res,
       bg=cols, col='black',
       pch=21, cex=1.5, lwd=2)
# color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
#              rect.col=colorRampPalette(c('brown', 'white'))(10),
#              align='r', gradient='y')
# text(3.39, 0.56, labels='Mean watershed')
# text(3.52, 0.54, labels='elevation (m)')
axis(1, padj=-.9, tck=-.02)
# if(covr %in% c(1,3)){axis(2)}else{axis(2, )}
}
lines(x=c(-27,13.8), y=c(-.115,-.115), xpd=NA, lwd=2, col='gray70')
# mtext(expression(paste(Delta,'Q ', Delta, theta^-1)), side=2, outer=TRUE, cex=1.3, line=1)
# par(defpar)
# dev.off()

# 9.2 - DISCHARGE loading regression ####

# png('02_loadings_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# pdf('02_loadings_reg.pdf', width=7, height=6)
# defpar <- par(mar=c(5,5,2,5))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
# plot(land$Ice06_11, dfa$Estimates$Z[,1],
#      xlab='Watershed % ice cover', ylab='Trend 1 loadings', type='n',
#      main='', cex.lab=1.3, cex.axis=1, font=2)
# abline(mod, col='gray', lty=2, lwd=3)
# points(land$Ice06_11, dfa$Estimates$Z[,1], col='black', bg=cols,
#        pch=21, cex=1.5, lwd=2)
# mod <- lm(dfa$Estimates$Z[,1] ~ land$Ice06_11)
# color.legend(xl=4,xr=4.4,yb=1.32, yt=1.82, legend=c('147', '1349'),
#              rect.col=colorRampPalette(c('brown', 'white'))(10),
#              align='r', gradient='y')
# text(3.39, 1.61, labels='Mean watershed')
# text(3.52, 1.53, labels='elevation (m)')
# par(defpar)
# dev.off()
plot(1,1,type='n',ann=FALSE,axes=FALSE)
legend(x=1,y=1.505, legend=c('Rain-dominated','Rain-and-snow','Snow-dominated'),
       xpd=NA, pt.bg=c('black','gray75','white'), pch=21, xjust=0.5,
       col='black', cex=1.3, horiz=FALSE, bty='o', box.col='gray70', box.lwd=2)

# pdf('02b_discharge_loadings_reg.pdf', width=7, height=7)
landvar = list(NULL,land$WsSlope, land$RunoffWs, land$WsAreaSqKm, land$RckDepWs)
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
par(mar=c(3.5,3.5,0,0))
# defpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,1,1))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
for(trnd in 2:5){
    plot(landvar[[trnd]], dfa$Estimates$Z[,trnd], type='n', yaxt='n',
         xlab='', main='', ylab='', cex.lab=1, font=2, xaxt='n',bty='l')
    mod <- lm(dfa$Estimates$Z[,trnd] ~ landvar[[trnd]])
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    if(trnd == 2){
        mtext(bquote(textstyle(bold('Mean slope '))*scriptstyle(plain('(% rise)'))), 1, cex=.8, font=2, line=2)
        print.letter(paste(letters[trnd+2],sig),c(.9,.9), cex=1.2, font=2, col='springgreen4')
        mtext(bquote(bold('Shared trend loadings')~
                         scriptstyle(plain('(')*plain(Delta)*plain('Q ')*plain(Delta)*plain('?'^-1)*plain(')'))),
                     side=2, cex=.8, line=1.8, adj=-14)
    }
    if(trnd == 3){
        mtext(bquote(bold('Total runoff ')*scriptstyle(plain('(mm ')*plain(mo^-1)*plain(')'))), 1, cex=.8, font=2, line=2)
        print.letter(paste(letters[trnd+2],sig),c(.1,.9), cex=1.2, font=2, col='springgreen4')
    }
    if(trnd == 4){
        mtext(bquote(bold('Total area ')*scriptstyle(plain('(')*plain(km^2)*plain(')'))), 1, cex=.8, font=2, line=2)
        print.letter(paste(letters[trnd+2],sig),c(.9,.9), cex=1.2, font=2, col='springgreen4')
    }
    if(trnd == 5){
        mtext(bquote(bold('Mean bedrock depth ')*scriptstyle(plain('(cm)'))), 1, cex=.8, font=2, line=2)
        print.letter(paste(letters[trnd+2],sig),c(.1,.9), cex=1.2, font=2, col='springgreen4')
    }
    # if(trnd == 3){
    #     mtext(bquote(bold('Mean slope (% rise)')~.(sig)), 1, cex=1, font=2, line=2.7)
    #     # mtext(bquote(plain(theta) == bold('Precip') ~.(sig)), 3, cex=1.1, font=2, line=1.2)
    # }
    abline(mod, col='springgreen4', lty=2, lwd=3)
    points(landvar[[trnd]], dfa$Estimates$Z[,trnd],
           bg=cols, col='black',
           pch=21, cex=1.5, lwd=2)
    # color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
    #              rect.col=colorRampPalette(c('brown', 'white'))(10),
    #              align='r', gradient='y')
    # text(3.39, 0.56, labels='Mean watershed')
    # text(3.52, 0.54, labels='elevation (m)')
    axis(1, padj=-.9, tck=-.02); axis(2, padj=.9, tck=-.02)
}

par(defpar)
dev.off()
shell('C:\\Users\\Mike\\git\\stream_nuts_DFA\\manuscript\\figures\\14_discharge_all_reg.pdf')

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
# saveRDS(list(covs, obs_ts), '../../saved_structures/air_water_discharge.rds') # from temp_due_4m_at_byMo_allMos.rda
# saveRDS(list(covs[,1,drop=FALSE], obs_ts), '../../saved_structures/air_water_discharge2.rds') # from temp_due_5m_atpcsn_byMo_allMos.rda
# saveRDS(obs_ts, '../../saved_structures/just_discharge.rds') #from discharge, for awd plot

# 16 - DISCHARGE vs. TEMP ####

#load discharge_due_5m_atpcsn_byMo_allMos.rda in the setup section.

disch_moInts <- readRDS('../../saved_structures/moInts_discharge_due_4m_at.rds')
temp_moInts <- readRDS('../../saved_structures/moInts_temp_due_4m_at.rds')
# "A"  "B"  "C"  "E"  "F"  "G"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "S"  "T"  "U"
# [20] "V"  "W"  "X"  "Z"  "ZA"
# "A"  "B"  "C"  "E"  "F"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "U"  "W"  "Z"  "ZA"

pdf('11a_disch_vs_air_byWsAreaOver1000.pdf', width=8, height=8)
# pdf('11b_disch_vs_air_byWsArea.pdf', width=8, height=8)
# pal <- colorRampPalette(c('white', 'black'))
pal <- colorRampPalette(c('goldenrod', 'blue'))
# pal <- colorRampPalette(c('black', 'white'))
# brks = as.numeric(cut(land$WsAreaSqKm, breaks=10))
brks = as.numeric(cut(land$WsAreaOver1000*100, breaks=10))
# brks = as.numeric(cut(land$ElevWs, breaks=10))
# brks = as.numeric(cut(land$Ice06_11, breaks=10))
cols = pal(length(unique(brks)))
cols = as.character(factor(brks, labels=cols))
defpar <- par(mfcol=c(6,2), oma=c(5,8,1,3), mar=c(0.5,0.5,0.5,0.5))
for(i in 1:12){
    mod <- lm(temp_moInts[,i][-c(6,17,18,20,22)] ~ disch_moInts[,i])
    plot(disch_moInts[,i], temp_moInts[,i][-c(6,17,18,20,22)], cex=1.5,
         yaxt='n', xaxt='n', bty='l',
         # pch=21, bg=cols, col='black',
         pch=land$siteCode, col=cols,
         xlim=c(min(disch_moInts), max(disch_moInts)),
         ylim=c(min(temp_moInts[-c(6,17,18,20,22)]), max(temp_moInts[-c(6,17,18,20,22)])))
    abline(mod, col='darkred', lty=2, lwd=2)
    abline(h=0, lty=3, col='steelblue4')
    abline(v=0, lty=3, col='steelblue4')
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

# 17 - DISCHARGE vs. TEMP (one column) ####

#load discharge_due_4m_atpc_byMo_allMos.rda in the setup section, or else names
#will be screwed up.

disch_moInts <- readRDS('../../saved_structures/moInts_discharge_due_4m_at.rds')
temp_moInts <- readRDS('../../saved_structures/moInts_temp_due_4m_at.rds')
# "A"  "B"  "C"  "E"  "F"  "G"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "S"  "T"  "U"
# [20] "V"  "W"  "X"  "Z"  "ZA"
# "A"  "B"  "C"  "E"  "F"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "U"  "W"  "Z"  "ZA"

pdf('13_disch_vs_air_column.pdf', width=6, height=12)
# pdf('11b_disch_vs_air_byWsArea.pdf', width=8, height=8)
# pal <- colorRampPalette(c('white', 'black'))
pal <- colorRampPalette(c('black', 'white'))
# pal <- colorRampPalette(c('black', 'white'))
# brks = as.numeric(cut(land$WsAreaSqKm, breaks=10))
# brks = as.numeric(cut(land$ElevWs, breaks=10))
brks = as.numeric(cut(land$Ice06_11, breaks=10))
# brks = as.numeric(cut(land$WsAreaOver1000*100, breaks=10))
# layout(matrix(c(1,1,2,3,4,5), ncol=2, byrow=T), heights=c(.3,1,1,1,1))
# defpar <- par(oma=c(5,3.5,0,4), mar=c(0,0,0,0))
defpar <- par(mfcol=c(13,1), oma=c(5,3.5,0,4), mar=c(0,0,0,0))
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray75'; cols[elev_hi] = 'white'
# cols = pal(length(unique(brks)))
# cols = as.character(factor(brks, labels=cols))
# for(i in c(0,2,5,8,11)){
for(i in 0:12){
    if(i == 0) {
        plot(1,1,axes=FALSE,,xlim=c(0,1),ylim=c(0,1),ann=FALSE,type='n')
        color.legend(xl=.4,xr=.6,yb=.5, yt=.7, legend=c('Rainfed', 'Snowfed'),
                 rect.col=colorRampPalette(c('black', 'white'))(6),
                 gradient='x', align='lt')
    }
    if(i != 0){
        mod <- lm(temp_moInts[,i][-c(6,17,18,20,22)] ~ disch_moInts[,i])
        sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
        plot(disch_moInts[,i], temp_moInts[,i][-c(6,17,18,20,22)],
             yaxt='n', xaxt='n', bty='l', type='n', axes=FALSE,
             xlim=c(min(disch_moInts), max(disch_moInts)),
             ylim=c(min(temp_moInts[-c(6,17,18,20,22)]), max(temp_moInts[-c(6,17,18,20,22)])))
        polygon(x=c(0,0.4,0.4,0), y=c(0,0,1.3,1.3),
                col=adjustcolor('darkred', alpha.f=0.15), border=FALSE)
        polygon(x=c(0,-0.5,-0.5,0), y=c(0,0,-0.65,-0.65),
                col=adjustcolor('navy', alpha.f=0.15), border=FALSE)
        abline(mod, col='gray30', lty=1, lwd=2)
        # abline(h=0, lty=3, col='gray50')
        points(mean(disch_moInts[elev_lo,i]),
                         mean(temp_moInts[-c(6,17,18,20,22),i][elev_lo]),
                         pch=24, cex=3, lwd=1, col='black', bg='black')
        points(mean(disch_moInts[elev_med,i]),
                         mean(temp_moInts[-c(6,17,18,20,22),i][elev_lo]),
                         pch=24, cex=3, lwd=1, col='black', bg='gray75')
        points(mean(disch_moInts[elev_hi,i]),
                         mean(temp_moInts[-c(6,17,18,20,22),i][elev_hi]),
                         # pch=4, cex=4, lwd=3, col='steelblue3')
                         pch=24, cex=3, lwd=1, col='black', bg='white')
        points(disch_moInts[,i], temp_moInts[,i][-c(6,17,18,20,22)], cex=1.5,
               # pch=21, bg=cols, col='black')
               pch=land$siteCode, col='black')
        mtext(paste0(month.abb[i], sig), 4, line=.5, las=2)
        if(i==1){
            # axis(4, las=2, line=.6, at=c(-.5,0,.5,1),
            #           labels=c('-0.5', sprintf('%4s', c('0.0','0.5','1.0'))))
            axis(2, las=2, line=0)
            lines(x=c(-.247,-.247), y=c(-.65,1.3))
        }
        # if(i %in% 7:12) mtext(month.abb[i], 4, line=1, las=2)
        if(i == 12) axis(1)
    }
}
# lines(x=c(0,0), y=c(-0.8, 23.04), xpd=NA, col='gray50', lty=3, lwd=1)
lines(x=c(-.247,.308), y=c(-.62,-.62), xpd=NA)
mtext(expression(paste(~italic('ln')~bold(Delta)~bold('Q (CFS) / ')
                       ~bold(Delta)~bold('air')~bold(degree)~bold('C'), sep='')),
      # ~bold(Delta)~bold('precip. (cm)'), sep='')),
      side=1, outer=TRUE, line=3)
mtext(expression(paste(bold(Delta)~bold('water')~bold(degree)~bold('C / ')
                       ~bold(Delta)~bold('air')~bold(degree)~bold('C'), sep='')),
                       # ~bold(Delta)~bold('air')~bold(degree)~bold('C')^-1, sep='')),
      side=2, outer=TRUE, line=1)
legend(x=-.02, y=24.1, legend='', pch=4, pt.cex=4, pt.lwd=3,
       col=c('orange2','steelblue3'), xpd=NA, horiz=TRUE, bty='n')
legend(x=.06, y=24.1, legend='', pch=4, pt.cex=4, pt.lwd=3,
       col=c('steelblue3'), xpd=NA, horiz=TRUE, bty='n')
legend(x=-.21, y=24.8, legend='', xpd=NA,
       fill=adjustcolor('navy', alpha.f=0.2), bty='n', border=NA, cex=2)
text(-.135,24.05,'+air -Q -water',xpd=NA)
legend(x=.21, y=24.8, legend='', xpd=NA,
       fill=adjustcolor('darkred', alpha.f=0.2), bty='n', border=NA, cex=2)
text(.189,24.05,'+air +Q +water',xpd=NA)
text(0.03,23.65,'(Centroids)',xpd=NA)
# legend.gradient(pts=cbind(x=c(-.04,.04,.04,-.04), y=c(22,22,24,24)), cols=c('black','white'),
#                 limits=c('Rainfed', 'Snowfed'), title='X = mean')
par(defpar)
dev.off()

# 18 - DISCHARGE vs. TEMP (big four) ####

#load discharge_due_4m_atpc_byMo_allMos.rda in the setup section, or else names
#will be screwed up.

disch_moInts <- readRDS('../../saved_structures/moInts_discharge_due_4m_at.rds')
temp_moInts <- readRDS('../../saved_structures/moInts_temp_due_4m_at.rds')
# "A"  "B"  "C"  "E"  "F"  "G"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "S"  "T"  "U"
# [20] "V"  "W"  "X"  "Z"  "ZA"
# "A"  "B"  "C"  "E"  "F"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "U"  "W"  "Z"  "ZA"

pdf('14_disch_vs_air_focal.pdf', width=6, height=7)
# pdf('11b_disch_vs_air_byWsArea.pdf', width=8, height=8)
# pal <- colorRampPalette(c('white', 'black'))
pal <- colorRampPalette(c('black', 'white'))
# pal <- colorRampPalette(c('black', 'white'))
# brks = as.numeric(cut(land$WsAreaSqKm, breaks=10))
# brks = as.numeric(cut(land$WsAreaOver1000*100, breaks=10))
# brks = as.numeric(cut(land$ElevWs, breaks=10))
brks = as.numeric(cut(land$Ice06_11, breaks=10))
# cols = pal(length(unique(brks)))
# cols = as.character(factor(brks, labels=cols))
# layout(matrix(c(1,2,3,4), ncol=2, byrow=T), heights=c(1,1,1,1))
# layout(matrix(c(1,1,2,3,4,5), ncol=2, byrow=T), heights=c(.2,1,1,1,1))
       # respect=c(matrix(1,ncol=2,nrow=3))
defpar <- par(mfrow=c(2,2), oma=c(5,4,5,2), mar=c(0,1,0,0))
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray75'; cols[elev_hi] = 'white'
for(i in c(2,5,8,11)){
    # if(i == 0) {
        # plot(1,1,axes=FALSE,,xlim=c(0,1),ylim=c(0,1),ann=FALSE,type='n')
        # color.legend(xl=.4,xr=.6,yb=.5, yt=.7, legend=c('Rainfed', 'Snowfed'),
        #          rect.col=colorRampPalette(c('black', 'white'))(6),
        #          gradient='x', align='lt')
    # }
    if(i != 0){
        mod <- lm(temp_moInts[,i][-c(6,17,18,20,22)] ~ disch_moInts[,i])
        sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
        plot(disch_moInts[,i], temp_moInts[,i][-c(6,17,18,20,22)],
             yaxt='n', xaxt='n', bty='l', type='n', axes=FALSE,
             xlim=c(-.15, max(disch_moInts)),
             ylim=c(min(temp_moInts[-c(6,17,18,20,22)]), max(temp_moInts[-c(6,17,18,20,22)])))
        polygon(x=c(0,0.4,0.4,0), y=c(0,0,1.3,1.3),
                col=adjustcolor('darkred', alpha.f=0.2), border=FALSE)
        polygon(x=c(0,-0.5,-0.5,0), y=c(0,0,-0.65,-0.65),
                col=adjustcolor('navy', alpha.f=0.15), border=FALSE)
        if(i == 5) abline(mod, col='gray30', lty=1, lwd=2)
        # abline(h=0, lty=3, col='gray50')
        # points(mean(disch_moInts[elev_lo,i]),
        #                  mean(temp_moInts[-c(6,17,18,20,22),i][elev_lo]),
        #                  pch=24, cex=3, lwd=1, col='black', bg='black')
        # points(mean(disch_moInts[elev_med,i]),
        #                  mean(temp_moInts[-c(6,17,18,20,22),i][elev_lo]),
        #                  pch=24, cex=3, lwd=1, col='black', bg='gray75')
        # points(mean(disch_moInts[elev_hi,i]),
        #                  mean(temp_moInts[-c(6,17,18,20,22),i][elev_hi]),
        #                  # pch=4, cex=4, lwd=3, col='steelblue3')
        #                  pch=24, cex=3, lwd=1, col='black', bg='white')
        points(disch_moInts[,i], temp_moInts[,i][-c(6,17,18,20,22)], cex=1.5,
               pch=21, bg=cols, col='black')
               # pch=land$siteCode, col=cols,
        text(-.1,1.1,paste0(month.abb[i], sig), cex=1.2)
        if(i %in% c(2,8)){
            # axis(4, las=2, line=.6, at=c(-.5,0,.5,1),
            #           labels=c('-0.5', sprintf('%4s', c('0.0','0.5','1.0'))))
            axis(2, las=2, line=0)
            lines(x=c(-.167,-.167), y=c(-.62,1.3))
        }
        # if(i %in% 7:12) mtext(month.abb[i], 4, line=1, las=2)
        if(i %in% c(8,11)){
            axis(1)
            lines(x=c(-.165,.308), y=c(-.62,-.62), xpd=NA)
        }
    }
}
# lines(x=c(0,0), y=c(-0.8, 23.04), xpd=NA, col='gray50', lty=3, lwd=1)
mtext(expression(paste(~bolditalic('ln')~bold(Delta)~bold('Q (CFS)'),sep='')),# / ')
                       # ~bold(Delta)~bold('air')~bold(degree)~bold('C'), sep='')),
      # ~bold(Delta)~bold('precip. (cm)'), sep='')),
      side=1, outer=TRUE, line=3)
mtext(expression(paste(bold(Delta)~bold('water')~bold(degree)~bold('C'))),# / ')
                       # ~bold(Delta)~bold('air')~bold(degree)~bold('C'), sep='')),
                       # ~bold(Delta)~bold('air')~bold(degree)~bold('C')^-1, sep='')),
      side=2, outer=TRUE, line=2)
# legend(x=-.02, y=24.1, legend='', pch=4, pt.cex=4, pt.lwd=3,
#        col=c('orange2','steelblue3'), xpd=NA, horiz=TRUE, bty='n')
# legend(x=.06, y=24.1, legend='', pch=4, pt.cex=4, pt.lwd=3,
#        col=c('steelblue3'), xpd=NA, horiz=TRUE, bty='n')
# legend(x=-.5, y=24.8, legend='', xpd=NA,
#        fill=adjustcolor('navy', alpha.f=0.2), bty='n', border=NA, cex=2)
# text(-.135,24.05,'+air -Q -water',xpd=NA)
# legend(x=.21, y=24.8, legend='', xpd=NA,
#        fill=adjustcolor('darkred', alpha.f=0.2), bty='n', border=NA, cex=2)
# text(.189,24.05,'+air +Q +water',xpd=NA)
# legend(x=-.68, y=3.7, legend=c('Rain-fed'), xpd=NA, pt.bg=c('black'), pch=21,
#        col='black', cex=1.5, horiz=TRUE, bty='n', adj=c(.1,.4))
# legend(x=-.52, y=3.7, legend=c('Rain-and-snow'), xpd=NA, pt.bg=c('gray75'), pch=21,
#        col='black', cex=1.5, horiz=TRUE, bty='n', adj=c(.05,.4))
# legend(x=-.28, y=3.7, legend=c('Snow-fed'), xpd=NA, pt.bg=c('white'), pch=21,
#        col='black', cex=1.5, horiz=TRUE, bty='n', adj=c(.1,.4))
legend(x=-.68, y=3.85, legend=c('Rain-dominated','Rain-and-snow','Snow-dominated'),
       xpd=NA, pt.bg=c('black','gray75','white'), pch=21,
       col='black', cex=1.3, horiz=FALSE, bty='n')
legend(x=0.13, y=4, legend=c(expression(paste(plain('Q, T'[water]))),
                            expression(paste(frac(1,Q), ' ,  ', frac(1,T[water])))),
       xpd=NA, fill=c(adjustcolor('darkred', alpha.f=0.2),
                      adjustcolor('navy', alpha.f=0.15)),
       border='black', cex=1.1, horiz=FALSE, bty='n', adj=c(0,.4))
# legend(x=0.08, y=3.75, legend=c(expression(paste(plain('' %up% 'Q, '%up% 'T'[water]))),
#                             expression(paste(plain('' %down% 'Q, '%down% 'T'[water])))),
#        xpd=NA, fill=c(adjustcolor('darkred', alpha.f=0.2),
#                       adjustcolor('navy', alpha.f=0.15)),
#        border='black', cex=1.2, horiz=FALSE, bty='n', adj=c(0,.4))
text(0.08, 3.58, expression(paste(plain('T'[air]) %prop% '')),
     xpd=NA, cex=1.6)
# text(0.01, 3.58, expression(paste(plain('' %up% 'T'[air]) %=>% '')),
#      xpd=NA, cex=1.6)
# text(0.03,23.65,'(Centroids)',xpd=NA)
# legend.gradient(pts=cbind(x=c(-.04,.04,.04,-.04), y=c(22,22,24,24)), cols=c('black','white'),
#                 limits=c('Rainfed', 'Snowfed'), title='X = mean')
par(defpar)
dev.off()

# 19 - DISCHARGE vs. TEMP (big five) ####

#load discharge_due_5m_atpcsn_byMo_allMos.rda in the setup section, or else names
#will be screwed up.
dams = read.csv('../../data/watershed_data/watershed_data_simp.csv', stringsAsFactors=FALSE)[,c('siteCode','dam_upstream')]
dams$siteCode[dams$siteCode=='AA'] <- 'ZA'
land = merge(land, dams, by='siteCode')


library(plotrix)
#these are generated in the effect size v month plots
disch_moInts <- readRDS('../../saved_structures/moInts_discharge_due_5m_atpcsn.rds')
temp_moInts <- readRDS('../../saved_structures/moInts_temp_due_5m_atpcsn.rds')
# "A"  "B"  "C"  "E"  "F"  "G"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "S"  "T"  "U"
# [20] "V"  "W"  "X"  "Z"  "ZA"
# "A"  "B"  "C"  "E"  "F"  "H"  "I"  "J"  "L"  "M"  "N"  "O"  "P"  "Q"  "R"  "U"  "W"  "Z"  "ZA"

pdf('15_disch_vs_air_focal5.pdf', width=5, height=6.4)
pdf('21_disch_vs_air_focal5_NAMES.pdf', width=5, height=6.4)
# pal <- colorRampPalette(c('white', 'black'))
pal <- colorRampPalette(c('black', 'white'))
# pal <- colorRampPalette(c('black', 'white'))
# brks = as.numeric(cut(land$WsAreaSqKm, breaks=10))
# brks = as.numeric(cut(land$WsAreaOver1000*100, breaks=10))
# brks = as.numeric(cut(land$ElevWs, breaks=10))
brks = as.numeric(cut(land$Ice06_11, breaks=10))
# cols = pal(length(unique(brks)))
# cols = as.character(factor(brks, labels=cols))
# layout(matrix(c(1,2,3,4), ncol=2, byrow=T), heights=c(1,1,1,1))
# layout(matrix(c(1,1,2,3,4,5), ncol=2, byrow=T), heights=c(.2,1,1,1,1))
# respect=c(matrix(1,ncol=2,nrow=3))
defpar <- par(mfrow=c(3,2), oma=c(5,5,1,1), mar=c(.5,.5,0,0))
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray75'; cols[elev_hi] = 'white'
# dam_col = rep('black',nrow(land))
# dam_col[land$dam_upstream != 0] <- 'chocolate2'
dam_pch = rep(FALSE,nrow(land))
dam_pch[land$dam_upstream != 0] <- TRUE
for(i in c(0,2,4,6,8,11)){
    if(i == 0) {
    plot(1,1,axes=FALSE,,xlim=c(0,1),ylim=c(0,1),ann=FALSE,type='n')
    # color.legend(xl=.4,xr=.6,yb=.5, yt=.7, legend=c('Rainfed', 'Snowfed'),
    #          rect.col=colorRampPalette(c('black', 'white'))(6),
    #          gradient='x', align='lt')
    }
    if(i != 0){
        ymin=-.07
        ymax=max(temp_moInts[-c(6,17,18,20,22)])
        mod <- lm(temp_moInts[,i][-c(6,17,18,20,22)] ~ disch_moInts[,i])
        sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
        sig = NULL #suppressing sig plotting
        plot(disch_moInts[,i], temp_moInts[,i][-c(6,17,18,20,22)],
             yaxt='n', xaxt='n', bty='o', type='n', fg='gray60',
             xlim=c(-.15, max(disch_moInts)),
             ylim=c(ymin, ymax), yaxs='i')
        polygon(x=c(0,0.4,0.4,0), y=c(0,0,ymax,ymax),
                col=adjustcolor('darkred', alpha.f=0.2), border=FALSE)
        polygon(x=c(0,-.5,-.5,0), y=c(0,0,ymin,ymin),
                col=adjustcolor('navy', alpha.f=0.15), border=FALSE)
        if(i %in% c(4,5)) abline(mod, col='gray30', lty=2, lwd=1)
        # abline(h=0, lty=3, col='gray50')
        # points(mean(disch_moInts[elev_lo,i]),
        #                  mean(temp_moInts[-c(6,17,18,20,22),i][elev_lo]),
        #                  pch=24, cex=3, lwd=1, col='black', bg='black')
        # points(mean(disch_moInts[elev_med,i]),
        #                  mean(temp_moInts[-c(6,17,18,20,22),i][elev_lo]),
        #                  pch=24, cex=3, lwd=1, col='black', bg='gray75')
        # points(mean(disch_moInts[elev_hi,i]),
        #                  mean(temp_moInts[-c(6,17,18,20,22),i][elev_hi]),
        #                  # pch=4, cex=4, lwd=3, col='steelblue3')
        #                  pch=24, cex=3, lwd=1, col='black', bg='white')
        points(disch_moInts[,i], temp_moInts[,i][-c(6,17,18,20,22)],# cex=1.5,
               # pch=21, bg=cols, col='black',
               pch=land$siteCode, col='black')
               # cex=1.3)
               # cex=rescale(log(land$elevation_),c(1.2,3.2))) #WtDepWs WsAreaSqKm BFIWs WsSlope elevation_
        points(disch_moInts[,i][dam_pch], temp_moInts[,i][-c(6,17,18,20,22)][dam_pch],
               pch='|', col='chocolate2', cex=1.3, lwd=1)
        print.letter(paste0(month.abb[i], sig), c(.9,.9), cex=1.2, font=2)
        if(i %in% c(2,4,8)){
            # axis(4, las=2, line=.6, at=c(-.5,0,.5,1),
            #           labels=c('-0.5', sprintf('%4s', c('0.0','0.5','1.0'))))
            axis(2, las=2, line=0, col='gray60')
            # lines(x=c(-.167,-.167), y=c(-.62,1.3))
        }
        # if(i %in% 7:12) mtext(month.abb[i], 4, line=1, las=2)
        if(i %in% c(8,11)){
            axis(1, col='gray60')
            # lines(x=c(-.165,.308), y=c(-.62,-.62), xpd=NA)
        }
        # if(i == 8){
        #     axis(1, col='gray60', at=c(-.1,0,.1,.2))
        # }
    }
}
# lines(x=c(0,0), y=c(-0.8, 23.04), xpd=NA, col='gray50', lty=3, lwd=1)
mtext(bquote(paste(bolditalic('ln')*bold(Delta)*bold('Q')~plain('(CFS)'))),# / ')
      # ~bold(Delta)~bold('air')~bold(degree)~bold('C'), sep='')),
      # ~bold(Delta)~bold('precip. (cm)'), sep='')),
      side=1, outer=TRUE, line=3)
mtext(bquote(paste(bold(Delta)*bold(T[water])~plain('(')*plain(degree)*plain('C)'))),# / ')
      # ~bold(Delta)~bold('air')~bold(degree)~bold('C'), sep='')),
      # ~bold(Delta)~bold('air')~bold(degree)~bold('C')^-1, sep='')),
      side=2, outer=TRUE, line=3)
# legend(x=-.02, y=24.1, legend='', pch=4, pt.cex=4, pt.lwd=3,
#        col=c('orange2','steelblue3'), xpd=NA, horiz=TRUE, bty='n')
# legend(x=.06, y=24.1, legend='', pch=4, pt.cex=4, pt.lwd=3,
#        col=c('steelblue3'), xpd=NA, horiz=TRUE, bty='n')
# legend(x=-.5, y=24.8, legend='', xpd=NA,
#        fill=adjustcolor('navy', alpha.f=0.2), bty='n', border=NA, cex=2)
# text(-.135,24.05,'+air -Q -water',xpd=NA)
# legend(x=.21, y=24.8, legend='', xpd=NA,
#        fill=adjustcolor('darkred', alpha.f=0.2), bty='n', border=NA, cex=2)
# text(.189,24.05,'+air +Q +water',xpd=NA)
# legend(x=-.68, y=3.7, legend=c('Rain-fed'), xpd=NA, pt.bg=c('black'), pch=21,
#        col='black', cex=1.5, horiz=TRUE, bty='n', adj=c(.1,.4))
# legend(x=-.52, y=3.7, legend=c('Rain-and-snow'), xpd=NA, pt.bg=c('gray75'), pch=21,
#        col='black', cex=1.5, horiz=TRUE, bty='n', adj=c(.05,.4))
# legend(x=-.28, y=3.7, legend=c('Snow-fed'), xpd=NA, pt.bg=c('white'), pch=21,
#        col='black', cex=1.5, horiz=TRUE, bty='n', adj=c(.1,.4))
legend(x=-.7, y=.67, legend=c('Rain-dominated','Rain-and-snow','Snow-dominated','Dam-influenced'),
       xpd=NA, pt.bg=c('black','gray75','white','white'), pch=c(21,21,21,124),
       col=c('black','black','black','chocolate2'), cex=1.3, horiz=FALSE, bty='n')
legend(x=-.58, y=.595, legend=c(expression(paste(plain('Q, T'[water]))),
                             expression(paste(frac(1,Q), ' ,  ', frac(1,T[water])))),
       xpd=NA, fill=c(adjustcolor('darkred', alpha.f=0.2),
                      adjustcolor('navy', alpha.f=0.15)),
       border='transparent', cex=1.3, horiz=FALSE, bty='n', adj=c(0,.4))
# legend(x=0.08, y=3.75, legend=c(expression(paste(plain('' %up% 'Q, '%up% 'T'[water]))),
#                             expression(paste(plain('' %down% 'Q, '%down% 'T'[water])))),
#        xpd=NA, fill=c(adjustcolor('darkred', alpha.f=0.2),
#                       adjustcolor('navy', alpha.f=0.15)),
#        border='black', cex=1.2, horiz=FALSE, bty='n', adj=c(0,.4))
text(-.66, .52, expression(paste(plain('T'[air]))),
     xpd=NA, cex=1.5)
text(-.59, .52, expression(paste(plain('' %prop% ''))),
     xpd=NA, cex=2.5)
# text(0.01, 3.58, expression(paste(plain('' %up% 'T'[air]) %=>% '')),
#      xpd=NA, cex=1.6)
# text(0.03,23.65,'(Centroids)',xpd=NA)
# legend.gradient(pts=cbind(x=c(-.04,.04,.04,-.04), y=c(22,22,24,24)), cols=c('black','white'),
#                 limits=c('Rainfed', 'Snowfed'), title='X = mean')
par(defpar)
dev.off()
shell('C:\\Users\\Mike\\git\\stream_nuts_DFA\\manuscript\\figures\\15_disch_vs_air_focal5.pdf')

# 20 - air water discharge (average) ####

#load temp_due_5m_atpcsn_byMo_allMos.rda for this one
awd <- readRDS('../../saved_structures/air_water_discharge2.rds')
dis <- readRDS('../../saved_structures/just_discharge.rds')
dams <- read.csv("C:/Users/Mike/git/stream_nuts_DFA/data/watershed_data/watershed_data_simp.csv")[-c(3,18),'dam_upstream']

elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
# cols = rep('black', nrow(land))
# cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'

dam_ind = rep(FALSE,19) #dont use this for other stuff
dam_ind[dams[-c(6,17,18,20,22)] != 0] <- TRUE #or this

# png('12_air_water_discharge.png', width=8, height=6, units='in', res=96, type='cairo')
pdf('12_air_water_discharge.pdf', width=7, height=6)
# pdf('12b_air_water_discharge_expand.pdf', width=7, height=6)
defpar <- par(mar=c(5,4,1,6))
airtemps = watertemps_rain =watertemps_rs = watertemps_snow = discharge = discharge_dam = numeric()
for(i in 1:12){
    airtemps[i] <- mean(awd[[1]][seq(from=i, by=12, to=nrow(awd[[1]]))])
    watertemps_rain[i] <- mean(as.matrix(awd[[2]][seq(from=i, by=12, to=nrow(awd[[1]])),elev_lo]), na.rm=TRUE)
    watertemps_rs[i] <- mean(as.matrix(awd[[2]][seq(from=i, by=12, to=nrow(awd[[1]])),elev_med]), na.rm=TRUE)
    watertemps_snow[i] <- mean(as.matrix(awd[[2]][seq(from=i, by=12, to=nrow(awd[[1]])),elev_hi]), na.rm=TRUE)
    discharge[i] <- mean(as.matrix(dis[seq(from=i, by=12, to=nrow(awd[[1]])),]), na.rm=TRUE)
    discharge_dam[i] <- mean(as.matrix(dis[seq(from=i, by=12, to=nrow(awd[[1]])),dam_ind]), na.rm=TRUE)
}
plot(1:12, discharge, type='n', pch=17, col='gray60', xaxt='n', yaxt='n', xlab='', ylab='',
     # yaxs='i', ylim=c(0,4480), bty='o')
     yaxs='i', ylim=c(0,8000), bty='o')
disch_col = adjustcolor('darkslategray4', alpha.f=0.35)
axis(4, las=1, tck=-.01, hadj=.3, cex.axis=.8, col.axis=adjustcolor('darkslategray4', alpha.f=1),
     col.ticks=adjustcolor('darkslategray4', alpha.f=.7))
# axis(4, tcl=0, col='white', labels='')
polygon(x=c(0,.57,1:12,12.5,13), y=c(0,mean(c(discharge[1],discharge[12])),discharge,mean(c(discharge[1],discharge[12])),0), col=disch_col, border=NA)
# lines(c(.58,1:12,12.4), c(mean(c(discharge_dam[1],discharge_dam[12])),
#                           discharge_dam,mean(c(discharge_dam[1],discharge_dam[12]))),
#       lwd=3, col='cadetblue', lty=1)
lines(x=c(12.44,12.44),y=c(0,4018),
                         col='white', xpd=NA, lwd=1)
lines(x=c(12.44,12.44),y=c(0,4018),
                         col=adjustcolor('darkslategray4', alpha.f=.35), xpd=NA, lwd=1)
mtext('Mean discharge (CFS)', 4, font=2, las=3, line=2.5, col=adjustcolor('darkslategray4', alpha.f=1))
mtext(expression(paste(~bold('Mean temperature')~bold('(')*bold(degree)*bold('C)'))),
      line=1.4, side=2)
for(i in 1:12) lines(x=c(i,i),y=c(0,discharge[i]), , col='white', lty='29')
par(new=T)
plot(1:12, airtemps, xaxt='n', xlab='', bty='c',
     ylab='', yaxs='i', type='n', ylim=c(3,17.6), yaxt='n')
mtext(expression(paste(~bold('Month'))), 1, line=1.5)
axis(2, las=2, at=c(10,12,14,16), cex.axis=.8, hadj=.3, tck=-.01)
axis(2, las=2, at=c(4,6,8), cex.axis=.8, hadj=-.5, tck=-.01)
# for(i in 1:12) lines(x=c(i,i),y=c(0,max(watertemps_rain[i],watertemps_rs[i],watertemps_snow[i])), lty='19')
# for(i in c(1:6,11,12)){
#     lines(x=c(i,i),y=c(max(watertemps_rain[i],watertemps_rs[i],watertemps_snow[i]),
#                        min(discharge[i])), lty='29', col='white')
# }
lines(c(.58,1:12,12.4), c(mean(c(airtemps[1],airtemps[12])),airtemps,mean(c(airtemps[1],airtemps[12]))),
      lwd=3, col='chocolate2', lty='33')
points(1:12, watertemps_rain, col='black', bg='black', pch=21, cex=1.3)
points(1:12, watertemps_rs, bg='gray75', col='black', pch=22, cex=1.3)
points(1:12, watertemps_snow, bg='white', col='black', pch=24, cex=1.3, lwd=1)
axis(1, at=1:12, labels=month.abb, cex.axis=.8, padj=-1.5, tck=-.01)
legend(x=1.6, y=17.4, legend=c('Rain-dominated', 'Rain-and-snow', 'Snow-dominated'),
       pch=c(21,22,24), col='black', bty='n', pt.bg=c('black','gray75','white'), cex=.8)
legend(2.2,15.9,legend='Air', lty='33', lwd=2, col='chocolate2', bty='n', seg.len=2.5, x.intersp=.8, cex=.8)

par(defpar)
dev.off()
shell('C:\\Users\\Mike\\git\\stream_nuts_DFA\\manuscript\\figures\\12_air_water_discharge.pdf')

#for results section:
# chili = data.frame(watertemps_rain,watertemps_rs,watertemps_snow)
# chili2 = apply(chili,1,range)
# chili3 = diff(chili2)
min(chili3); max(chili3)

# 21 - air water discharge (over time [never actually worked on this]) ####

#load temp_due_4m_at_byMo_allMos.rda for this one
awd <- readRDS('../../saved_structures/air_water_discharge.rds')

# png('12_air_water_discharge.png', width=8, height=6, units='in', res=96, type='cairo')
# pdf('12_air_water_discharge.pdf', width=8, height=6)
defpar <- par(mar=c(5,4,1,6))
airtemps = watertemps = discharge = numeric()
for(i in 1:12){
    airtemps[i] <- mean(awd[[1]][seq(from=i, by=12, to=nrow(awd[[1]]))])
    watertemps[i] <- mean(as.matrix(awd[[2]][seq(from=i, by=12, to=nrow(awd[[1]])),]), na.rm=TRUE)
    discharge[i] <- mean(as.matrix(obs_ts[seq(from=i, by=12, to=nrow(awd[[1]])),]), na.rm=TRUE)
}
plot(1:12, airtemps, xaxt='n', xlab=expression(paste(~bold('Month'))),
     ylab='')
points(1:12, watertemps, pch=20)
axis(1, at=1:12, labels=month.abb)
legend('left', legend=c('Air', 'Water', 'Discharge'), pch=c(1,20,17), bty='n',
       col=c('black','black','gray60'))
par(new=T)
plot(1:12, discharge, type='b', pch=17, col='gray60', xaxt='n', yaxt='n', xlab='', ylab='')
axis(4, las=1)
mtext('Mean discharge (CFS)', 4, font=2, las=3, line=3)
mtext(expression(paste(~bold('Mean temperature')~bold(degree)~bold('C'))),
      line=2, side=2)
par(defpar)
# dev.off()

# 22 - discharge eff size and loading regressions separated [dont think i worked on this either] ####
# pdf('01b_discharge_effect_size_reg.pdf', width=8, height=3.75)
# layout(matrix(c(1:5,8,6,7,9),nrow=3,byrow=TRUE))
# par(oma=c(
landvar=land$ElevWs/100
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
defpar <- par(oma=c(2,4,2,1), mar=c(4,2,2,1))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
for(covr in 1:3){
    plot(landvar, rescaled_effect_size[,covr], type='n', yaxt='n',
         # xlab='Watershed area over 1000m (%)', main='',
         xlab='', main='', ylab='', cex.lab=1.3, cex.axis=1, font=2, xaxt='n',bty='o')
    mod <- lm(rescaled_effect_size[,covr] ~ I(landvar))
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    if(covr == 1) mtext(bquote(plain(theta) == bold(T[air]) ~.(sig)), 3, cex=1.1, font=2, line=1.2)
    if(covr == 2){
        mtext(bquote(plain(theta) == bold('Precip') ~.(sig)), 3, cex=1.1, font=2, line=1.2)
        abline(mod, col='steelblue', lty=2, lwd=3)
    }
    if(covr == 3){
        mtext(bquote(plain(theta) == bold('Snowmelt') ~.(sig)), 3, cex=1.1, font=2, line=1.2)
        abline(mod, col='steelblue', lty=2, lwd=3)
    }
    points(landvar, rescaled_effect_size[,covr],
           bg=cols, col='black',
           pch=21, cex=2, lwd=2)
    # color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
    #              rect.col=colorRampPalette(c('brown', 'white'))(10),
    #              align='r', gradient='y')
    # text(3.39, 0.56, labels='Mean watershed')
    # text(3.52, 0.54, labels='elevation (m)')
    axis(1); axis(2)
}
mtext('Mean watershed elevation (100 m)', side=1, outer=TRUE, cex=1.2)
mtext(expression(paste(Delta,'Q ', Delta, theta^-1)), side=2, outer=TRUE, cex=1.3, line=1)
# par(defpar)
# dev.off()

# 10 - DISCHARGE loading regression

# png('02_loadings_reg.png', width=7, height=6, units='in', res=96, type='cairo')
# pdf('02_loadings_reg.pdf', width=7, height=6)
# defpar <- par(mar=c(5,5,2,5))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
# plot(land$Ice06_11, dfa$Estimates$Z[,1],
#      xlab='Watershed % ice cover', ylab='Trend 1 loadings', type='n',
#      main='', cex.lab=1.3, cex.axis=1, font=2)
# abline(mod, col='gray', lty=2, lwd=3)
# points(land$Ice06_11, dfa$Estimates$Z[,1], col='black', bg=cols,
#        pch=21, cex=1.5, lwd=2)
# mod <- lm(dfa$Estimates$Z[,1] ~ land$Ice06_11)
# color.legend(xl=4,xr=4.4,yb=1.32, yt=1.82, legend=c('147', '1349'),
#              rect.col=colorRampPalette(c('brown', 'white'))(10),
#              align='r', gradient='y')
# text(3.39, 1.61, labels='Mean watershed')
# text(3.52, 1.53, labels='elevation (m)')
# par(defpar)
# dev.off()

# pdf('02b_discharge_loadings_reg.pdf', width=7, height=7)
landvar = list(NULL,land$WsSlope, land$RunoffWs, land$WsAreaSqKm, land$RckDepWs)
land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
#% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)
# png('01_effect_size_reg.png', width=7, height=6, units='in', res=96, type='cairo')
defpar <- par(mfrow=c(2,2), oma=c(0,0,0,0), mar=c(4,4,1,1))
# pal <- colorRampPalette(c('brown', 'white'))
# cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
for(trnd in 2:5){
    plot(landvar[[trnd]], dfa$Estimates$Z[,trnd], type='n', yaxt='n',
         xlab='', main='', ylab='', cex.lab=1, font=2, xaxt='n',bty='o')
    mod <- lm(dfa$Estimates$Z[,trnd] ~ landvar[[trnd]])
    sig <- ifelse(summary(mod)$coefficients[2,4]<=0.1, '*', '')
    if(trnd == 2) mtext(bquote(bold('Mean slope ')*plain('(% rise)')~.(sig)), 1, cex=1, font=2, line=2.7)
    if(trnd == 3) mtext(bquote(bold('Total runoff ')*plain('(mm ')*plain(mo^-1)*plain(')')~.(sig)), 1, cex=1, font=2, line=2.7)
    if(trnd == 4) mtext(bquote(bold('Total area ')*plain('(')*plain(km^2)*plain(')')~.(sig)), 1, cex=1, font=2, line=2.7)
    if(trnd == 5) mtext(bquote(bold('Mean bedrock depth ')*plain('(cm)')~.(sig)), 1, cex=1, font=2, line=2.7)
    # if(trnd == 3){
    #     mtext(bquote(bold('Mean slope (% rise)')~.(sig)), 1, cex=1, font=2, line=2.7)
    #     # mtext(bquote(plain(theta) == bold('Precip') ~.(sig)), 3, cex=1.1, font=2, line=1.2)
    # }
    abline(mod, col='springgreen4', lty=2, lwd=3)
    points(landvar[[trnd]], dfa$Estimates$Z[,trnd],
           bg=cols, col='black',
           pch=21, cex=2, lwd=2)
    # color.legend(xl=4,xr=4.4,yb=0.5, yt=0.6, legend=c('147', '1349'),
    #              rect.col=colorRampPalette(c('brown', 'white'))(10),
    #              align='r', gradient='y')
    # text(3.39, 0.56, labels='Mean watershed')
    # text(3.52, 0.54, labels='elevation (m)')
    axis(1); axis(2)
}
mtext('Shared trends', side=2, outer=TRUE, cex=1, line=1, font=2)

par(defpar)
dev.off()

# 23 - elwha before and after (temp) ####

pre = obs_ts[1:396,23] #elwha temp data through 2010 (restoration began sept 2011)
post2011 = obs_ts[397:408,23]
post2012 = obs_ts[409:420,23]
post2013 = obs_ts[421:432,23]
post2014 = obs_ts[433:444,23]
post2015 = obs_ts[445:456,23] #2015 only (ended aug 2014)
prevec  = 1:12
for(i in 1:12){
    prevec[i] = mean(pre[seq(i,396-(12-i),12)], na.rm=TRUE)
}
plot(prevec, type='l')
for(i in 1:5){
    # browser()
    chili = get(paste0('post201',i))
    lines(chili, col=heat.colors(5)[i])
}

# 24 - elwha before and after (discharge) ####

pre = obs_ts[1:396,18] #elwha temp data through 2010 (restoration began sept 2011)
post2011 = obs_ts[397:408,18]
post2012 = obs_ts[409:420,18]
post2013 = obs_ts[421:432,18]
post2014 = obs_ts[433:444,18]
post2015 = obs_ts[445:456,18] #2015 only (ended aug 2014)
prevec  = 1:12
for(i in 1:12){
    prevec[i] = mean(pre[seq(i,396-(12-i),12)], na.rm=TRUE)
}
plot(prevec, type='l', ylim=c(0,4000))
for(i in 1:5){
    browser()
    chili = get(paste0('post201',i))
    lines(chili, col=heat.colors(5)[i])
}
