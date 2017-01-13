#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 1/13/17
#plots of best temperature model: temp_due_2m_fixed_at_1978-2015

rm(list=ls()); cat('\014') #clear env and console

# 0 - setup ####

setwd('C:/Users/Mike/git/stream_nuts_DFA/manuscript/')
setwd('~/git/puget_sound_rivers_DFA/manuscript')
setwd('Z:/stream_nuts_DFA/manuscript/')

library(viridis)

load('Renv_image2.rda') #load all objects generated during the creation and postprocessing
                       #of the best model
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
land$siteCode

# 1 - effect size regression ####

# land_sub <- land[,landcols] #subset landscape variables by those used in the analysis

#% of watershed area classified as ice/snow land cover (NLCD 2011 class 12)
#% of watershed area classified as ice/snow land cover (NLCD 2006 class 12)

pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(10)[as.numeric(cut(land$ElevWs, breaks=10))]
plot(land$Ice06_11, rescaled_effect_size,
     xlab='Watershed % ', main='', col=cols,
     ylab=expression(paste(Delta, ' stream temp (', degree, 'C) / ',
                            Delta, ' air temp (', degree, 'C)')),
     pch=colnames(trans$trans))
mod <- lm(rescaled_effect_size ~ land_sub[,names(best)[1:4][i]])
abline(mod, col='gray', lty=2)

#why is water table depth such a strong factor? what else is it correlated with?
rev(tail(sort(abs(apply(land[,43:ncol(land)], 2, function(x) cor(land$WtDepWs, x)))), 15))
#it's just elevation (note the returned correlations have been abs()'d)

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

#plot best cors along with fitted models
defpar <- par(mfrow=c(3,2))

pal <- colorRampPalette(c('blue', 'red'))
cols <- pal(10)[as.numeric(cut(land_sub$ElevWs, breaks=10))]
# full_names <- c('organic matter', 'base flow', 'coastal alluvium', 'runoff',
#                 'riparian urbanization (high)', 'riparian urbanization (low)')
for(i in 1:6){
    plot(land_sub[,names(best1)[i]], dfa$Estimates$Z[,1],
         xlab=names(best1)[i], ylab='loading on common trend 1',
         main='blue=low elev, red=high elev', col=cols, pch=colnames(trans$trans))
    mod <- lm(dfa$Estimates$Z[,1] ~ land_sub[,names(best1)[i]])
    abline(mod, col='gray', lty=2)
}

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
    mod <- lm(rescaled_seas[,i] ~ land$Ice06_11, weights=land$watershedA)
    slope <- round(unname(mod$coefficients[2]), 2)
    plot(land$Ice06_11, rescaled_seas[,i], main=paste(month.abb[i], 'slope =', slope),
         ylab=paste(month.abb[i], 'change in water temp'), xlab='% ice',
         ylim=c(min(rescaled_seas), max(rescaled_seas)), col=cols, cex=1,
         pch=colnames(trans$trans))
    abline(mod, col='gray', lty=2, lwd=2)
    abline(h=0, col='red')
}
par(defpar)
# dev.off()
