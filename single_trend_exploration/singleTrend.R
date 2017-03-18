#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 3/17/17
rm(list=ls()); cat('\014')

setwd('~/git/puget_sound_rivers_DFA/single_trend_exploration')

# load('1trend.rda')
# load('1trendNoSeas.rda')
load('2trendNoSeas.rda')

ice2006 <- read.csv('../manuscript/figures/PctIce2006Ws.csv', stringsAsFactors=FALSE)
land <- merge(land,ice2006[,2:3], by.x='siteCode', by.y='site')
land$Ice06_11 <- rowMeans(cbind(land$PctIce2006Ws, land$PctIce2011Ws))

# pcascores <- read.csv('../data/watershed_data/pca_scores.csv', stringsAsFactors=FALSE)
# land <- merge(land,pcascores, by.x='siteCode', by.y='X')

process_plotter_TMB(dfa, mm, chunk=1)
loading_plotter_TMB(dfa, mm)

defpar = par(mfrow=c(2,1))
ccf(cov_and_seas[1,], as.vector(dfa$Estimates$u[1,])) #13
ccf(cov_and_seas[1,], as.vector(dfa$Estimates$u[2,])) #13
ccf(colSums(cov_and_seas[3:14,]), as.vector(dfa$Estimates$u[1,])) #15:26
ccf(colSums(cov_and_seas[3:14,]), as.vector(dfa$Estimates$u[2,])) #15:26
par(defpar)

load_regress_plotter(mm, 'indiv', 'PC2', 'WsAreaOver1000')
load_regress_plotter(mm, 'exploration', , 'ElevWs')
