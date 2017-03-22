#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 5/17/16
rm(list=ls()); cat('\014') #clear env and console

# 0 - setup and preprocessing ####
setwd('C:/Users/Mike/git/stream_nuts_DFA')

# load packages and some functions
library(vegan)
library(cluster)
source('http://www.umass.edu/landeco/teaching/ecodata/labs/biostats.R')

#load desired model output and accessories
# load('1trend.rda')
# load('2trend.rda')
# load('1trendNoSeas.rda')
load('single_trend_exploration/2trendNoSeas.rda')

#add % ice and dams to watershed data
ice2006 = read.csv('manuscript/figures/PctIce2006Ws.csv', stringsAsFactors=FALSE)
land = merge(land,ice2006[,2:3], by.x='siteCode', by.y='site')
land$Ice06_11 = rowMeans(cbind(land$PctIce2006Ws, land$PctIce2011Ws))

dams = read.csv('data/watershed_data/watershed_data_simp.csv', stringsAsFactors=FALSE)[,c('siteCode','dam_upstream')]
dams$siteCode[dams$siteCode=='AA'] <- 'ZA'
land = merge(land, dams, by='siteCode')

#compile watershed variables for pca
ordnames = c('BFIWs','ElevWs','PctImp2006WsRp100','RunoffWs',
             'WtDepWs','PermWs', 'WsAreaSqKm','WsAreaOver1000','WsSlope','Ice06_11','dam_upstream')

ordlandcols = rep(NA, length(ordnames))
for(i in 1:length(ordnames)){
    ordlandcols[i] = which(colnames(land) == ordnames[i])
}

ordvars = land[,ordlandcols]
row.names(ordvars) = land$siteCode

#transform proportion and percent data
ordvars[,'PctImp2006WsRp100'] = asin(sqrt(ordvars[,'PctImp2006WsRp100']/100))
ordvars[,'WsAreaOver1000'] = asin(sqrt(ordvars[,'WsAreaOver1000']))
ordvars[,'Ice06_11'] = asin(sqrt(ordvars[,'Ice06_11']/100))

#get Gower's distance (standardization is done internally)
ordvars$dam_upstream = factor(ordvars$dam_upstream)
orddist = daisy(ordvars, metric='gower')

# 1 - run PCoA ####

#perform pcoa
ord_out = cmdscale(orddist, eig=TRUE, add=FALSE)
#get scores (PCs) and eigenvalues
scores = ord_out$points
colnames(scores) = c('PCo1','PCo2')
ord_out$eig
#get percent variation explained by each principal coordinate
opt = options(scipen=100)
ord_out$eig/sum(ord_out$eig)*100 #(the first 5 are for the first 5 principal coords, etc.)
options(opt)
#broken stick (compare eigenvalues to expectations)
plot(ord_out$eig[1:35]/sum(ord_out$eig)*100,type="b",lwd=2,col="blue",
     xlab= "Principal Component from PCoA", ylab="% variation explained",
     main="% variation explained by PCoA (blue) vs. random expectation (red)")
lines(bstick(35)*100,type="b",lwd=2,col="red")

#loadings (i.e. principal weights in the eigenvectors) on each principal coord
ord_load = envfit(ord_out$points, k=45, ordvars, perm=1000)
ord_load$vectors$arrows

#plot pcoa
pdf('../18_PCoA.pdf', width=5, height=5)
ordiplot(ord_out, choices = c(1, 2), type="text", display='sites',
         xlab='PC 1 (51.0%)', ylab='PC 2 (30.0%)')
plot(ord_load, p.max=.05, col='blue')
dev.off()

#export sample scores
write.csv(scores, 'data/watershed_data/pcoa_scores.csv')
