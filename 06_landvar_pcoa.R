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
# load('1trend.rda') #doesnt have aspect added to landscape vars
# load('2trend.rda') #doesnt have aspect added to landscape vars
# load('1trendNoSeas.rda') #doesnt have aspect added to landscape vars
load('single_trend_exploration/2trendNoSeasNoSnow.rda')

#add % ice and dams to watershed data
ice2006 = read.csv('manuscript/figures/PctIce2006Ws.csv', stringsAsFactors=FALSE)
land = merge(land,ice2006[,2:3], by.x='siteCode', by.y='site')
land$Ice06_11 = rowMeans(cbind(land$PctIce2006Ws, land$PctIce2011Ws))

# dams = read.csv('data/watershed_data/watershed_data_simp.csv', stringsAsFactors=FALSE)[,c('siteCode','dam_upstream')]
# dams$siteCode[dams$siteCode=='AA'] <- 'ZA'
# land = merge(land, dams, by='siteCode')

#compile watershed variables for pca
ordnames = c('BFIWs','ElevWs','PctImp2006WsRp100','RunoffWs',
             'WtDepWs','PermWs', 'WsAreaSqKm','WsAreaOver1000','WsSlope',
             'Ice06_11','dam_upstream','aspect')

ordlandcols = rep(NA, length(ordnames))
for(i in 1:length(ordnames)){
    ordlandcols[i] = which(colnames(land) == ordnames[i])
}

ordvars = land[,ordlandcols]
row.names(ordvars) = land$siteCode

#transform proportion and percent data (choosing asin transform over logit to avoid unrealistic -inf and inf
ordvars$PctImp2006WsRp100 = asin(sqrt(ordvars$PctImp2006WsRp100/100))
ordvars$WsAreaOver1000 = asin(sqrt(ordvars$WsAreaOver1000))
ordvars$Ice06_11 = asin(sqrt(ordvars$Ice06_11/100))
x = ordvars$aspect
ordvars$aspect = asin(sqrt((x-min(x))/(max(x)-min(x))))
ordvars$dam_upstream = factor(ordvars$dam_upstream)

#get Gower's distance (standardization is done internally)
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

#rename variable IDs
row.names(ord_load$vectors$arrows) <- c('BFI','elev.','imperv.','runoff','water table depth',
                                        'soil perm.','area','% area > 1000 m','slope',
                                        '% ice cover','aspect')
row.names(ord_load$factors$centroids) <- c('','','(dam-influenced)')
# row.names(ord_load$factors$centroids) <- c('','','')

#plot pcoa
pdf('manuscript/figures/18_PCoA.pdf', width=5, height=5)
elev_hi = which(land$Ice06_11 >= 0.7)
elev_med = which(land$Ice06_11 < 0.7 & land$ElevWs > 600)
elev_lo = which(land$Ice06_11 < 0.7 & land$ElevWs <= 600)
cols = rep('black', nrow(land))
cols[elev_med] = 'gray60'; cols[elev_hi] = 'white'
dam_col = rep('black',nrow(land))
dam_col[land$dam_upstream != 0] <- 'chocolate2'
fig = ordiplot(ord_out, choices = c(1, 2), type="none", display='sites', bty='l',
         xlab=bquote(textstyle(bold('PC 1 '))*textstyle(plain('(50.8%)'))), 
         ylab=bquote(textstyle(bold('PC 2 '))*textstyle(plain('(28.6%)'))),
         xlim=c(-.4,.5))
plot(ord_load, p.max=.05, col='steelblue', font=2)
points(fig, 'sites', pch=21, col=dam_col, bg=cols, cex=2, lwd=2)
arrows(x0=0, x1=-.025, y0=0, y1=0.16, length=.07, col='steelblue', lty=2)
dev.off()

#export sample scores
write.csv(scores, 'data/watershed_data/pcoa_scores.csv')
