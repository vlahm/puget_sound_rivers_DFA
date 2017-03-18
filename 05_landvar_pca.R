#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 5/17/16
rm(list=ls()); cat('\014') #clear env and console

# 0 - setup ####
setwd('C:/Users/Mike/git/stream_nuts_DFA/manuscript/figures/diagnostic_plots')
setwd('~/git/puget_sound_rivers_DFA/manuscript/figures/diagnostic_plots')

# load packages and some functions
library(vegan)
source('http://www.umass.edu/landeco/teaching/ecodata/labs/biostats.R')

#load desired model output and accessories
load('1trend.rda')
# load('2trend.rda')
# load('1trendNoSeas.rda')
# load('2trendNoSeas.rda')

#add % ice to watershed data
ice2006 = read.csv('../PctIce2006Ws.csv', stringsAsFactors=FALSE)
land = merge(land,ice2006[,2:3], by.x='siteCode', by.y='site')
land$Ice06_11 = rowMeans(cbind(land$PctIce2006Ws, land$PctIce2011Ws))

#compile watershed variables for pca
pcanames = c('BFIWs','ElevWs','PctImp2006WsRp100','RunoffWs',
              'WtDepWs','PermWs', 'WsAreaSqKm','WsAreaOver1000','WsSlope','Ice06_11')

pcalandcols = rep(NA, length(pcanames))
for(i in 1:length(pcanames)){
    pcalandcols[i] = which(colnames(land) == pcanames[i])
}

pcavars = land[,pcalandcols]
row.names(pcavars) = land$siteCode

#transform proportion and percent data
pcavars[,'PctImp2006WsRp100'] = asin(sqrt(pcavars[,'PctImp2006WsRp100']/100))
pcavars[,'WsAreaOver1000'] = asin(sqrt(pcavars[,'WsAreaOver1000']))
pcavars[,'Ice06_11'] = asin(sqrt(pcavars[,'Ice06_11']/100))

# 1 - run pca ####

pca_out = prcomp(pcavars, scale=TRUE)
summary(pca_out)

#variance explained by each PC
pca.eigenval(pca_out)

#check significance of eigenvalues (of PCs)
ordi.monte(pcavars, ord='pca', dim=5) #p vals < 0.05 are sig

#check eigenvectors (variable loadings on each PC) - high abs vals show high contribution
pca.eigenvec(pca_out,dim=5,digits=3,cutoff=.3)

#convert eigenvector coefficients to correlation coefficients
pca.structure(pca_out,pcavars,dim=5,cutoff=.4)
    #these are correlations between the original variables and the PC scores
    #squaring these gives percentage of variance in each original var accounted for by each PC
    #use these or the eigenvectors to generate an ecological interpretation of each PC

#get sample scores
samp_scores = pca_out$x
    #standardized scores for each sample on each PC axis - they represent the position
    #of each sample on each standardized PC axis

# 3 - output ####

#plot
pdf('../17_pca.pdf', width=5, height=5)
ordiplot(pca_out, choices = c(1, 2), type="text", display='sites',
         xlab='PC 1 (50.5%)', ylab='PC 2 (21.7%)', ylim=c(-3,3.5))
arrows(0,0,pca_out$rotation[,1]*5,pca_out$rotation[,2]*5,col='blue',length=.05)
text(pca_out$rotation[,1]*5.4,pca_out$rotation[,2]*5.4,
     row.names(pca_out$rotation), col='blue')
dev.off()

#export sample scores
write.csv(samp_scores[,1:2], '../pca_scores.csv')
