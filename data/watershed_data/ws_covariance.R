#this script was produced after making the following choices
y_choice = 'TEMP'
cov_choices = c('meantemp')
region = '3_4'
average_regions = TRUE
method = 'fixed_individual'
startyr = 1978
endyr = 2015
#and running 03_setup_bymonth.R down to where obs_ts is generated



#organize landscape data
land <- read.csv('watershed_data/watershed_data_simp.csv', stringsAsFactors=FALSE)
land$siteCode[land$siteCode == 'AA'] <- 'ZA' #rename sammamish @ bothell sitecode
land <- land[land$siteCode %in% names(obs_ts),] #remove sites not in analysis
land <- land[match(names(obs_ts), land$siteCode),] #sort landscape data by site order in model

#round 1 ####
#identify potentially informative landscape variables
landvars <- c('BFIWs','ElevWs','PctImp2006Ws','PctImp2006WsRp100',
              'PctGlacTilLoamWs','PctGlacLakeFineWs',
              'PctAlluvCoastWs','PctIce2011Ws','PctUrbOp2011Ws','PctUrbLo2011Ws',
              'PctUrbMd2011Ws','PctUrbHi2011Ws','PctHay2011Ws','PctCrop2011Ws',
              'PctIce2011WsRp100','PctUrbOp2011WsRp100','PctUrbLo2011WsRp100',
              'PctUrbMd2011WsRp100','PctUrbHi2011WsRp100','PctHay2011WsRp100',
              'PctCrop2011WsRp100','RdDensWs','RdDensWsRp100','RunoffWs','OmWs',
              'RckDepWs','WtDepWs','PermWs','HUDen2010Ws','PopDen2010Ws','PopDen2010WsRp100')
landcols <- which(colnames(land) %in% landvars)

#get abs correlations between landscape vars
landcors <- round(abs(cor(land[landcols])), 2)
#replace low correlatios with NA for visual inspection
landcors[landcors < .6] <- NA
landcors
#remove some of the less interesting variables that covary strongly with others and repeat

#round 2 ####
#identify potentially informative landscape variables
landvars <- c('BFIWs','ElevWs','PctImp2006WsRp100',
              'PctGlacTilLoamWs','PctGlacLakeFineWs',
              'PctAlluvCoastWs','PctIce2011Ws',
              'PctHay2011Ws','PctCrop2011Ws',
              'PctIce2011WsRp100','PctUrbOp2011WsRp100','PctUrbLo2011WsRp100',
              'PctUrbMd2011WsRp100','PctUrbHi2011WsRp100','PctHay2011WsRp100',
              'PctCrop2011WsRp100','RdDensWsRp100','RunoffWs','OmWs',
              'RckDepWs','WtDepWs','PermWs','PopDen2010Ws','PopDen2010WsRp100')
landcols <- which(colnames(land) %in% landvars)

#get abs correlations between landscape vars
landcors <- round(abs(cor(land[landcols])), 2)
#replace low correlatios with NA for visual inspection
landcors[landcors < .6] <- NA
landcors
#remove some of the less interesting variables that covary strongly with others and repeat

#round 3 ####
#identify potentially informative landscape variables
landvars <- c('BFIWs','ElevWs','PctImp2006WsRp100',
              'PctGlacLakeFineWs',
              'PctAlluvCoastWs','PctIce2011Ws',
              'PctCrop2011Ws',
              'PctUrbOp2011WsRp100','PctUrbLo2011WsRp100',
              'PctUrbMd2011WsRp100','PctUrbHi2011WsRp100',
              'RdDensWsRp100','RunoffWs','OmWs',
              'RckDepWs','WtDepWs','PermWs','PopDen2010Ws')
landcols <- which(colnames(land) %in% landvars)

#get abs correlations between landscape vars
landcors <- round(abs(cor(land[landcols])), 2)
#replace low correlatios with NA for visual inspection
landcors[landcors < .6] <- NA
landcors

#done. these are the vars that end up in the analysis
