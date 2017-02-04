setwd('C:/Users/Mike/git/stream_nuts_DFA/data/watershed_data')

x <- read.csv('wsAreas.csv')

elev_inds <- which(!grepl('COUNT', colnames(x)))
names <- colnames(x)[elev_inds]

out <- data.frame(siteCode=names, WsAreaOver1000=NA)

for(i in 1:length(elev_inds)){
    over1000 <- which(x[,elev_inds[i]] > 1000)
    areaOver1000 <- sum(x[over1000,elev_inds[i]+1])
    propOver1000 <- areaOver1000/sum(x[,elev_inds[i]+1], na.rm=TRUE)
    out$WsAreaOver1000[i] <- propOver1000
}

write.csv(out, 'WsAreaOver1000.csv', row.names=FALSE)
