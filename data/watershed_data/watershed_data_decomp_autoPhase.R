sites <- read.csv('../file_ffs2_csv.csv')

files <- dir(pattern='.csv')
obj_names <- vector(length=length(files))
for(i in 1:length(files)){
  obj_names[i] <- substr(files[i], 1, nchar(files[i])-4)
  temp <- read.csv(files[i])
  assign(obj_names[i], temp)
}

out <- Reduce(function(x,y){merge(x,y, by='COMID')}, list(sites,
        BFI_Region17,Elevation_Region17,EPA_FRS_Region17,ImperviousSurfaces2006_Region17,
        ImperviousSurfaces2006RipBuf100_Region17,Lithology_Region17,NADP_Region17,NLCD2011_Region17,
        NLCD2011RipBuf100_Region17,RoadDensity_Region17,RoadDensityRipBuf100_Region17,
        Runoff_Region17,STATSGO_Set2_Region17,USCensus2010_Region17,USCensus2010RipBuf100_Region17))

write.csv(t(out), "C:/Users/vlahm/Desktop/sites.csv")

#in excel I:
#removed .x from first instance of 
#CatAreaSqKm.x WsAreaSqKm.x CatPctFull.x WsPctFull.x AND CatAreaSqKmRp100.x WsAreaSqKmRp100.x CatPctFullRp100.x WsPctFullRp100.x
#deleted all the rest


colnames(out)[duplicated(colnames(out))]
