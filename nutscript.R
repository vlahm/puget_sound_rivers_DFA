setwd("C:/Users/Mike/Desktop/Grad/Projects/Thesis/stream_nuts_DFA/data")

nuts <- read.csv('nutsdata.csv')
str(nuts)

nuts$dateTime[head(which(!is.na(nuts$NH3_N)), 1)] #first NH3 reading
nuts$dateTime[head(which(!is.na(nuts$TP_P)), 1)] #first TP reading

plot(nuts$NH3_N, type='l', col='red', na.rm=T)
plot(nuts$TP_P, type='l', col='blue')
