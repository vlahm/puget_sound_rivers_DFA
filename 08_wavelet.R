#Puget Sound rivers DFA
#Mike Vlah (vlahm13@gmail.com)
#created: 5/28/17
#wavelet analysis and determining variance explained by trends

rm(list=ls()); cat('\014') #clear env and console

# 0 - setup ####

setwd('C:/Users/Mike/git/stream_nuts_DFA/')
setwd('~/git/puget_sound_rivers_DFA/')

if(!require(WaveletComp)) install.packages('WaveletComp')
if(!require(imputeTS)) install.packages('imputeTS')
library(WaveletComp)
library(imputeTS)

# load('discharge_due_5m_atpcsn_byMo_allMos.rda') #best discharge model
# load('temp_due_5m_atpcsn_byMo_allMos.rda') #best temp model
load('single_trend_exploration/2trendNoSeasNoSnow.rda') #simplified 2-trend temp model
# load('manuscript/figures/temp_due_5m_atpcsn_byMo_allMos.rda') #best temp model
# load('manuscript/figures/discharge_due_5m_atpcsn_byMo_allMos.rda') #best discharge model

# 1 - determine variance explained by trends and full model ####

trendFit = dfa$Estimates$Z %*% dfa$Estimates$u
covFit = dfa$Estimates$D %*% cov_and_seas
fullFit = dfa$Fits
data = dat_z

get_R2 <- function(fit){
    R2 <- rep(NA, nrow(data))
    for(i in 1:nrow(data)){
        mod <- lm(data[i,] ~ fit[i,])
        R2[i] <- summary(mod)$r.squared
    }
    return(list(mean=mean(R2), sd=sd(R2)))
}

fullR2 = get_R2(fullFit)
trendR2 = get_R2(trendFit)
covR2 = get_R2(covFit)

# 2 - wavelet ####

#testing on raw temp data
# testdata = data.frame('date'=yy$date, 'test'=na.seasplit(data[1,], 'interpolation'))
# wav_test = analyze.wavelet(testdata, 'test', dt=1)
# wt.image(wav_test, n.levels=100)

#for reals
trends = data.frame(date=yy$date, trnd1=dfa$Estimates$u[1,], trnd2=dfa$Estimates$u[2,])
wav1 = analyze.wavelet(trends, 'trnd1', dj=1/200)
wav2 = analyze.wavelet(trends, 'trnd2', dj=1/200)

#plot wavelet power spectra
# png('24_wav1.png', width=8, height=6, units='in', res=300, type='cairo')
pdf('manuscript/figures/24_wav1.pdf', width=7, height=6)
wt.image(wav1, show.date=TRUE)
dev.off()
# png('24b_wav2.png', width=8, height=6, units='in', res=300, type='cairo')
pdf('manuscript/figures/24b_wav2.pdf', width=7, height=6)
wt.image(wav2, show.date=TRUE)
dev.off()

#clean trend 2
rec = reconstruct(wav2, sel.period=12, legend.coords = "bottomleft")
cleaned = rec$series$trnd2.r
plot(cleaned, type='l')

pdf('manuscript/figures/25_wavClean.pdf', width=7, height=6)
ind = rep(1:12, times=38)
year_mean = tapply(cleaned, ind, mean)
year_sd = tapply(cleaned, ind, sd)
plot(year_mean, type='n', xaxs='i',
     ylim=c(min(year_mean-year_sd),max(year_mean+year_sd)),
     ylab=expression(paste(bold('Standardized')~bold(T[water])~bold('(')*bold(degree)*bold('C)'))),
     xaxt='n', las=2, xlab=expression(bold('Time')))
axis(1, 1:12, labels=month.abb)
polygon(x=c(1:12,12:1), y=c(year_mean+year_sd, rev(year_mean-year_sd)),
        col='gray85', border=NA)
lines(year_mean, col='blue', lwd=2)
dev.off()
