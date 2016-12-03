#run 03_setup_bymonth.R through step 3.1 before using this script

setwd('C:/Users/Mike/git/stream_nuts_DFA/boxcox_testing')
# source('../00_tmb_uncor_Rmat.R')

#if there are any timepoints with 0 or 1 observation, remove them
propNA <- apply(dat_z, 2, function(x) sum(is.na(x))/length(x))
highNA_ind <- as.numeric(names(propNA[propNA >= 1-1/nrow(dat_z)]))

#run several DFAs to determine the relationship between effect sizes with and without boxcox transform
#here testing COND, TEMP
if(length(highNA_ind)){
    mod_covs <- t(as.matrix(cov_and_seas[13,-highNA_ind]))

    dfa <- runDFA(obs=dat_z[,-highNA_ind], NumStates=1, ErrStruc='DUE',
              EstCovar=TRUE, Covars=mod_covs)
} else {
    mod_covs <- t(as.matrix(cov_and_seas[13,]))

    dfa <- runDFA(obs=dat_z, NumStates=1, ErrStruc='DUE',
              EstCovar=TRUE, Covars=mod_covs)
}

modified_eff_rescaler <- function(all_cov, seas){
    #get covariate effect sizes (D coefficients) from model, isolated from seasonal effects
    if(nrow(all_cov) > 2){
        z_effect_size <- as.matrix(dfa$Estimates$D[,(nrow(seas)+1):ncol(dfa$Estimates$D)])
    } else {
        z_effect_size <- dfa$Estimates$D
    }

    nstream <- nrow(z_effect_size)
    ncov <- ncol(z_effect_size)

    #convert effects back to original scale (units response/units covar)
    #x/sd(response) = 1/sd(covar); x is the effect size scale factor
    rescaled_effect_size <- matrix(NA, nrow=nstream, ncol=ncov)
    for(i in 1:nstream){
        sd_response <- trans$sds[i]
        for(j in 1:ncov){
            sd_covar <- sd(covs[,j], na.rm=TRUE)
            rescaled_effect_size[,j] <- z_effect_size[,j] * (sd_response/sd_covar)
        }
    }

    return(list(z=z_effect_size, rescaled=rescaled_effect_size))
}
out <- modified_eff_rescaler(mod_covs, NULL)

# saveRDS(out, 'oxygen_raw.rds')
# saveRDS(trans$lambdas, 'oxygen_lambdas.rds')

#read in objects. no need to run the above unless testing beyond the follwoing conditions
#these response variables were chosen because they're close to being normally distributed:
# boxcoxables() #run this to see for yourself

# press_lambdas <- readRDS('press_lambdas.rds')
# press_raw <- readRDS('press_raw.rds')
# press_trans <- readRDS('press_trans.rds')
cond_lambdas <- readRDS('cond_lambdas.rds')
cond_raw <- readRDS('cond_raw.rds')
cond_trans <- readRDS('cond_trans.rds')
temp_lambdas <- readRDS('temp_lambdas.rds')
temp_raw <- readRDS('temp_raw.rds')
temp_trans <- readRDS('temp_trans.rds')
oxygen_lambdas <- readRDS('oxygen_lambdas.rds')
oxygen_raw <- readRDS('oxygen_raw.rds')
oxygen_trans <- readRDS('oxygen_trans.rds')

#temp comparisons (z-scored)
plot(temp_lambdas, temp_trans[[1]])
plot(temp_lambdas, temp_raw[[1]], pch=20, cex=2)
plot(temp_lambdas, temp_trans[[1]]-temp_raw[[1]])
#(original scale)
plot(temp_lambdas, temp_trans[[2]])
plot(temp_lambdas, temp_raw[[2]], pch=20, cex=2)
plot(temp_lambdas, temp_trans[[2]]-temp_raw[[2]])

#cond comparisons (z-scored)
plot(cond_lambdas, cond_trans[[1]])
plot(cond_lambdas, cond_raw[[1]], pch=20, cex=2)
plot(cond_lambdas, cond_trans[[1]]-cond_raw[[1]])
#(original scale)
plot(cond_lambdas, cond_trans[[2]])
plot(cond_lambdas, cond_raw[[2]], pch=20, cex=2)
plot(cond_lambdas, cond_trans[[2]]-cond_raw[[2]])

#oxygen comparisons (z-scored)
plot(oxygen_lambdas, oxygen_trans[[1]])
plot(oxygen_lambdas, oxygen_raw[[1]], pch=20, cex=2)
plot(oxygen_lambdas, oxygen_trans[[1]]-oxygen_raw[[1]])
#(original scale)
plot(oxygen_lambdas, oxygen_trans[[2]])
plot(oxygen_lambdas, oxygen_raw[[2]], pch=20, cex=2)
plot(oxygen_lambdas, oxygen_trans[[2]]-oxygen_raw[[2]])
