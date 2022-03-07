setwd('~/Dropbox/Phoebe_Files_Folder/hierarchical_modeling')

library(ggplot2)
library(parallel)

####These dataframes will be converted into the output csv's holding the results of the models. 

#meaningful effect: % of draws the survival rate for a given treatment is higher than the mean across all the outcomes 
#meaningful effect lower: % of draws the survival rate for a given treatment is lower than the mean across all the outcomes 
#meanignful_effect_all: % of draws the survival rate for a given treatment is higher than all other outcomes 
 
meaningful_effect <- data.frame('CCSR',"A_Thiazide_Diuretic", "A_ARB",   "A_CCB", "A_Beta_Blocker","A_ACE")
colnames(meaningful_effect) <-c('CCSR', "A_Thiazide_Diuretic", "A_ARB", "A_CCB", "A_Beta_Blocker","A_ACE")

meaningful_effect_lower <- data.frame('CCSR',"A_Thiazide_Diuretic", "A_ARB",   "A_CCB", "A_Beta_Blocker","A_ACE")
colnames(meaningful_effect_lower)  <-  c('CCSR', "A_Thiazide_Diuretic", "A_ARB", "A_CCB", "A_Beta_Blocker","A_ACE")


meaningful_effect_all <- data.frame('CCSR',"A_Thiazide_Diuretic", "A_ARB",   "A_CCB", "A_Beta_Blocker","A_ACE")
colnames(meaningful_effect_all)  <-  c('CCSR', "A_Thiazide_Diuretic", "A_ARB", "A_CCB", "A_Beta_Blocker","A_ACE")



#### This is where you define the directory that holds the survival information files. 
#### Each file contains the 100 bootstrap survival results at a given timepoint for the outcome specificed in the filename. 
#### Example File: BLD001_one_year_survival.csv : outcome CCSR BLD001, timepoint 365 days 

dir = paste(getwd(),'hypertensive_IPW_one_year_survival', sep='/')
print(dir)
list.files(dir)

skip = c('BLD005', 'CIR002', 'CIR007', 'CIR008', 'INF001','INF011','MAL001','MAL002','MAL003', 'MAL004','MAL005','MAL006', 'MAL007','MAL008', 'MAL009','MAL010', 'MUS008', 'MUS032', 'MUS035', 'MUS036', 'NEO008', 'NEO014', 'NEO016','NEO020', 'NEO021', 'NEO023','NEO036', 'NEO038', 'NEO049',  'NEO050', 'NEO060', 'NEO061', 'NEO064', 'NEO066', 'NEO069','NVS003',  'NVS005')


files = list.files(dir)
files = files[!substr(files,1,6) %in% skip]
ccsr_table = read.csv('~/Dropbox/Phoebe_Files_Folder/Hierarchical_Modeling/ccsr_to_ICD10.csv')
ccsr_table = unique(ccsr_table[,c("CSSR.Code","CCSR.Category.Description")])

for (i in 1:length(files)) { #loop through all files in directory to run model for each file in one go 
  filename = files[i]
  ccsr = sapply(strsplit(filename,"_"), `[`, 1) #get the CCSR (BLD001) from the filename BLD001_one_year_survival.csv

  file_dir = paste('hypertensive_IPW_one_year_survival',filename, sep='/') 
  year_survival <- read.csv(file_dir)
  year_survival$Unnamed..0 <- NULL

  year_survival  <- year_survival[c("A_Thiazide_Diuretic", "A_ARB", "A_CCB", "A_Beta_Blocker","A_ACE")] #reorder columns so treatments in same order each time model is run 

  year_effects <- 1-t(colMeans(year_survival,na.rm = TRUE))  #1x5 dataframe where each column now contains the mean survival from the 100 boostrapps 

  treatments <- colnames(year_survival) #list of column names 1x5 

  cov_mat <- solve(cov(year_survival,use = "complete.obs")) #5x5 covariance matrix where the diagonal is the variance of a single treatment 
  
  
  #inputs of the stan model 
  dataList <- list( 
    I = 5,
    all_cov = cov_mat,
    delta_est = year_effects[1,]
  )

  modelString ="model{
    for (i in 1:I) {
        delta_real[i] ~ dnorm(mu,tau)
    }
    delta_est ~ dmnorm(delta_real,all_cov)
    mu ~ dunif(0,.3)
    sigma ~ dunif(0,.1)
    tau<-pow(sigma,-2)
  }
  "
writeLines( modelString , con="TEMPmodel.txt")

  # Run the chains:
  jagsModel =jags.model( file="TEMPmodel.txt", data=dataList , #inits=initsList ,
                       n.chains=3, n.adapt=1000)     
  update( jagsModel , n.iter=2000)
  codaSamples =coda.samples( jagsModel , variable.names=c("mu","delta_real","sigma") ,n.iter=2000,thin = 5)
  samps =data.frame(do.call('rbind',codaSamples))
  names(samps)[1:5] = c("Thiazide_Diuretic", "ARB", "CCB", "Beta_Blocker","ACE")
  ##Creates Interval plot of resulting confidence intervals from model 
  save(samps,file=paste0('one_outcome_posterior_samples/one_outcome_posterior_draws_',ccsr,'.RData'))
  jpeg(file=paste0('one_outcome_posterior_plots/one_outcome_posterior_plot_',ccsr,'.jpg'))
  posterior_plot <- mcmc_intervals(samps,pars=names(samps)[1:5],prob_outer = .95)+labs(title =paste0(ccsr_table$CCSR.Category.Description[ccsr_table$CCSR.Code==ccsr],", Single Outcome Pooling"),
                                                    x="Cumulative Incidence Rate at 1 Year After Baseline")
  print(posterior_plot)
  dev.off()
}


effects_1out = data.frame(matrix(NA,nrow=length(files),ncol=9))
names(effects_1out) = c('ccsr',"Thiazide_Diuretic", "ARB", "CCB", "Beta_Blocker","ACE",'min','min_RR','min_prob')

for (i in 1:length(files)) { #loop through all files in directory to run model for each file in one go 
  filename = files[i]
  ccsr = sapply(strsplit(filename,"_"), `[`, 1) #get the CCSR (BLD001) from the filename BLD001_one_year_survival.csv
  effects_1out$ccsr[i] = sapply(strsplit(filename,"_"), `[`, 1) #get the CCSR (BLD001) from the filename BLD001_one_year_survival.csv
  load(paste0('one_outcome_posterior_samples/one_outcome_posterior_draws_',ccsr,'.RData'))
  year_effects = apply(samps[,1:5],2,mean)
  effects_1out[i,2:6] = year_effects
  ord = order(year_effects)
  effects_1out$min[i] = names(effects_1out[ord[1]+1])
  effects_1out$min_RR[i] = mean(samps[,ord[1]]/samps[,ord[2]])
  effects_1out$min_prob[i] = mean(samps[,ord[1]]==apply(samps[,1:5],1,min)) #pnorm(mean(year_survival[,ord[1]]-year_survival[,ord[2]],na.rm=T)/sd(year_survival[,ord[1]]-year_survival[,ord[2]],na.rm=T))
  effects_1out$max[i] = names(effects[ord[5]+1])
  effects_1out$max_RR[i] = mean(samps[,ord[5]]/samps[,ord[4]],na.rm=T)
  effects_1out$max_prob[i] = mean(samps[,ord[5]]==apply(samps[,1:5],1,max))
  # second_samps = apply(samps[,1:5],1,function(x) x[order(x)[2]])
  # for(j in 1:5){
  #   prob_best_1out[i,j] = mean(samps[,j] / second_samps <= .8)
  # }
}

names(effects_1out)[1] = 'CCSR.Code'
effects_1out = merge(effects_1out,ccsr_table)
View(effects_1out[effects_1out$min_prob>=.9 & effects_1out$min_RR<=.8,])
View(effects_1out[effects_1out$max_prob>=.9 & effects_1out$max_RR>=1.2,])


ccsrs = substr(files,1,6)
prob_best_1out$CCSR.Code = ccsrs  
prob_best_1out = merge(prob_best_1out,ccsr_table) 
keepers = which(apply(prob_best_1out[,2:6],1,max)>=.5)
View(prob_best_1out[keepers,])

