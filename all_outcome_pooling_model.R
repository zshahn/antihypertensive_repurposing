setwd('~/Dropbox/Phoebe_Files_Folder/hierarchical_modeling')
library(devtools)
library(parallel)
library(stringr)
#files = c('BLD001_year_survival_100.csv', 'CIR001_year_survival_100.csv')

dir = paste(getwd(),'hypertensive_IPW_one_year_survival', sep='/')
#dir = paste(getwd(),'one_year_survival', sep='/')
print(dir)
list.files(dir)
print(length(list.files(dir)))
skip = c('BLD005', 'CIR002', 'CIR007', 'CIR008', 'INF001','INF011','MAL001','MAL002','MAL003', 'MAL004','MAL005','MAL006', 'MAL007','MAL008', 'MAL009','MAL010', 'MUS008', 'MUS032', 'MUS035', 'MUS036', 'NEO008', 'NEO014', 'NEO016','NEO020', 'NEO021', 'NEO023','NEO036', 'NEO038', 'NEO049',  'NEO050', 'NEO060', 'NEO061', 'NEO064', 'NEO066', 'NEO069','NVS003',  'NVS005')

for (filename in list.files(dir)){
  #print(substr(filename,1,6))
  if (substr(filename,1,6) %in% skip) {
    print(filename)
  }
}

year_list = list()
cov_list = list()

files = list.files(dir)
files = files[!substr(files,1,6) %in% skip]
for (i in 1:length(files)){ 
  filename = files[i]
  outcome_code = sapply(strsplit(filename,"_"), `[`, 1)
  #assign(outcome_code,read.csv(year_survival_file))
  
  year_surv <- read.csv(paste('hypertensive_IPW_one_year_survival',filename, sep='/'))
  year_surv$Unnamed..0 <- NULL
  year_surv  <- 1-year_surv[c("A_Thiazide_Diuretic", "A_ARB", "A_CCB", "A_Beta_Blocker","A_ACE")]
  
  cov_mat <- cov(year_surv,use = "complete.obs")
  
  year_surv_average <- t(colMeans(year_surv,na.rm = TRUE))
  year_list[[i]] <- year_surv_average[1,]
  cov_list[[i]] <- cov_mat
}

N = length(files) #length(length(list.files(dir)))  #number of outcomes 
I = 5  #number of treatments

cov_matrix <- array(NA_real_, dim = c(N,I,I))
for(n in 1:N) cov_matrix[n,,] <- cov_list[[n]]
dim(cov_matrix)


year_array <- array(NA_real_, dim = c(N,I))
for(n in 1:N) {
  year_array[n,]  <- year_list[[n]]
}

dim(year_array)


print(cov_matrix)

cov_matrix <- array(NA_real_, dim = c(N,I,I))
for(n in 1:N) cov_matrix[n,,] <- solve(cov_list[[n]])
dim(cov_matrix)

dataList =list(N=N, I=I,delta_est =year_array,
               all_cov=cov_matrix)

#log RR pooling model
# modelString ="model{
# for (i in 1:N) {
#    for (j in 1:I){
#     alpha[i,j] ~ dnorm(0,tau)
#     delta_real[i,j] <- mu[i]*exp(alpha[i,j])
#    }
#   }
# for(i in 1:N){
#   mu[i] ~ dunif(0,.3)
#   delta_est[i,] ~ dmnorm(delta_real[i,],all_cov[i,,])
# }
# sigma ~ dunif(0,.25)
# tau<-pow(sigma,-2)
# }
# "

#RR pooling model
modelString ="model{
for (i in 1:N) {
   for (j in 1:I){
    alpha[i,j] ~ dnorm(1,tau)
    delta_real[i,j] <- mu[i]*alpha[i,j]
   }
  }
for(i in 1:N){
  mu[i] ~ dunif(0,.3)
  delta_est[i,] ~ dmnorm(delta_real[i,],all_cov[i,,])
}
sigma ~ dunif(0,.25)
tau<-pow(sigma,-2)
}
"

writeLines( modelString , con="TEMPmodel.txt")

library(rjags)
# Run the chains:
jagsModel =jags.model( file="TEMPmodel.txt", data=dataList , #inits=initsList ,
                       n.chains=3, n.adapt=1000)     
update( jagsModel , n.iter=2000)
codaSamples2 =coda.samples( jagsModel , variable.names=c("mu","delta_real","sigma","alpha") ,n.iter=2000,thin = 5)
s2 =data.frame(do.call('rbind',codaSamples2))
hist(s2$sigma)
save(s2,file='multi_outcome_hypertensive_RR_1.RData')
load('multi_outcome_hypertensive_RR_1.RData')

post_probs2 = matrix(NA,nrow = N,ncol=I)
for(i in 1:N){
  post_delta_i = s2[,grep(paste0("delta_real\\.",i,"\\."),names(s2))]
  seconds = apply(post_delta_i,1,function(x)x[order(x)[2]])
  for(j in 1:I){
    # post_probs2[i,j] = mean(post_delta_i[,j] == apply(post_delta_i,1,min))
    post_probs2[i,j] = mean(post_delta_i[,j] <= .8*seconds)
  }
}
colnames(post_probs2) = c("A_Thiazide_Diuretic", "A_ARB", "A_CCB", "A_Beta_Blocker","A_ACE")
post_probs2 = data.frame(post_probs2)
post_probs2$CCSR.Code = substr(files,1,6)

ccsr_table = read.csv('~/Dropbox/Phoebe_Files_Folder/Hierarchical_Modeling/ccsr_to_ICD10.csv')
ccsr_table = unique(ccsr_table[,c("CSSR.Code","CCSR.Category.Description")])

names(ccsr_table)[1] = 'CCSR.Code'
post_probs2 = merge(post_probs2,ccsr_table)

#make forest plots
library(bayesplot)
library(ggplot2)

library(forestplot)
library(ggplot2)
library(dplyr)

thiaz_irr_means = rep(NA,N)
arb_irr_means = rep(NA,N)
ccb_irr_means = rep(NA,N)
beta_irr_means = rep(NA,N)
ace_irr_means = rep(NA,N)

thiaz_irr_lowers = rep(NA,N)
arb_irr_lowers = rep(NA,N)
ccb_irr_lowers = rep(NA,N)
beta_irr_lowers = rep(NA,N)
ace_irr_lowers = rep(NA,N)

thiaz_irr_uppers = rep(NA,N)
arb_irr_uppers = rep(NA,N)
ccb_irr_uppers = rep(NA,N)
beta_irr_uppers = rep(NA,N)
ace_irr_uppers = rep(NA,N)

for(i in 1:N){
  post_delta_i = s2[,grep(paste0("delta_real\\.",i,"\\."),names(s2))]
  thiaz = post_delta_i[,1]
  not_thiaz = apply(post_delta_i[,-1],1,mean)
  arb = post_delta_i[,2]
  not_arb = apply(post_delta_i[,-2],1,mean)
  ccb = post_delta_i[,3]
  not_ccb = apply(post_delta_i[,-3],1,mean)
  beta = post_delta_i[,4]
  not_beta = apply(post_delta_i[,-4],1,mean)
  ace = post_delta_i[,5]
  not_ace = apply(post_delta_i[,-5],1,mean)
  
  thiaz_irr_means[i] = mean(thiaz/not_thiaz)
  arb_irr_means[i] = mean(arb/not_arb)
  ccb_irr_means[i] = mean(ccb/not_ccb)
  beta_irr_means[i] = mean(beta/not_beta)
  ace_irr_means[i] = mean(ace/not_ace)
  
  thiaz_irr_lowers[i] = quantile(thiaz/not_thiaz,.025)
  arb_irr_lowers[i] = quantile(arb/not_arb,.025)
  ccb_irr_lowers[i] = quantile(ccb/not_ccb,.025)
  beta_irr_lowers[i] = quantile(beta/not_beta,.025)
  ace_irr_lowers[i] = quantile(ace/not_ace,.025)
  
  thiaz_irr_uppers[i] = quantile(thiaz/not_thiaz,.975)
  arb_irr_uppers[i] = quantile(arb/not_arb,.975)
  ccb_irr_uppers[i] = quantile(ccb/not_ccb,.975)
  beta_irr_uppers[i] = quantile(beta/not_beta,.975)
  ace_irr_uppers[i] = quantile(ace/not_ace,.975)
}
  
thiaz_irr_summary = cbind(mean = thiaz_irr_means, lower=thiaz_irr_lowers, upper=thiaz_irr_uppers)
arb_irr_summary = cbind(mean = arb_irr_means, lower=arb_irr_lowers, upper=arb_irr_uppers)
ccb_irr_summary = cbind(mean = ccb_irr_means, lower=ccb_irr_lowers, upper=ccb_irr_uppers)
beta_irr_summary = cbind(mean = beta_irr_means, lower=beta_irr_lowers, upper=beta_irr_uppers)
ace_irr_summary = cbind(mean = ace_irr_means, lower=ace_irr_lowers, upper=ace_irr_uppers)


png(filename='forest_thiazide_irr.png', height=4000, width=1000, pointsize=12)
thiaz_irr_summary %>% 
  forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
             labeltext=data.frame(post_probs2$CCSR.Category.Description),
            title="Thiazide IRR 95% Credible Intervals")
dev.off()

png(filename='forest_arb_irr.png', height=4000, width=1000, pointsize=12)
arb_irr_summary %>% 
  forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
             labeltext=data.frame(post_probs2$CCSR.Category.Description),
             title="ARB IRR 95% Credible Intervals")
dev.off()

png(filename='forest_ccb_irr.png', height=4000, width=1000, pointsize=12)
ccb_irr_summary %>% 
  forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
             labeltext=data.frame(post_probs2$CCSR.Category.Description),
             title="CCB IRR 95% Credible Intervals")
dev.off()

png(filename='forest_beta_irr.png', height=4000, width=1000, pointsize=12)
beta_irr_summary %>% 
  forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
             labeltext=data.frame(post_probs2$CCSR.Category.Description),
             title="Beta Blocker IRR 95% Credible Intervals")
dev.off()

png(filename='forest_ace_irr.png', height=4000, width=1000, pointsize=12)
ace_irr_summary %>% 
  forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
             labeltext=data.frame(post_probs2$CCSR.Category.Description),
             title="ACE Inhibitor IRR 95% Credible Intervals")
dev.off()



# arb_alpha_means = apply(s2[,grep('alpha\\..*\\.2',names(s2))],2,mean)
# arb_alpha_lowers = apply(s2[,grep('alpha\\..*\\.2',names(s2))],2,function(x)quantile(x,.025))
# arb_alpha_uppers = apply(s2[,grep('alpha\\..*\\.2',names(s2))],2,function(x)quantile(x,.975))
# arb_alpha_summary = cbind(mean = arb_alpha_means, lower=arb_alpha_lowers, upper=arb_alpha_uppers)
# 
# png(filename='forest_arb.png', height=4000, width=1000, pointsize=12)
# arb_alpha_summary %>% 
#   forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
#              labeltext=data.frame(post_probs2$CCSR.Category.Description),
#              title=expression(paste("ARB ",alpha," 95% Credible Intervals")))
# dev.off()
# 
# ccb_alpha_means = apply(s2[,grep('alpha\\..*\\.3',names(s2))],2,mean)
# ccb_alpha_lowers = apply(s2[,grep('alpha\\..*\\.3',names(s2))],2,function(x)quantile(x,.025))
# ccb_alpha_uppers = apply(s2[,grep('alpha\\..*\\.3',names(s2))],2,function(x)quantile(x,.975))
# ccb_alpha_summary = cbind(mean = ccb_alpha_means, lower=ccb_alpha_lowers, upper=ccb_alpha_uppers)
# 
# png(filename='forest_ccb.png', height=4000, width=1000, pointsize=12)
# ccb_alpha_summary %>% 
#   forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
#              labeltext=data.frame(post_probs2$CCSR.Category.Description),
#              title=expression(paste("CCB ",alpha," 95% Credible Intervals")))
# dev.off()
# 
# beta_alpha_means = apply(s2[,grep('alpha\\..*\\.4',names(s2))],2,mean)
# beta_alpha_lowers = apply(s2[,grep('alpha\\..*\\.4',names(s2))],2,function(x)quantile(x,.025))
# beta_alpha_uppers = apply(s2[,grep('alpha\\..*\\.4',names(s2))],2,function(x)quantile(x,.975))
# beta_alpha_summary = cbind(mean = beta_alpha_means, lower=beta_alpha_lowers, upper=beta_alpha_uppers)
# 
# png(filename='forest_beta.png', height=4000, width=1000, pointsize=12)
# beta_alpha_summary %>% 
#   forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
#              labeltext=data.frame(post_probs2$CCSR.Category.Description),
#              title=expression(paste("Beta Blocker ",alpha," 95% Credible Intervals")))
# dev.off()
# 
# ace_alpha_means = apply(s2[,grep('alpha\\..*\\.5',names(s2))],2,mean)
# ace_alpha_lowers = apply(s2[,grep('alpha\\..*\\.5',names(s2))],2,function(x)quantile(x,.025))
# ace_alpha_uppers = apply(s2[,grep('alpha\\..*\\.5',names(s2))],2,function(x)quantile(x,.975))
# ace_alpha_summary = cbind(mean = ace_alpha_means, lower=ace_alpha_lowers, upper=ace_alpha_uppers)
# 
# png(filename='forest_ace.png', height=4000, width=1000, pointsize=12)
# ace_alpha_summary %>% 
#   forestplot(clip=c(.5,1.5), xticks = c(.5, .75, 1,1.25,1.5),zero=1,
#              labeltext=data.frame(post_probs2$CCSR.Category.Description),
#              title=expression(paste("ACE Inhibitor ",alpha," 95% Credible Intervals")))
# dev.off()

#make posterior probability tables to rank repurposing opportunities
mcmc_intervals(s2,pars=c('alpha.1.1.','alpha.1.2.','alpha.1.3.','alpha.1.4.','alpha.1.5.'))
treatments = c("Thiazide_Diuretic", "ARB", "CCB", "Beta_Blocker","ACE")
for(i in 1:N){
  for(j in 1:I){
    # s2[,paste0('surv_prob.',i,'.',j,'.')] = exp(s2[,paste0('alpha.',i,'.',j,'.')])*s2[,paste0('mu.',i,'.')]
    s2[,paste0('surv_prob_',ccsrs[i],'_',treatments[j])] = s2[,paste0('alpha.',i,'.',j,'.')]*s2[,paste0('mu.',i,'.')]
    }
}

write.csv(s2,file='posterior_samples_full_pooling.csv',row.names = F)

for (i in 1:length(files)) { #loop through all files in directory to run model for each file in one go 
  filename = files[i]
  ccsr = sapply(strsplit(filename,"_"), `[`, 1) #get the CCSR (BLD001) from the filename BLD001_one_year_survival.csv
  temp_samps = s2[,grep(paste0('surv_prob_',ccsrs[i]),names(s2))]
  names(temp_samps) = c("Thiazide_Diuretic", "ARB", "CCB", "Beta_Blocker","ACE")
  jpeg(file=paste0('all_outcome_posterior_plots/all_outcome_posterior_plot_',ccsr,'_rr1.jpg'))
  post_plot = mcmc_intervals(temp_samps,pars=names(temp_samps),prob_outer = .95) + labs(title =paste0(ccsr_table$CCSR.Category.Description[ccsr_table$CCSR.Code==ccsr],", Full Pooling"),
                                                           x="Cumulative Incidence Rate at 1 Year After Baseline")
  print(post_plot)
  dev.off()
}

effects_all_out = data.frame(matrix(NA,nrow=length(files),ncol=9))
names(effects_all_out) = c('ccsr',"Thiazide_Diuretic", "ARB", "CCB", "Beta_Blocker","ACE",'min','min_RR','min_prob')

for(i in 1:length(files)) { #loop through all files in directory to run model for each file in one go 
  filename = files[i]
  ccsr = sapply(strsplit(filename,"_"), `[`, 1) #get the CCSR (BLD001) from the filename BLD001_one_year_survival.csv
  samps = s2[,c(paste0('surv_prob.',i,'.1.'),paste0('surv_prob.',i,'.2.'),paste0('surv_prob.',i,'.3.'),paste0('surv_prob.',i,'.4.'),paste0('surv_prob.',i,'.5.'))]
  names(temp_samps) = c("Thiazide_Diuretic", "ARB", "CCB", "Beta_Blocker","ACE")
  effects_all_out$ccsr[i] = sapply(strsplit(filename,"_"), `[`, 1) #get the CCSR (BLD001) from the filename BLD001_one_year_survival.csv
  year_effects = apply(samps[,1:5],2,mean)
  effects_all_out[i,2:6] = year_effects
  ord = order(year_effects)
  effects_all_out$min[i] = names(effects_all_out[ord[1]+1])
  effects_all_out$min_RR[i] = mean(samps[,ord[1]]/samps[,ord[2]])
  effects_all_out$min_prob[i] = mean(samps[,ord[1]]==apply(samps[,1:5],1,min)) #pnorm(mean(year_survival[,ord[1]]-year_survival[,ord[2]],na.rm=T)/sd(year_survival[,ord[1]]-year_survival[,ord[2]],na.rm=T))
  effects_all_out$max[i] = names(effects[ord[5]+1])
  effects_all_out$max_RR[i] = mean(samps[,ord[5]]/samps[,ord[4]],na.rm=T)
  effects_all_out$max_prob[i] = mean(samps[,ord[5]]==apply(samps[,1:5],1,max))
  # second_samps = apply(samps[,1:5],1,function(x) x[order(x)[2]])
  # for(j in 1:5){
  #   prob_best_all_out[i,j] = mean(samps[,j] / second_samps <= .8)
  # }
}

names(effects_all_out)[1] = 'CCSR.Code'
effects_all_out = merge(effects_all_out,ccsr_table)
View(effects_all_out[effects_all_out$min_prob>=.9 & effects_all_out$min_RR<=.85,])
View(effects_all_out[effects_all_out$max_prob>=.9 & effects_all_out$max_RR>=1.2,])


