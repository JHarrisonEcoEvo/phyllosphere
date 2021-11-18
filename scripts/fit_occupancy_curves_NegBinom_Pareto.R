rm(list=ls())

##############################
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG_OCCUPANCY.csv",
                stringsAsFactors = F)
dat$sample <- gsub("X","", dat$sample)

meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", 
                 stringsAsFactors = F)

dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG_OCCUPANCY.csv",
                stringsAsFactors = F)
dat$sample <- gsub("X","", dat$sample)

meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv", 
                 stringsAsFactors = F)

merged_dat <- merge(dat, meta, by.x = "sample", by.y = "sample")

#combine rows by host, rebinarize, and redetermine prevalence

#First aggregate by taxon and compartment for a given microbe, thus combining info
#from multiple replicates. 
counts <- merged_dat[,grep("Zotu",names(merged_dat))]

agg_host <- data.frame(matrix(nrow = length(unique(merged_dat$taxon_final))*2, #bc of two compartments
                                ncol = 1))
k <- 1
for(i in 1:length(counts)){
  agg_host[,k] <- aggregate(counts[,i] ~ 
                              merged_dat$taxon_final + merged_dat$compartment, 
                              FUN = sum)[,3]
  colnames(agg_host)[k] <- names(counts)[i]
  k <- k + 1
}
#taking the first couple columns from our call to aggregate for labeling purposes and
#combining with our dataframe.
agg_host <- data.frame(aggregate(counts[,i] ~ 
                                   merged_dat$taxon_final + merged_dat$compartment, 
                                 FUN = sum)[,1:2],
                       agg_host)

#rebinarize our dataframe.
agg_host[,3:length(agg_host)] <- ifelse(agg_host[,3:length(agg_host)] > 0, 1, 0)
#

#add grwoth habit
agg_host$habit <- NA
for(i in 1:length(agg_host[,1])){
  agg_host$habit[i] <- unique(meta$habit[meta$taxon_final ==
                                    agg_host$merged_dat.taxon_final[i]])
}

#TEST CODE####################################
#COMBINE all the taxa so that we have one vector
#that conveys the number of hosts that each microbe was present within. 

#IMPORTANT: choose some subset of the data here
# hostUse <- colSums(agg_host[agg_host$merged_dat.compartment == "EN",
#                             grep("Zotu", names(agg_host))])
# length(hostUse)
#hostUse <- hostUse[hostUse != 0]

# z <- rpareto(1000, scale = pareto.MLE(c(hostUse))[1], 
#               shape = pareto.MLE(c(hostUse))[2]) 
# # plot((density(z)))
#  hist(z)
 
#Alternatively we can use ptsuite if we wanted to use non MLE estimators for some reason.
#install.packages("ptsuite")
# library(ptsuite)
# generate_all_estimates(hostUse)
# 
# #a function to make a qq plot to see if the Pareto is an ok fit
# pareto_qq_test(hostUse)
# 
# #Another way to test if the Pareto is a good fit. 
# pareto_test(hostUse)
# 
# library(EnvStats)
# epareto(hostUse, method = "mle")


#TRY negative binomial
# library(fitdistrplus)    # fits distributions using maximum likelihood
# library(gamlss)          # defines pdf, cdf of ZIP
# library(MASS)
# library(actuar)
# 
# fit_zipp <- fitdist(hostUse,
#                     'pareto',
#                     start=list(shape = 0.5, scale = 1))
# #plot(fit_zipp)
# out_stats <- gofstat(fit_zipp)
# out_stats

# 
# fit_zip2 <- fitdist(hostUse, 
#                     'nbinom',
#                     start = list(mu = 5, size = 0.1)) 

# Zero-inflated Poisson fails for some reason
#fit_zip = fitdist(i.vec, 'ZIP', start = list(mu = 2, sigma = 0.5))

# VISUALIZE TEST AND COMPUTE GOODNESS OF FIT    
#plot(fit_zip2)
# out_statsnb <- gofstat(fit_zip2)
# out_statsnb

#######################
#######################
#It seems that the Pareto has lower AIC, 
#despite having more parameters. But it requires zeros to be omitted.
#######################
#######################

############
#Get to modeling


nbfitter <- function(x){
  fit <- fitdist(x,
                    'nbinom',
                    method="mle")
#can cconvert from this parameterization to the standard version via this: 
#μ = n(1-p)/p 
#see Linden et al. 2011

#prob = size/(size+mu)

#can sum log probabilities to compare models if desired
ll <- sum(log(
  pnbinom(x, 
            # size = length(hostUseT),
          size = fit$estimate[names(fit$estimate) == "size"], 
             mu = fit$estimate[names(fit$estimate) == "mu"]
          )))
return(list(ll=ll,
             gof = gofstat(fit),
             fit = fit))
}

hostUseENTREE <- colSums(agg_host[agg_host$merged_dat.compartment == "EN" &
                                    agg_host$habit == "tree",
                                  grep("Zotu", names(agg_host))])

outENTREE <- nbfitter(hostUseENTREE)

hostUseEPTREE <- colSums(agg_host[agg_host$merged_dat.compartment == "EP" &
                                    agg_host$habit == "tree",
                                  grep("Zotu", names(agg_host))])

outEPTREE <-nbfitter(hostUseEPTREE)

#QC and interesting thing
# sum(log(
  # pnbinom(hostUseENTREE,
  #         size = outEPTREE$fit$estimate[names(outEPTREE$fit$estimate) == "size"],
  #         mu = outEPTREE$fit$estimate[names(outEPTREE$fit$estimate) == "mu"]
  # )))
# #Using the EN parameters leads to a lower log likelihood, as expected
# sum(log(
#   pnbinom(hostUseENTREE,
#           size = outENTREE$fit$estimate[names(outENTREE$fit$estimate) == "size"], 
#           mu = outENTREE$fit$estimate[names(outENTREE$fit$estimate) == "mu"]
#   )))



hostUseENForb<- colSums(agg_host[agg_host$merged_dat.compartment == "EN" &
                                    agg_host$habit == "forb",
                                  grep("Zotu", names(agg_host))])

nbfitter(hostUseENForb)

table(
  pnbinom(hostUseENForb,
                  size = outEPTREE$fit$estimate[names(outEPTREE$fit$estimate) == "size"], 
                  mu = outEPTREE$fit$estimate[names(outEPTREE$fit$estimate) == "mu"])
  < 0.25)
        
        
hostUseEN<- colSums(agg_host[agg_host$merged_dat.compartment == "EN" ,
                                 grep("Zotu", names(agg_host))])

outnben <- nbfitter(hostUseEN)

hostUseEP<- colSums(agg_host[agg_host$merged_dat.compartment == "EP" ,
                             grep("Zotu", names(agg_host))])

nbfitter(hostUseEP)

#Sort of strange way to see how a model does on another data set. 
t1 <- table(
  pnbinom(hostUseEP,
          size = outnben$fit$estimate[names(outnben$fit$estimate) == "size"], 
          mu = outnben$fit$estimate[names(outnben$fit$estimate) == "mu"])
  < 0.15)

t1[2] / sum(t1) 

t2 <- table(
  pnbinom(hostUseEN,
          size = outnben$fit$estimate[names(outnben$fit$estimate) == "size"], 
          mu = outnben$fit$estimate[names(outnben$fit$estimate) == "mu"])
  < 0.15)

t2[2] / sum(t1) 
