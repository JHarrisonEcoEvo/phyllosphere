#Two part script. First part figures out taxa to model at landscape level (across all hosts and samples) based on prevalence
#Part two calculates prevalence among samples from only a SINGLE host and outputs combinations
#of host and microbe taxa to model. 

#########
# part 1#
#########
rm(list=ls())
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG", 
                fill = T, header = T, stringsAsFactors = F)

dat[,3:length(dat)] <- dat[,3:length(dat)] - 1

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good

#Figure out how many things were in 100 or more samples.
#This should be how many things I tried to model with the random forest
prev <- dat[,colSums(dat[,3:length(dat)]) >=100]

#Omit the isd, duds etc. 
prev <- prev[,grep("Zotu*", names(prev))]
dim(prev)

write.csv(names(prev), file = "./processedData/ITS_taxa_to_model_via_randomforest.csv", row.names = F)


############
# Bacteria
rm(list=ls())
dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG", 
                fill = T, header = T, stringsAsFactors = F)


dat[,3:length(dat)] <- dat[,3:length(dat)] - 1

#convert data to qualitative data, 1 if present 0 otherwise
dat[,3:length(dat)][dat[,3:length(dat)] > 0] <- 1
dim(dat) #still good

#Figure out how many things were in 100 or more samples.
#This should be how many things I tried to model with the random forest
prev <- dat[,colSums(dat[,3:length(dat)]) >=100]

#Omit the isd, duds etc. 
prev <- prev[,grep("Zotu*", names(prev))]
dim(prev)

write.csv(names(prev), file = "./processedData/sixteenS_taxa_to_model_via_randomforest.csv", row.names = F)

##########
# Part 2 #
###########

#Choose combinations of host and microbial taxa to model based on microbial prevalence and how much I sampled
#the host taxon. 

rm(list=ls())

#load data
possibles <- read.csv("./processedData/ITS_taxa_to_model_via_randomforest.csv", stringsAsFactors = F)
X <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", stringsAsFactors = F)
dat <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG", stringsAsFactors = F)

#Subset metadata to the most sampled hosts
tomodel <- names(which(table(X$taxon_final) > 30))
#Remove unknown stuff
tomodel <- tomodel[-grep("nknown", tomodel)]
X <- X[X$taxon_final %in% tomodel,]


##################################################
# wrangle the otu table as was done for modeling #
##################################################

#Note that we did considerable wrangling of these data using the
#"combine_pcr_dupes.*" script
dim(dat)

#Remove the one that was addedd

dat[,3:length(dat)] <- dat[,3:length(dat)] -1
dat$sample <- gsub("X","", dat$sample)
dat$region_site_plant <-  gsub("(\\d+_\\d+_\\d+)_.*","\\1", dat$sample)

# Loop through by host species and find all OTUs that are fairly prevalent. 
# Write the host taxon and the OTU to two fields in an output file to feed into subsequent modeling. 

keepers <- data.frame(matrix(ncol = 2, nrow = 1))
names(keepers) <- c("host","microbe")
k <- 1

for(i in unique(X$taxon_final)){
  #Find the region/site/plant prefixs for a given host
  #and extract the count data for those samples
  hits <- X$region_site_plant[X$taxon_final == i]
  workingdf <- dat[dat$region_site_plant %in% hits, ]
  
  #calculate prevalence for each taxon within this plant
  for(j in grep("Zotu", names(workingdf))){
    #figure out prevalence of the focal microbe in this plant
    table_microbe <- table(workingdf[,j] > 0)
      if("TRUE" %in% names(table_microbe)){
        #For microbes that are present in this taxon, extract the proportion of samples in which they occur
        proportion <- table_microbe["TRUE" == names(table_microbe)] / sum(table_microbe)
        #If microbes occur in 30% or more of samples, then write the combination of host and microbial taxon to the output
        if(proportion >= 0.3){
          keepers[k,] <- "NA" #clunky hack to avoid row number and data input amount discrepencies
          keepers$host[k] <- i
          keepers$microbe[k] <- names(workingdf)[j]
          k <- k + 1
        }
      }else{
        next
      }
    }
}
#sanity check. Lets look at one of the combos and see if it seems right.
#Pseudotsuga menziesii  Zotu1136
# 
# test <- X[X$taxon_final == "Pseudotsuga menziesii",]
# testdat <- dat[dat$region_site_plant %in% test$region_site_plant,]
# testdat[, grep("Zotu1136", names(testdat))] #pretty prevalent.

write.csv(keepers, file = "./processedData/combination_hosts_microbes_to_analyze_ITS.csv")
unique(keepers$microbe)

##############
# DO FOR 16S #
##############

rm(list=ls())

#load data
possibles <- read.csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv", stringsAsFactors = F)
X <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv", stringsAsFactors = F)
dat <- read.csv("./processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG", stringsAsFactors = F)

#Subset metadata to the most sampled hosts
tomodel <- names(which(table(X$taxon_final) > 30))
#Remove unknown stuff
tomodel <- tomodel[-grep("nknown", tomodel)]
X <- X[X$taxon_final %in% tomodel,]


##################################################
# wrangle the otu table as was done for modeling #
##################################################

#Note that we did considerable wrangling of these data using the
#"combine_pcr_dupes.*" script
dim(dat)

#Remove the one that was addedd

dat[,3:length(dat)] <- dat[,3:length(dat)] -1
dat$sample <- gsub("X","", dat$sample)
dat$region_site_plant <-  gsub("(\\d+_\\d+_\\d+)_.*","\\1", dat$sample)

# Loop through by host species and find all OTUs that are fairly prevalent. 
# Write the host taxon and the OTU to two fields in an output file to feed into subsequent modeling. 

keepers <- data.frame(matrix(ncol = 2, nrow = 1))
names(keepers) <- c("host","microbe")
k <- 1

for(i in unique(X$taxon_final)){
  #Find the region/site/plant prefixs for a given host
  #and extract the count data for those samples
  hits <- X$region_site_plant[X$taxon_final == i]
  workingdf <- dat[dat$region_site_plant %in% hits, ]
  
  #calculate prevalence for each taxon within this plant
  for(j in grep("Zotu", names(workingdf))){
    #figure out prevalence of the focal microbe in this plant
    table_microbe <- table(workingdf[,j] > 0)
    if("TRUE" %in% names(table_microbe)){
      #For microbes that are present in this taxon, extract the proportion of samples in which they occur
      proportion <- table_microbe["TRUE" == names(table_microbe)] / sum(table_microbe)
      #If microbes occur in 30% or more of samples, then write the combination of host and microbial taxon to the output
      if(proportion >= 0.3){
        keepers[k,] <- "NA" #clunky hack to avoid row number and data input amount discrepencies
        keepers$host[k] <- i
        keepers$microbe[k] <- names(workingdf)[j]
        k <- k + 1
      }
    }else{
      next
    }
  }
}

#sanity check. see ITS

write.csv(keepers, file = "./processedData/combination_hosts_microbes_to_analyze_16S.csv")
unique(keepers$microbe)



