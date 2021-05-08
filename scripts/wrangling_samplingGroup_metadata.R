#This script condenses the metadata by sampmling group. 
#this makes it easier to assign colors /subset the data by group during 
#other analyses

rm(list=ls())
dat <- read.csv("processedData/treatments_metadata.csv", stringsAsFactors = F)
meta <- read.csv("processedData/metadata.csv", stringsAsFactors = F)
meta <- meta[meta$substrate == "plant",]
#add compartment column
dat$compartment <- "EP"
dat$compartment[grep("EN", dat$x)] <- "EN"
table(dat$compartment)

dat$plant <- gsub("(\\d+_\\d+_\\d+)E.*", "\\1", dat$x)

dat$taxon <- NA
dat$location <- NA
dat$lifehistory <- NA
for(i in unique(dat$plant)){
  if(length(unique(meta$taxon.y[grep(i, meta$region_site_plant)])) > 1){
    print(i)
  }
  dat$taxon[dat$plant == i] <- 
    unique(meta$taxon.y[grep(i, meta$region_site_plant)])
  dat$location[dat$plant == i] <- 
    unique(meta$region_site[grep(i, meta$region_site_plant)])
  dat$lifehistory[dat$plant == i] <- 
    unique(meta$habit[grep(i, meta$region_site_plant)])
}

write.csv(dat, file = "./processedData/treatments_metadata.csv")
