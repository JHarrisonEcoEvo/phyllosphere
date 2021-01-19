#Calculate richness and diversity by host. Do an ordination by host.
#Samples are lumped within a site because the comparison of interest is site and host
#and OTU clustering threshold. 

#We will compare shifts in a site across OTU thresholds, so sample number
#will not confound the analysis.

dat <- read.table("data/its90_otuTable", header = T, stringsAsFactors = F)


#First, check that the dimensions of the table show the same number of samples 
#for each dataset.
if(dim(dat)[1] != 3857){
  print("Warning: dimensions of data do not appear correct")
}


  
#Calculate richness
#Output richness as a table

#Calculate diversity
#Output diversity as a table


#Output ordination of data

