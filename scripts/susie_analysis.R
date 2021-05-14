library(susieR)
library(Matrix) 
##Matrix is used to make matrix sparse for susieR
set.seed(666)
dat <- read.csv("processedData/16sp_estimates_wrangled_for_post_modeling_analysis.csv",
                stringsAsFactors = F)
meta <- read.csv("processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                 stringsAsFactors = F)

#Check which fields should be removed
tail(names(dat))
names(dat)[length(dat)-4]

#Will need to do this for each sampling group of interest!!!!
#Need to omit NAs in the metadata and make the meta and dat match afterward
dat_no_na <- as.matrix(dat[!is.na(meta$taxon.x), 5:(length(dat)-4)])
meta <- meta[!is.na(meta$taxon.x),]

res <- list()
taxa <- unique(meta$taxon.x)

for(i in unique(na.omit(taxa))){
  dummy <- rep(0, length(dat_no_na[,1]))
  dummy[meta$taxon.x == i] <- 1

  #Test run. Need to make sure spelling is correct for all ataxa DELTE THIS
  res$i <- susie(Matrix(dat_no_na, 
                         sparse = T),
             dummy,L=20)
}

hist(res$i$pip)

#Names should be in same order as pip.
names(dat)[5:(length(dat)-4)][res$i$pip > 0.05]

#Would need to then use these taxa in a predictive model to see how well
#they predict the host taxon, with k fold validation. 
#The question would become can specific taxa predict which host 
#
