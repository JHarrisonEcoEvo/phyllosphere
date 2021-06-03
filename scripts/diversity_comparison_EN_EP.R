#J. Harrison.
#Univ. of Wyoming.
rm(list=ls())
load("./processedData/ITS_div_isd_normalized.Rdata")

str(div)

library(viridis)
color_vec <- viridis_pal(option = "viridis")(250)  

meta <- read.csv("./processedData/treatments_metadata.csv", stringsAsFactors = F)
dim(meta)
colvir <- viridis_pal(option = "viridis")(4)
color_vec[meta$lifehistory == "tree"] <- colvir[1]
color_vec[meta$lifehistory == "shrub"] <- colvir[2]
color_vec[meta$lifehistory == "forb"] <- colvir[3]
color_vec[meta$lifehistory == "graminoid"] <- colvir[4]

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

meta$div_means <- unlist(lapply(div$entropy_pi, FUN = mean))

comparison <- data.frame(matrix(ncol = 3))
k <- 1
for(taxon in unique(meta$taxon)){
  subsetMeta <- meta[meta$taxon == taxon,]
  for(j in unique(subsetMeta$location)){
    comparison[k,1] <- taxon
    comparison[k,2] <- j 
    doublesubset <- subsetMeta[subsetMeta$location == j, ]
    comparison[k,3] <- doublesubset$div_means[doublesubset$compartment == "EN"] -
      doublesubset$div_means[doublesubset$compartment == "EP"]
    k <- 1 + k
  }
}

#bring in elevation of site
sampleMeta <- read.csv("processedData/metadata.csv")

comparison$elevation <- NA
comparison$lifehistory <- NA

for(i in unique(meta$location)){
  comparison$elevation[comparison$X2 == i]  <- unique(sampleMeta$elev_m[sampleMeta$region_site == i])
  
}

for(i in unique(meta$taxon)){
  comparison$lifehistory[comparison$X1 == i]  <- unique(meta$lifehistory[meta$taxon == i])
}

library(lme4)
reg <- lmer(comparison$X3 ~ comparison$elevation + (1|comparison$X1))
summary(reg)
plot(comparison$X3, comparison$elevation)
summary(lm(comparison$X3 ~ comparison$elevation))
#Elevation not a good predictor of shift in diversity from EP to EN

comparison$lifehistory[79] <- "shrub"
reg <- aov(comparison$X3 ~ comparison$lifehistory)
summary(reg)
TukeyHSD(reg)
boxplot(comparison$X3 ~ comparison$lifehistory)
