#plot rank abundance curves for fungal pis. Turns out this wasn't very interesting. 
rm(list=ls())

#Bring in data, packages, set up color scheme
load("processedData/ests.Rdata")
library(viridis)
color_vec <- viridis_pal(option = "viridis")(250)  
meta <- read.csv("./processedData/treatments_metadata.csv", stringsAsFactors = F)


dim(meta)
colvir <- viridis_pal(option = "viridis")(4)
color_vec[meta$lifehistory == "tree"] <- colvir[1]
color_vec[meta$lifehistory == "shrub"] <- colvir[2]
color_vec[meta$lifehistory == "forb"] <- colvir[3]
color_vec[meta$lifehistory == "graminoid"] <- colvir[4]

#Function from Mage, for making transparency easier
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#IMPORTANT the meta file is in the same order as the estimates. 
#See the CNVRG modeling and wrangling_samplingGroup scripts

finalindex <- (length(ests$pointEstimates_pi)-2)


toplot <- sort(unlist(ests$pointEstimates_pi[1,2:finalindex]))
plot(density(toplot))

for(i in 2:250){
  toplot <- sort(unlist(ests$pointEstimates_pi[i,2:finalindex]))
  lines(density(toplot))
}

