rm(list=ls())
load("processedData/ITS_div.Rdata")

diversity_plotter <- function(div, color_vec = "black"){
  #Extract the largest density value to make plot dimensions
  densities <- NA
  for(i in 1:length(div$entropy_pi)){
    densities[i] <- max(density(div$entropy_pi[[i]])$y)
  }
  
  plot(NULL, 
       xlim = c(0, 1.1*max(unlist(div$entropy_pi))),
       ylim = c(0, 1.1*max(densities)),
       las = 2, 
       ylab = "density",
       xlab = "entropy",
       frame.plot = F)
  if(length(color_vec) !=length(div$entropy_pi)){
    print("WARNING: The color vector that was passed in does not have the same number of elements as there were sampling groups.")
    print("CNVRG will use red for plotting. Sorry if this isn't what you wanted!")

      for(i in 1:length(div$entropy_pi)){
        lines(density(div$entropy_pi[[i]]), col = "red")
      }
    }
  
  for(i in 1:length(div$entropy_pi)){
    lines(density(div$entropy_pi[[i]]), col = color_vec[i])
  }
}

library(viridis)
color_vec <- viridis_pal(option = "viridis")(250)  

Need to update meta so it is in the correct order
meta <- read.csv("./processedData/treatments_metadata.csv", stringsAsFactors = F)
color_vec[meta$lifehistory == "blue"]
diversity_plotter(div, color_vec = color_vec)
