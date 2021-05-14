rm(list=ls())
load("processedData/ITS_div.Rdata")

diversity_plotter <- function(div, color_vec = "black", xlab = "entropy"){
  #Extract the largest density value to make plot dimensions
  densities <- NA
  for(i in 1:length(div$entropy_pi)){
    densities[i] <- max(density(div$entropy_pi[[i]])$y)
  }
  
  plot(NULL, 
       xlim = c(0, 1.1*max(unlist(div$entropy_pi))),
       ylim = c(0,0.12),#c(0, 1.1*max(densities)),
       las = 2, 
       cex.lab = 1.2,
       ylab = "density",
       xaxt = "n",
       xlab = xlab,
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

div_means <- unlist(lapply(div$entropy_pi, FUN = mean))

############
############Plot

pdf(width = 7, height = 12, file = "./visuals/diversity_density_lifehistoryITS.pdf")
par(mfrow = c(2,1), mar = c(3, 4, 0, 0), oma = c(4, 1, 1,1))

diversity_plotter(div, color_vec = color_vec, xlab = "Shannon's equivalents")
axis(1, at = seq(0,4000, 1000), lab = seq(0,4000,1000))
legend(x = 0, y = 0.1, legend = c("tree", "shrub", "forb", "graminoid"),
       pch = 15, cex = 2, col = colvir, bty = "n")

text("fungi", x =2000, y = .125, xpd = NA, cex = 2)

stripchart(div_means ~ meta$compartment + meta$lifehistory,
           vertical = TRUE, 
       
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colvir[c(1,1,2,2,3,3,4,4)], alpha=0.7),
           cex=1,
           xaxt="n",
           yaxt="n",
           cex.lab = 1.2,
           ylab="Shannon's equivalents",
           xlim=c(0,8.5)
           ,ylim=c(1000,4500)
           ,frame.plot=F
)
boxplot(div_means ~ meta$compartment + meta$lifehistory, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        axes=F,
        col=c(add.alpha(colvir[1], alpha=0.3),
              add.alpha(colvir[1], alpha=0.6),
              add.alpha(colvir[2], alpha=0.3),
              add.alpha(colvir[2], alpha=0.6),
              add.alpha(colvir[3], alpha=0.3),
              add.alpha(colvir[3], alpha=0.6),
              add.alpha(colvir[4], alpha=0.3),
              add.alpha(colvir[4], alpha=0.6)),
        ylab=""
)

axis(side =2,
     at = seq(1000,4500,500),
     lab= seq(1000,4500,500),
     las = 2
)

axis(side =1,
     at = seq(1,8,1),
     lab= rep("", 8)
)
text(x = seq(1,8,1),
     y = 400,
     cex = 1.3,
     xpd = NA,
     srt = -50, c("EN forb","EP forb","EN gram.","EP gram.","EN shrub","EP shrub","EN tree","EP tree"))

dev.off()
pdf(width = 7, height = 12, file = "./visuals/diversity_density_lifehistory16S.pdf")
par(mfrow = c(2,1), mar = c(3, 4, 0, 0), oma = c(4, 1, 1,1))

#DO FOR BACTERIA
load("processedData/16s_div.Rdata")

div_means <- unlist(lapply(div$entropy_pi, FUN = mean))

par(mfrow = c(2,1), mar = c(3, 4, 0, 0), oma = c(4, 1, 1,1))

diversity_plotter(div, color_vec = color_vec, xlab = "Shannon's equivalents")
axis(1, at = c(seq(0,4000, 1000), 4500), lab =c(seq(0,4000, 1000), 4500))


text("bacteria", x =2000, y = .125, xpd = NA, cex = 2)

stripchart(div_means ~ meta$compartment + meta$lifehistory,
           vertical = TRUE, 
           
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colvir[c(1,1,2,2,3,3,4,4)], alpha=0.7),
           cex=1,
           xaxt="n",
           yaxt="n",
           cex.lab = 1.2,
           ylab="Shannon's equivalents",
           xlim=c(0,8.5)
           ,ylim=c(1000,4500)
           ,frame.plot=F
)
boxplot(div_means ~ meta$compartment + meta$lifehistory, 
        na.omit = T,
        add=T,
        las=2,
        outline=F,
        axes=F,
        col=c(add.alpha(colvir[1], alpha=0.3),
              add.alpha(colvir[1], alpha=0.6),
              add.alpha(colvir[2], alpha=0.3),
              add.alpha(colvir[2], alpha=0.6),
              add.alpha(colvir[3], alpha=0.3),
              add.alpha(colvir[3], alpha=0.6),
              add.alpha(colvir[4], alpha=0.3),
              add.alpha(colvir[4], alpha=0.6)),
        ylab=""
)

axis(side =2,
     at = seq(1000,4500,500),
     lab= seq(1000,4500,500),
     las = 2
)

axis(side =1,
     at = seq(1,8,1),
     lab= rep("", 8)
)
text(x = seq(1,8,1),
     y = 400,
     cex = 1.3,
     xpd = NA,
     srt = -50, c("EN forb","EP forb","EN gram.","EP gram.","EN shrub","EP shrub","EN tree","EP tree"))

dev.off()