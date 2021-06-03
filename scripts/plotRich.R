rm(list=ls())

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
load("./processedData/rich.Rdata")

rich_means0.001 <- apply(rich0.0001$rich_pi, MARGIN = 2, FUN=mean)
rich_means0.002 <- apply(rich0.0002$rich_pi, MARGIN = 2, FUN=mean)
rich_means0.005 <- apply(rich0.0005$rich_pi, MARGIN = 2, FUN=mean)

############
############Plot

pdf(width = 7, height = 5, file = "./visuals/richness0.0002_lifehistoryITS.pdf")

stripchart(rich_means0.002~ meta$compartment + meta$lifehistory,
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
           ,ylim=c(1000,3000)
           ,frame.plot=F
)
boxplot(rich_means0.002 ~ meta$compartment + meta$lifehistory, 
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
     at = seq(1000,3000,500),
     lab= seq(1000,3000,500),
     las = 2
)

axis(side =1,
     at = seq(1,8,1),
     lab= rep("", 8)
)
text(x = seq(1,8,1),
     y = 700,
     cex = 1.3,
     pos = 1,
     xpd = NA,
     srt = -50, c("EN forb","EP forb","EN gram.","EP gram.","EN shrub","EP shrub","EN tree","EP tree"))

text("fungi", x =4.5, y = 3500, xpd = NA, cex = 2)

dev.off()


pdf(width = 7, height = 5, file = "./visuals/richness0.0005_lifehistoryITS.pdf")

stripchart(rich_means0.005~ meta$compartment + meta$lifehistory,
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
           ,ylim=c(0,100)
           ,frame.plot=F
)
boxplot(rich_means0.005 ~ meta$compartment + meta$lifehistory, 
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
     at = seq(0,100,20),
     lab= seq(0,100,20),
     las = 2
)

axis(side =1,
     at = seq(1,8,1),
     lab= rep("", 8)
)
text(x = seq(1,8,1),
     y = -15,
     pos = 1,
     cex = 1.3,
     xpd = NA,
     srt = -50, c("EN forb","EP forb","EN gram.","EP gram.","EN shrub","EP shrub","EN tree","EP tree"))

text("fungi", x =4.5, y = 100, xpd = NA, cex = 2)

dev.off()

