rm(list=ls())
meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv")

library(wesanderson)
colors <- wes_palette("FantasticFox1", n = 5, type = "discrete")

color_vec <- vector(length = length(meta$lifehistory))
color_vec[meta$lifehistory == "tree"] <- colors[1]
color_vec[meta$lifehistory == "shrub"] <- colors[2]
color_vec[meta$lifehistory == "forb"] <- colors[3]
color_vec[meta$lifehistory == "graminoid"] <- colors[5]

colvir<- colors[c(1,2,3,5)]

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

#PLOT raw_isd data

pdf(width = 7, height = 7, file = "./visuals/diversity_density_lifehistoryITS_raw_isd.pdf")

par(oma = c(4,4,2,2))
stripchart((meta$div_raw) ~ meta$compartment + meta$lifehistory,
           vertical = TRUE, 
           
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colvir[c(1,1,2,2,3,3,4,4)], alpha=0.7),
           cex=1,
           xaxt="n",
           yaxt="n",
           cex.lab = 1.2,
           ylab="Shannon's equivalents",
           line = 4,
           xlim=c(0,8.5)
           ,ylim=c(0,3500)
           #,ylim=c(1000,3500)
           ,frame.plot=F
)
boxplot(meta$div_raw  ~ meta$compartment + meta$lifehistory, 
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
     at = seq(0,3500,500),
     lab= seq(0,3500,500),
     las = 2
)

axis(side =1,
     at = seq(1,8,1),
     lab= rep("", 8)
)
text(x = seq(1,8,1),
     y = -450,
     cex = 1.3,
     xpd = NA,
     pos = 1,
     srt = -50, c("EN forb","EP forb","EN gram.","EP gram.","EN shrub","EP shrub","EN tree","EP tree"))

text(x = 4.5,
     y = 3500,
     cex = 1,
     xpd = NA,
     "Fungi (raw_isd, ISD normalized)"
)

dev.off()

table(is.na(meta$lifehistory))
meta <- meta[!is.na(meta$lifehistory),]
TukeyHSD(aov(meta$div_raw ~ meta$compartment + meta$lifehistory))
xtable::xtable(TukeyHSD(aov(meta$div_raw ~ paste(meta$compartment, meta$lifehistory)))$`paste(meta$compartment, meta$lifehistory)`)

out <- list()
k <- 1
for(i in na.omit(unique(meta$lifehistory))){
  out[[k]] <- t.test(meta$div_raw[meta$lifehistory == i & meta$compartment == "EN"],
    meta$div_raw[meta$lifehistory == i & meta$compartment == "EP"])
    k <- k + 1
  }

############
############Plot modeled data
# 
# pdf(width = 7, height = 7, file = "./visuals/diversity_density_lifehistoryITS_modeled.pdf")
# par(oma = c(4,4,2,2))
# 
# stripchart((meta$shannonsISD) ~ meta$compartment + meta$lifehistory,
#            vertical = TRUE, 
#        
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha(colvir[c(1,1,2,2,3,3,4,4)], alpha=0.7),
#            cex=1,
#            xaxt="n",
#            yaxt="n",
#            cex.lab = 1.2,
#            ylab="Shannon's equivalents",
#            xlim=c(0,8.5)
#            ,ylim=c(0,3500)
#            #,ylim=c(1000,3500)
#            ,frame.plot=F
# )
# boxplot(meta$shannonsISD  ~ meta$compartment + meta$lifehistory, 
#         na.omit = T,
#         add=T,
#         las=2,
#         outline=F,
#         axes=F,
#         col=c(add.alpha(colvir[1], alpha=0.3),
#               add.alpha(colvir[1], alpha=0.6),
#               add.alpha(colvir[2], alpha=0.3),
#               add.alpha(colvir[2], alpha=0.6),
#               add.alpha(colvir[3], alpha=0.3),
#               add.alpha(colvir[3], alpha=0.6),
#               add.alpha(colvir[4], alpha=0.3),
#               add.alpha(colvir[4], alpha=0.6)),
#         ylab=""
# )
# 
# axis(side =2,
#      at = seq(0,3500,500),
#      lab= seq(0,3500,500),
#      las = 2
# )
# 
# axis(side =1,
#      at = seq(1,8,1),
#      lab= rep("", 8)
# )
# text(x = seq(1,8,1),
#      y = -500,
#      pos = 1,
#      cex = 1.3,
#      xpd = NA,
#      srt = -50, c("EN forb","EP forb","EN gram.","EP gram.","EN shrub","EP shrub","EN tree","EP tree"))
# 
# text(x = 4.5,
#      y = 3500,
#      cex = 1,
#      xpd = NA,
#      "Fungi (modeled, ISD normalized)"
# )
# 
# dev.off()


#################
#DO FOR BACTERIA#
#################
rm(list=ls())
meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv")

library(wesanderson)
colors <- wes_palette("FantasticFox1", n = 5, type = "discrete")
color_vec <- vector(length = length(meta$lifehistory))

color_vec[meta$lifehistory == "tree"] <- colors[1]
color_vec[meta$lifehistory == "shrub"] <- colors[2]
color_vec[meta$lifehistory == "forb"] <- colors[3]
color_vec[meta$lifehistory == "graminoid"] <- colors[5]

colvir<- colors[c(1,2,3,5)]
#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}
#PLOT raw_isd data

pdf(width = 7, height = 7, file = "./visuals/diversity_density_lifehistory16S_raw_isd.pdf")

par(oma = c(4,4,2,2))
stripchart((meta$div_raw) ~ meta$compartment + meta$lifehistory,
           vertical = TRUE, 
           
           method = "jitter", 
           pch = 20, 
           col = add.alpha(colvir[c(1,1,2,2,3,3,4,4)], alpha=0.7),
           cex=1,
           xaxt="n",
           yaxt="n",
           cex.lab = 1.2,
           ylab="",
           line = 3.5,
           #xpd = NA,
           xlim=c(0,8.5)
           ,ylim=c(2000,2400)
           #,ylim=c(1000,3500)
           ,frame.plot=F
)
boxplot(meta$div_raw  ~ meta$compartment + meta$lifehistory, 
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
     at = seq(2000,2400,100),
     lab= seq(2000,2400,100),
     las = 2
)

axis(side =1,
     at = seq(1,8,1),
     lab= rep("", 8)
)
text(x = seq(1,8,1),
     y = 1950,
     cex = 1.3,
     pos = 1,
     xpd = NA,
     srt = -50, c("EN forb","EP forb","EN gram.","EP gram.","EN shrub","EP shrub","EN tree","EP tree"))

text(x = 4.5,
     y = 2400,
     cex = 1,
     xpd = NA,
     "Bacteria (raw_isd, ISD normalized)"
)

text(
  "Shannon's equivalents",
  y = 2175,
  x = -2,
  srt = 90,
  xpd = NA
)

dev.off()

out <- list()
k <- 1
for(i in na.omit(unique(meta$lifehistory))){
  out[[k]] <- t.test(meta$div_raw[meta$lifehistory == i & meta$compartment == "EN"],
                     meta$div_raw[meta$lifehistory == i & meta$compartment == "EP"])
  k <- k + 1
}

table(is.na(meta$lifehistory))
meta <- meta[!is.na(meta$lifehistory),]
TukeyHSD(aov(meta$div_raw ~ meta$compartment + meta$lifehistory))
xtable::xtable(TukeyHSD(aov(meta$div_raw ~ paste(meta$compartment, meta$lifehistory)))$`paste(meta$compartment, meta$lifehistory)`)

############
############Plot modeled data
# 
# pdf(width = 7, height = 7, file = "./visuals/diversity_density_lifehistory16S_modeled.pdf")
# 
# par(oma = c(4,4,2,2))
# stripchart((meta$shannonsISD) ~ meta$compartment + meta$lifehistory,
#            vertical = TRUE, 
#            
#            method = "jitter", 
#            pch = 20, 
#            col = add.alpha(colvir[c(1,1,2,2,3,3,4,4)], alpha=0.7),
#            cex=1,
#            xaxt="n",
#            yaxt="n",
#            cex.lab = 1.2,
#            ylab="",
#            line = 3.5,
#            #xpd = NA,
#            xlim=c(0,8.5)
#            ,ylim=c(2100,2360)
#            #,ylim=c(1000,3500)
#            ,frame.plot=F
# )
# boxplot(meta$shannonsISD  ~ meta$compartment + meta$lifehistory, 
#         na.omit = T,
#         add=T,
#         las=2,
#         outline=F,
#         axes=F,
#         col=c(add.alpha(colvir[1], alpha=0.3),
#               add.alpha(colvir[1], alpha=0.6),
#               add.alpha(colvir[2], alpha=0.3),
#               add.alpha(colvir[2], alpha=0.6),
#               add.alpha(colvir[3], alpha=0.3),
#               add.alpha(colvir[3], alpha=0.6),
#               add.alpha(colvir[4], alpha=0.3),
#               add.alpha(colvir[4], alpha=0.6)),
#         ylab=""
# )
# 
# axis(side =2,
#      at = seq(2100,2350,50),
#      lab= seq(2100,2350,50),
#      las = 2
# )
# 
# axis(side =1,
#      at = seq(1,8,1),
#      lab= rep("", 8)
# )
# text(x = seq(1,8,1),
#      y = 2070,
#      cex = 1.3,
#      pos = 1,
#      xpd = NA,
#      srt = -50, c("EN forb","EP forb","EN gram.","EP gram.","EN shrub","EP shrub","EN tree","EP tree"))
# 
# text(x = 4.5,
#      y = 2365,
#      cex = 1,
#      xpd = NA,
#      "Bacteria (modeled, ISD normalized)"
# )
# 
# text(
#   "Shannon's equivalents",
#   y = 2225,
#   x = -2,
#   srt = 90,
#   xpd = NA
# )
# 
# dev.off()

