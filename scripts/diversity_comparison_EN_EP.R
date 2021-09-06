#J. Harrison.
#Univ. of Wyoming.
rm(list=ls())

meta_samples <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                         stringsAsFactors = F)

reg <- lm(meta_samples$shannonsISD ~ 
            meta_samples$altitude +
            meta_samples$MEM1 +
            meta_samples$MEM2 +
            meta_samples$shannons_flora +
            meta_samples$latitude +
            meta_samples$densitometer +
            meta_samples$compartment +
            meta_samples$plant #This is the only thing that matters.
)
summary(reg) #adj R squared very low, as expected

#assign colors to the metadata
library(viridis)
color_vec <- viridis_pal(option = "viridis")(250)  

colvir <- viridis_pal(option = "viridis")(4)
color_vec[meta_samples$lifehistory == "tree"] <- colvir[1]
color_vec[meta_samples$lifehistory == "shrub"] <- colvir[2]
color_vec[meta_samples$lifehistory == "forb"] <- colvir[3]
color_vec[meta_samples$lifehistory == "graminoid"] <- colvir[4]

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

boxplot(log(meta_samples$shannonsISD) ~
          meta_samples$compartment +meta_samples$lifehistory,
        outline = F)

#Raw count data
raw <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)

sum(raw[,3:(length(raw)-4)] )
raw$sample <- gsub("X","", raw$sample)

raw[,3:(length(raw)-4)] <- vegan::decostand(raw[,3:(length(raw)-4)],
                                            method = "hellinger")

#rearrange dat to match the meta_samples ordre
table(raw$sample %in% meta_samples$sample)
raw <- raw[match(meta_samples$sample,raw$sample),]
table(raw$sample == meta_samples$sample)

write.csv(raw, file = "./processedData/hellinger_standardized_16s.csv", row.names = F)

#calculate diversity 
meta_samples$div_raw <- exp(
  vegan::diversity(raw[,3:(length(raw)-4)], index = "shannon",MARGIN = 1)
)
#QC: we are calculating by samples, which are rows (see MARGIN argument above)
length(meta_samples$div_raw) == dim(raw)[1]

boxplot(meta_samples$div_raw ~
          meta_samples$compartment + meta_samples$lifehistory,
        outline = F)


reg <- lm(meta_samples$div_raw ~ 
            meta_samples$altitude +
            meta_samples$MEM1 +
            meta_samples$MEM2 +
            meta_samples$shannons_flora +
            meta_samples$latitude +
            meta_samples$densitometer +
            meta_samples$compartment +
           meta_samples$plant 
)
summary(reg)


################
# Do for fungi #
################

rm(list=ls())

meta_samples <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                         stringsAsFactors = F)

reg <- lm(meta_samples$shannonsISD ~ 
            meta_samples$altitude +
            meta_samples$MEM1 +
            meta_samples$MEM2 +
            meta_samples$shannons_flora +
            meta_samples$latitude +
            meta_samples$densitometer +
            meta_samples$compartment +
            meta_samples$plant #This is the only thing that matters.
)
summary(reg) #adj R squared very low, but not horrible

#assign colors to the metadata
library(viridis)
color_vec <- viridis_pal(option = "viridis")(250)  

colvir <- viridis_pal(option = "viridis")(4)
color_vec[meta_samples$lifehistory == "tree"] <- colvir[1]
color_vec[meta_samples$lifehistory == "shrub"] <- colvir[2]
color_vec[meta_samples$lifehistory == "forb"] <- colvir[3]
color_vec[meta_samples$lifehistory == "graminoid"] <- colvir[4]

#Function from Mage
add.alpha <- function(col, alpha=1){
  if(missing(col))
    stop("Please provide a vector of colours.")
  apply(sapply(col, col2rgb)/255, 2, 
        function(x) 
          rgb(x[1], x[2], x[3], alpha=alpha))  
}

boxplot((meta_samples$shannonsISD) ~
          meta_samples$compartment + meta_samples$lifehistory,
        outline = F)

#Raw count data
raw <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
raw$sample <- gsub("X","", raw$sample)
sum(raw[,3:(length(raw)-2)])
raw[,3:(length(raw)-2)] <- vegan::decostand(raw[,3:(length(raw)-2)],
                                            method = "hellinger")

#rearrange dat to match the meta_samples ordre
table(raw$sample %in% meta_samples$sample)
raw <- raw[match(meta_samples$sample,raw$sample),]
table(raw$sample == meta_samples$sample)

write.csv(raw, file = "./processedData/hellinger_standardized_its.csv", row.names = F)

#calculate diversity 
meta_samples$div_raw <- exp(
  vegan::diversity(raw[,3:(length(raw)-2)], index = "shannon",MARGIN = 1)
)
#QC: we are calculating by samples, which are rows (see MARGIN argument above)
length(meta_samples$div_raw) == dim(raw)[1]

boxplot(meta_samples$div_raw ~
          meta_samples$compartment + meta_samples$lifehistory,
        outline = F)

reg <- lm(meta_samples$div_raw ~ 
            meta_samples$altitude +
            meta_samples$MEM1 +
            meta_samples$MEM2 +
            meta_samples$shannons_flora +
            meta_samples$latitude +
            meta_samples$densitometer +
            meta_samples$compartment +
            meta_samples$plant 
)
summary(reg)
