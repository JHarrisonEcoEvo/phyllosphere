#J. Harrison.
#Univ. of Wyoming.

#NOTE: shannons diversity is the SAME whether it is calculated using ISD normalized data or not. 

#FIRST: make Hellinger table without normalizing via the ISD. This is relic code.
# rm(list=ls())
# 
# meta_samples <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
#                          stringsAsFactors = F)
# meta_samples[which(meta_samples$div_raw != meta_samples$div_raw_norm),]
# 
# #Not normalizing via ISD
# raw <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
#                 stringsAsFactors = F)
# 
# sum(raw[,3:(length(raw)-4)] )
# raw$sample <- gsub("X","", raw$sample)
# 
# raw[,3:(length(raw)-4)] <- vegan::decostand(raw[,3:(length(raw)-4)],
#                                             method = "hellinger")
# 
# #rearrange dat to match the meta_samples ordre
# table(raw$sample %in% meta_samples$sample)
# raw <- raw[match(meta_samples$sample,raw$sample),]
# table(raw$sample == meta_samples$sample)
# # 
# # #calculate diversity 
# # meta_samples$div_raw <- exp(
# #   vegan::diversity(raw[,3:(length(raw)-4)], index = "shannon",MARGIN = 1)
# # ) 
# 
# write.csv(raw, file = "./processedData/hellinger_standardized_16s.csv", row.names = F)

################
#SECOND: normalize by ISD, write table, calculate diversity and write new metadata
rm(list=ls())

meta_samples <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                         stringsAsFactors = F)
meta_samples[which(meta_samples$div_raw != meta_samples$div_raw_norm),]

#Raw count data
raw <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)

sum(raw[,3:(length(raw)-4)] )
raw$sample <- gsub("X","", raw$sample)

#Normalize via ISD
raw[,3:length(raw)] <- raw[,3:length(raw)] / raw$ISD

raw[,3:(length(raw)-4)] <- vegan::decostand(raw[,3:(length(raw)-4)],
                                            method = "hellinger")

#rearrange dat to match the meta_samples ordre
table(raw$sample %in% meta_samples$sample)
raw <- raw[match(meta_samples$sample,raw$sample),]
table(raw$sample == meta_samples$sample)

write.csv(raw, file = "./processedData/hellinger_standardized_16s_normalized.csv", row.names = F)

#calculate diversity 
meta_samples$div_raw <- exp(
  vegan::diversity(raw[,3:(length(raw)-4)], index = "shannon",MARGIN = 1)
) 
#cor.test(meta_samples$div_raw_norm, meta_samples$div_raw) #sanity check, ignoreable

write.csv(meta_samples, file = "./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv", row.names = F)

#QC: we are calculating by samples, which are rows (see MARGIN argument above)
length(meta_samples$div_raw) == dim(raw)[1]

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

#Raw count data
raw <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
raw$sample <- gsub("X","", raw$sample)
sum(raw[,3:(length(raw)-2)])

#Normalize via ISD
raw[,3:length(raw)] <- raw[,3:length(raw)] / raw$ISD

raw[,3:(length(raw)-2)] <- vegan::decostand(raw[,3:(length(raw)-2)],
                                            method = "hellinger")

#rearrange dat to match the meta_samples ordre
table(raw$sample %in% meta_samples$sample)
raw <- raw[match(meta_samples$sample,raw$sample),]
table(raw$sample == meta_samples$sample)

write.csv(raw, file = "./processedData/hellinger_standardized_its_normalized.csv", row.names = F)

#calculate diversity 
meta_samples$div_raw <- exp(
  vegan::diversity(raw[,3:(length(raw)-2)], index = "shannon",MARGIN = 1)
)
#QC: we are calculating by samples, which are rows (see MARGIN argument above)
length(meta_samples$div_raw) == dim(raw)[1]

write.csv(meta_samples, file = "./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", row.names = F)

