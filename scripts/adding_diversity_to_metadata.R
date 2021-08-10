rm(list=ls())

#Bring in data, packages, set up color scheme
its <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv")
#its_nonISd <- read.csv("./processedData/ITS_multinomial_p_estimates.csv")
meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", stringsAsFactors = F, header = T)

dim(its)
dim(meta)

for(i in 1:dim(its)[1]){
  meta$shannonsISD[i] <- exp(vegan::diversity(its[i,grep("Zotu", names(its))], index = "shannon"))
  meta$simpsonsISD[i] <- 1/(vegan::diversity(its[i,grep("Zotu", names(its))], index = "simpson"))
 # meta$shannons[i] <- exp(vegan::diversity(its_nonISd[i,grep("Zotu", names(its_nonISd))], index = "shannon"))
 # meta$simpsons[i] <- 1/(vegan::diversity(its_nonISd[i,grep("Zotu", names(its_nonISd))], index = "simpson"))
}
write.csv(meta, file = "./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",row.names = F)


#######
# 16S #
#######

rm(list=ls())

#Bring in data, packages, set up color scheme
sixs <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv")
#sixs_nonISd <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis.csv")
meta <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv", stringsAsFactors = F, header = T)

for(i in 1:dim(sixs)[1]){
  meta$shannonsISD[i] <- exp(vegan::diversity(sixs[i,grep("Zotu", names(sixs))], index = "shannon"))
  meta$simpsonsISD[i] <- 1/(vegan::diversity(sixs[i,grep("Zotu", names(sixs))], index = "simpson"))
#  meta$shannons[i] <- exp(vegan::diversity(sixs_nonISd[i,grep("Zotu", names(sixs_nonISd))], index = "shannon"))
#  meta$simpsons[i] <- 1/(vegan::diversity(sixs_nonISd[i,grep("Zotu", names(sixs_nonISd))], index = "simpson"))
}
write.csv(meta, file = "./processedData/16Smetadat_wrangled_for_post_modeling_analysis.csv",row.names = F)