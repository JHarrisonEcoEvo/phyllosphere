rm(list=ls())

X<- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv")
X <- X[X$substrate == "plant",]
X$taxon_final<- gsub("var.*","",X$taxon_final)
en <- X[X$compartment == "EN",]
ep <- X[X$compartment == "EP",]


pdf(width = 10, height = 18, file = "./visuals/its_diversity_host.pdf")
par(mfrow = c(2,1), mar = c(12,4,1,1), oma = c(14, 4,4,4))
boxplot(en$div_raw ~ en$taxon_final, main = "EN", xlab = "", las = 2, ylab = "Shannon's equivalents")
boxplot(ep$div_raw ~ ep$taxon_final, main = "EP", xlab = "", las = 2,ylab = "Shannon's equivalents")
dev.off()

#bacteria

rm(list=ls())

X<- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv")
X <- X[X$substrate == "plant",]
X$taxon_final<- gsub("var.*","",X$taxon_final)
en <- X[X$compartment == "EN",]
ep <- X[X$compartment == "EP",]


pdf(width = 10, height = 18, file = "./visuals/16s_diversity_host.pdf")
par(mfrow = c(2,1), mar = c(12,4,1,1), oma = c(14, 4,4,4))
boxplot(en$div_raw ~ en$taxon_final, main = "EN", xlab = "", las = 2, ylab = "Shannon's equivalents")
boxplot(ep$div_raw ~ ep$taxon_final, main = "EP", xlab = "", las = 2,ylab = "Shannon's equivalents")
dev.off()
