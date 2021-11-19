
#Calculate prevalence of modeled microbes
rm(list=ls())
dat <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG_OCCUPANCY.csv",
                stringsAsFactors = F)
focal_taxa <- read.csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv",
                       stringsAsFactors = F)

prevalence <- colSums(dat[,names(dat) %in% focal_taxa$x])

mcc <- read.csv("modelingResults/results_landscape_isd/all_16s_isd.csv",
                stringsAsFactors = F)
mcc <- mcc[mcc$taxon != "taxon",]

#checking arrangement
mcc <- mcc[match(mcc$taxon, names(prevalence)),]
table(mcc$taxon == names(prevalence))

pdf(width = 8, height = 4, file = "./visuals/rsq_prevalence.pdf")
par(mfrow = c(1,2))
plot(as.numeric(mcc$rsq_nested_resampling), 
     prevalence, main = "bacteria", xlab = "R^2", yaxt = "n", pch = 16, frame.plot = F)
axis(side = 2, las = 2)

corout1 <- cor.test(as.numeric(mcc$rsq_nested_resampling), 
                    prevalence)

#DO FOR FUNGI
#Calculate prevalence of modeled microbes
dat <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG_OCCUPANCY.csv",
                stringsAsFactors = F)
focal_taxa <- read.csv("./processedData/ITS_taxa_to_model_via_randomforest.csv",
                       stringsAsFactors = F)

prevalence <- colSums(dat[,names(dat) %in% focal_taxa$x])

mcc <- read.csv("modelingResults/results_landscape_isd/all_its_isd.csv",
                stringsAsFactors = F)
mcc <- mcc[mcc$taxon != "taxon",]

#checking arrangement
mcc <- mcc[match(names(prevalence),mcc$taxon),]
table(mcc$taxon == names(prevalence))

plot(as.numeric(mcc$rsq_nested_resampling), 
     prevalence, main = "fungi", xlab = "R^2", yaxt = "n",pch = 16, frame.plot = F)
axis(side = 2, las = 2)

corout <- cor.test(as.numeric(mcc$rsq_nested_resampling), 
                   prevalence)
text(paste("R = ", round(corout$estimate,2), sep = ""), x = -1500, y = 1400)
text(paste("R = ", round(corout1$estimate,2), sep = ""), x = -3000, y = 1400, xpd = NA)

dev.off()
