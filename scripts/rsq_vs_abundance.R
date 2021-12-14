
#Calculate abundance of modeled microbes
rm(list=ls())
dat <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

focal_taxa <- read.csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv",
                       stringsAsFactors = F)

abundance <- apply(dat[,names(dat) %in% focal_taxa$x], 2, FUN=mean, na.rm = T)
abundance[is.infinite(abundance)]<- NA

rsq <- read.csv("modelingResults/results_landscape_isd/all_16S.csv",
                stringsAsFactors = F)
rsq <- rsq[rsq$taxon != "taxon",]

#checking arrangement
rsq <- rsq[match(names(abundance),rsq$taxon ),]
table(rsq$taxon == names(abundance))

pdf(width = 8, height = 5, file = "./visuals/rsq_abundance.pdf")
par(mfrow = c(1,2), oma = c(1,1,1,1))
plot(as.numeric(rsq$rsq_nested_resampling), 
     abundance, main = "bacteria", xlab = "R^2", ylab = "",yaxt = "n", pch = 16, frame.plot = F)
axis(side = 2, las = 2)
title(ylab="abundance", line=3.5, cex.lab=1.2, xpd = NA)

corout1 <- cor.test(as.numeric(rsq$rsq_nested_resampling), 
                    abundance)

#DO FOR FUNGI
#Calculate prevalence of modeled microbes
dat <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

focal_taxa <- read.csv("./processedData/ITS_taxa_to_model_via_randomforest.csv",
                       stringsAsFactors = F)

abundance <- apply(dat[,names(dat) %in% focal_taxa$x], 2, FUN=mean, na.rm = T)
abundance[is.infinite(abundance)]<- NA

rsq <- read.csv("modelingResults/results_landscape_isd/all_ITS.csv",
                stringsAsFactors = F)
rsq <- rsq[rsq$taxon != "taxon",]

#checking arrangement

rsq <- rsq[match(names(abundance),rsq$taxon),]
table(rsq$taxon == names(abundance))

plot(as.numeric(rsq$rsq_nested_resampling), 
     abundance, main = "fungi", ylab = "", xlab = "R^2", yaxt = "n",pch = 16, frame.plot = F)
axis(side = 2, las = 2)
title(ylab="abundance", line=3.5, cex.lab=1.2, xpd = NA)

corout <- cor.test(as.numeric(rsq$rsq_nested_resampling), 
                    abundance)

text(paste("R = ", round(corout$estimate,2), sep = ""), x = -600, y = 0.2, xpd = NA)
text(paste("R = ", round(corout1$estimate,2), sep = ""), x = -2100, y = 0.2, xpd = NA)

dev.off()

###################
######################################
###################
###################
###################

#Do over for hellinger data
#Calculate abundance of modeled microbes
rm(list=ls())
dat <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

focal_taxa <- read.csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv",
                       stringsAsFactors = F)

abundance <- apply(dat[,names(dat) %in% focal_taxa$x], 2, FUN=mean, na.rm = T)
abundance[is.infinite(abundance)]<- NA

rsq <- read.csv("modelingResults/results_hella_noduds///all_16s_hella_noduds.csv",
                stringsAsFactors = F)
rsq <- rsq[rsq$taxon != "taxon",]

#checking arrangement
rsq <- rsq[match(names(abundance),rsq$taxon ),]
table(rsq$taxon == names(abundance))

pdf(width = 8, height = 5, file = "./visuals/rsq_abundance_hellinger16s.pdf")
par(mfrow = c(1,2), oma = c(1,2,1,1))
plot(as.numeric(rsq$rsq_nested_resampling), 
     abundance, main = "bacteria", xlab = "R^2", ylab = "",yaxt = "n", pch = 16, frame.plot = F)
axis(side = 2, las = 2)
title(ylab="abundance", line=3.5, cex.lab=1.2, xpd = NA)

corout1 <- cor.test(as.numeric(rsq$rsq_nested_resampling), 
                    abundance)

#DO FOR FUNGI
#Calculate prevalence of modeled microbes
dat <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

focal_taxa <- read.csv("./processedData/ITS_taxa_to_model_via_randomforest.csv",
                       stringsAsFactors = F)

abundance <- apply(dat[,names(dat) %in% focal_taxa$x], 2, FUN=mean, na.rm = T)
abundance[is.infinite(abundance)]<- NA

rsq <- read.csv("modelingResults/results_hella_noduds/all_its_hella_noduds.csv",
                stringsAsFactors = F)
rsq <- rsq[rsq$taxon != "taxon",]

#checking arrangement
rsq <- rsq[match(names(abundance),rsq$taxon),]
table(rsq$taxon == names(abundance))

plot(as.numeric(rsq$rsq_nested_resampling), 
     abundance, main = "fungi", ylab = "", xlab = "R^2", yaxt = "n",pch = 16, frame.plot = F)
axis(side = 2, las = 2)
title(ylab="abundance", line=3.5, cex.lab=1.2, xpd = NA)

corout <- cor.test(as.numeric(rsq$rsq_nested_resampling), 
                   abundance)

text(paste("R = ", round(corout$estimate,2), sep = ""), x = 0.3, y = 0.2, xpd = NA)
text(paste("R = ", round(corout1$estimate,2), sep = ""), x = -0.23, y = 0.2, xpd = NA)

dev.off()

