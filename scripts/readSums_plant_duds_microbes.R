
#Calculate abundance of modeled microbes
rm(list=ls())
dat <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

focal_taxa <- read.csv("./processedData/sixteenS_taxa_to_model_via_randomforest.csv",
                       stringsAsFactors = F)

#xtract total for the focal taxa
colsumsfocal <- colSums(dat[,names(dat) %in% focal_taxa$x], na.rm = T)
colsumsfocal[is.infinite(colsumsfocal)] <- NA
colsumsfocal <- sum(colsumsfocal, na.rm = T)

#xtract total for the other taxa
dat_nonfocal <- dat[,!names(dat) %in% focal_taxa$x]
dat_nonfocal <- dat_nonfocal[,!names(dat_nonfocal) %in% c("plant", "duds","mtDNA","ISD")]

colsumsUNfocal <- colSums(dat_nonfocal[,3:length(dat_nonfocal)], na.rm = T)
colsumsUNfocal[is.infinite(colsumsUNfocal)] <- NA
colsumsUNfocal <- sum(colsumsUNfocal, na.rm = T)

dat$plant[is.infinite(dat$plant)] <- NA
plantcounts <- sum(dat$plant, na.rm = T)

dat$mtDNA[is.infinite(dat$mtDNA)] <- NA
mtDNAcounts <- sum(dat$mtDNA, na.rm = T)

dat$duds[is.infinite(dat$duds)] <- NA
dudscounts <- sum(dat$duds, na.rm = T)



#DO FOR FUNGI
#Calculate prevalence of modeled microbes
datf <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
datf[,3:length(datf)] <- datf[,3:length(datf)] - 1
datf[,3:length(datf)] <- datf[,3:length(datf)] / datf$ISD

datf$duds[is.infinite(datf$duds)] <- NA
dudscountsf <- sum(datf$duds, na.rm = T)

focal_taxa <- read.csv("./processeddata/ITS_taxa_to_model_via_randomforest.csv",
                       stringsAsFactors = F)


#xtract total for the focal taxa
colsumsfocalf <- colSums(datf[,names(datf) %in% focal_taxa$x], na.rm = T)
colsumsfocalf[is.infinite(colsumsfocalf)] <- NA
colsumsfocalf <- sum(colsumsfocalf, na.rm = T)

#xtract total for the other taxa
datf_nonfocalf <- datf[,!names(datf) %in% focal_taxa$x]
datf_nonfocalf <- datf_nonfocalf[,!names(datf_nonfocalf) %in% c("plant", "duds","mtDNA","ISD")]

colsumsUNfocalf <- colSums(datf_nonfocalf[,3:length(datf_nonfocalf)], na.rm = T)
colsumsUNfocalf[is.infinite(colsumsUNfocalf)] <- NA
colsumsUNfocalf <- sum(colsumsUNfocalf, na.rm = T)

#make a stacked barplot

library(ggplot2)
library(gridExtra)

# create a dataset
taxa <- c("bacteria","fungi")
group_count <- c("modeled" , "rare" , "plant", "Unknown", "mtDNA")
value <- (c(colsumsfocal, colsumsUNfocal, plantcounts, dudscounts, mtDNAcounts, #bacteria
           colsumsfocalf, colsumsUNfocalf, 0, dudscountsf, 0)) #fungi
value[is.infinite(value)] <- NA
data <- data.frame(taxa, group_count, value)


# Stacked
p1 <- ggplot(data, aes(fill=group_count, y=value, x=taxa)) + 
  geom_bar(position="stack", stat="identity",show.legend = F)+
  xlab("")+
  ylab("abund (div. by ISD)")

#logged version
p2 <- ggplot(data, aes(fill=group_count, y=log(value), x=taxa)) + 
  geom_bar(position="stack", stat="identity")+
  xlab("")+
  ylab("ln abund.")
p2 <- p2 + guides(fill=guide_legend(title=""))

#proportion version
p3 <- ggplot(data, aes(fill=group_count, y=value, x=taxa)) + 
  geom_bar(position="fill", stat="identity",show.legend = F) +
  xlab("")+
  ylab("proportion")


library(egg)
pdf(width = 8, height = 7, file = "./visuals/stackedBarplots_readAbund.pdf", onefile = F)
ggarrange(p1, p2, p3, widths = c(1.5,2))
dev.off()


###################################
#plant microbe ratio by host taxon#
rm(list=ls())
dat <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

dat$sample <- gsub("X","", dat$sample)
meta <- read.csv("processedData/16smetadat_wrangled_for_post_modeling_analysis.csv",
                 stringsAsFactors = F)
meta <- meta[match(dat$sample,meta$sample),]
table(dat$sample == meta$sample)

dat$allbact <- rowSums(dat[,3:(length(dat)-4)] )

dat$b_p_ratio <- dat$allbact / dat$plant

#meta is in the same order as dat, see above
out <- aggregate(dat$b_p_ratio ~ meta$taxon_final, FUN = "median")
pdf(file = "visuals/16s_cpdna.pdf", width = 12, height = 8)
par(mar = c(10,3,3,3), oma = c(15,1,1,1))

plot(out$`dat$b_p_ratio`, pch = 16, col = "orange", cex = 2, ylab = "Ratio of bacteria to plant",
     xlab = "",xaxt = "n", frame.plot = F, xpd = NA)
axis(side = 1, at = c(1:length(out[,1])),
     labels = out$`meta$taxon_final`, las = 2)
dev.off()

#FUNGI
rm(list=ls())
dat <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                stringsAsFactors = F)
dat[,3:length(dat)] <- dat[,3:length(dat)] - 1
dat[,3:length(dat)] <- dat[,3:length(dat)] / dat$ISD

dat$sample <- gsub("X","", dat$sample)
meta <- read.csv("processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                 stringsAsFactors = F)
meta <- meta[match(dat$sample,meta$sample),]
table(dat$sample == meta$sample)

dat$allbact <- rowSums(dat[,3:(length(dat)-2)] )

dat$b_p_ratio <- dat$allbact / dat$duds

#meta is in the same order as dat, see above
out <- aggregate(dat$b_p_ratio ~ meta$taxon_final, FUN = "median")
pdf(file = "visuals/its_duds.pdf", width = 12, height = 8)
par(mar = c(10,3,3,3), oma = c(15,1,1,1))

plot(out$`dat$b_p_ratio`, pch = 16, col = "orange", cex = 2, ylab = "Ratio of fungi to duds",
     xlab = "",xaxt = "n", frame.plot = F, xpd = NA)
axis(side = 1, at = c(1:length(out[,1])),
     labels = out$`meta$taxon_final`, las = 2)
dev.off()


