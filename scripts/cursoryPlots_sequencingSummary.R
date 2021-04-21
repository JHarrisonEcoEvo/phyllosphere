dat <- read.csv("processedData/otuTables/smallmem97_16s_for_modeling", stringsAsFactors = F)

bacteria <- colSums(dat[!(dat$otus %in% c("duds", "ISD","mtDNA","plant")),
            3:length(dat)])
table(bacteria > 1000)
extraneous <- colSums(dat[(dat$otus %in% c("duds", "ISD","mtDNA","plant")),
                        3:length(dat)])
pdf(width = 10, height = 8, file = "./visuals/bacterial_proportion_reads.pdf")
par(mfrow = c(2,3))
hist(bacteria / (bacteria + extraneous), main = "Bacterial prop. of reads",
     xlab = "proportion")
hist(log(bacteria + 1), 
     main = "bacteria (ln + 1)",
     xlab = "ln")
hist(log(unlist(dat[dat$otus == "ISD", 3:length(dat)])+1), 
     main = "ISD (ln + 1)",
     xlab = "ln")
hist(log(unlist(dat[dat$otus == "duds", 3:length(dat)])+1), 
     main = "Unknown (ln + 1)",
     xlab = "ln")
hist(log(unlist(dat[dat$otus == "plant", 3:length(dat)])+1), 
     main = "cpDNA (ln + 1)",
     xlab = "ln")
hist(log(unlist(dat[dat$otus == "mtDNA", 3:length(dat)])+1), 
     main = "mtDNA (ln + 1)",
     xlab = "ln")
dev.off()

epdat <- dat[,c(2, grep("EP", names(dat)))]
endat <- dat[,c(2, grep("EN", names(dat)))]

bacteria_en <- colSums(endat[!(endat$otus %in% c("duds", "ISD","mtDNA","plant")),
                        2:length(endat)])
extraneous_en <- colSums(endat[(endat$otus %in% c("duds", "ISD","mtDNA","plant")),
                          2:length(endat)])
bacteria_ep <- colSums(epdat[!(epdat$otus %in% c("duds", "ISD","mtDNA","plant")),
                             2:length(epdat)])
extraneous_ep <- colSums(epdat[(epdat$otus %in% c("duds", "ISD","mtDNA","plant")),
                               2:length(epdat)])
pdf(width = 10, height = 8, file = "./visuals/en_vs_ep_bacterialReads.pdf")

par(mfrow = c(1,2))
hist(bacteria_en / (bacteria_en + extraneous_en), 
     main = "Endophyte - bacteria",
     xlab = "proportion")
hist(bacteria_ep / (bacteria_ep + extraneous_ep), 
     main = "Epiphyte - bacteria",
     xlab = "proportion")
dev.off()

#######
# ITS #
#######
rm(list = ls())
dat <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling", stringsAsFactors = F)

fungi <- colSums(dat[!(dat$otus %in% c("duds", "ISD","mtDNA","plant")),
                        3:length(dat)])
extraneous <- colSums(dat[(dat$otus %in% c("duds", "ISD","mtDNA","plant")),
                          3:length(dat)])
pdf(width = 10, height = 8, file = "./visuals/fungal_proportion_reads.pdf")
par(mfrow = c(2,2))
hist(fungi / (fungi + extraneous), main = "Fungal prop. of reads",
     xlab = "proportion")
hist(log(fungi + 1), 
     main = "fungi (ln + 1)",
     xlab = "ln")
hist(log(unlist(dat[dat$otus == "ISD", 3:length(dat)])+1), 
     main = "ISD (ln + 1)",
     xlab = "ln")
hist(log(unlist(dat[dat$otus == "duds", 3:length(dat)])+1), 
     main = "Unknown (ln + 1)",
     xlab = "ln")
dev.off()

epdat <- dat[,c(2, grep("EP", names(dat)))]
endat <- dat[,c(2, grep("EN", names(dat)))]
table(fungi > 1000)

fungi_en <- colSums(endat[!(endat$otus %in% c("duds", "ISD","mtDNA","plant")),
                             2:length(endat)])
table(fungi_en > 1000)
table(fungi_en < 100)
extraneous_en <- colSums(endat[(endat$otus %in% c("duds", "ISD","mtDNA","plant")),
                               2:length(endat)])
fungi_ep <- colSums(epdat[!(epdat$otus %in% c("duds", "ISD","mtDNA","plant")),
                             2:length(epdat)])
table(fungi_ep > 1000)
table(fungi_ep < 100)

extraneous_ep <- colSums(epdat[(epdat$otus %in% c("duds", "ISD","mtDNA","plant")),
                               2:length(epdat)])
pdf(width = 10, height = 8, file = "./visuals/en_vs_ep_fungalReads.pdf")

par(mfrow = c(1,2))
hist(fungi_en / (fungi_en + extraneous_en), 
     main = "Endophyte - fungi",
     xlab = "proportion")
hist(fungi_ep / (fungi_ep + extraneous_ep), 
     main = "Epiphyte - fungi",
     xlab = "proportion")
dev.off()
