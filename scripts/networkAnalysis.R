rm(list=ls())
datf <- read.csv("processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                 stringsAsFactors = F)
datf[,3:length(datf)] <- datf[,3:length(datf)] - 1
datf[,3:length(datf)] <- datf[,3:length(datf)] / datf$ISD

for(i in 3:length(datf)){
  datf[is.infinite(datf[,i]),i] <- NA
}
datf$sample <- gsub("X","",datf$sample)

X <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                 stringsAsFactors = F)

X <- X[X$substrate == "plant",]

#Make a leaf density variable
X$sla = X$mass_extracted_g / X$area_cm2

X$sla[is.infinite(X$sla)] <- NA

#Convert partial leaf to 0.5
X$leaves_extracted <- (as.character(X$leaves_extracted))
X$leaves_extracted[X$leaves_extracted == "partial"] <- 0.5
X$leaves_extracted <- (as.numeric(X$leaves_extracted))

X$waterRetention <- (as.character(X$waterRetention))
indices <- which(X$waterRetention == "90+")
X$waterRetention[indices] <- 90
X$waterRetention <- (as.numeric(X$waterRetention))

#remove spaces in taxon_final for ease of output file manipulation. 
X$taxon_final <- gsub(" ", "", X$taxon_final)

#the missing stuff were things I removed due to poor sequencing
table(X$sample %in% datf$sample)

merged_dat <- merge(X, datf, 
                    by.x = "sample", 
                    by.y = "sample", 
                    all.y = T)
newdf <- data.frame(
  merged_dat$sample, 
  merged_dat$habit,
  merged_dat$taxon_final, 
  merged_dat[,206:(length(merged_dat)-2)]
  )

library(tidyr)
#I f'ng hate switching from wide to long and back.
data_long <- gather(newdf, 
                    Microbe, 
                    abundance, 
                    Zotu10034:Zotu9997)
data_long
table(data_long$merged_dat.sample)
dim( merged_dat[,206:(length(merged_dat)-2)])

#QC look for a specific sample and make sure it is the same in the wide and long data. 
# tdat <-t(data_long[data_long$merged_dat.sample == "3_2_6_8_EP",3:4])
# test <- newdf[newdf$merged_dat.sample == "3_2_6_8_EP", grep("Zotu", names(newdf))]
# test <- test[match(names(test), tdat[1,])]
# 
# table((tdat[1,]) == names(test))
# table(as.numeric(tdat[2,]) == as.numeric(test))
# 
# #inexplicably the stuff that is supposed to be different doesn't appear to be. 
# #must be a rounding thing.
# (cbind(test[test != as.numeric(tdat[2,]) ],
#       as.numeric(tdat[2,])[test != as.numeric(tdat[2,]) ]))
# 
# test[test != as.numeric(tdat[2,]) ] ==
#       as.numeric(tdat[2,])[test != as.numeric(tdat[2,]) ]
# 
# #try rounding. yup. fuck.
# round(test[test != as.numeric(tdat[2,]) ],10) ==
#   round(as.numeric(tdat[2,])[test != as.numeric(tdat[2,]) ],10)

#SUBSET THE DATA, so this doesn't take forever and crash my computer
#data_long <- data_long[
#  data_long$merged_dat.habit == "tree",]
data_long_abbrev <- data_long[
  as.numeric(data_long$abundance) != 0,]
dim(data_long_abbrev)

library(bipartite)
data_long <- data.frame(data_long_abbrev[,3:4],
                        data_long_abbrev[,2],
                       data_long_abbrev[,5])
names(data_long) <- c("higher","lower","webID", "freq")
stupidNet <- frame2webs(data_long, type.out = "array", emptylist = F)
plotweb(stupidNet)

#This fails and I think the bipartite package may have reached the end of its useful life. IT is very slow.
