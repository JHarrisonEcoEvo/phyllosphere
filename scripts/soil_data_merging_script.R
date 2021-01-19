dat <- read.csv("~/Downloads/soil_sample_data - Sheet1.csv", stringsAsFactors = F)
dim(dat)
dat2 <- read.csv("~/Downloads/Harrison Soil Samples 2018.xlsx - Sheet1.csv")
dim(dat2)
dat3 <- read.csv("~/Downloads/EPSCOR_K2SO4.xlsx - Sheet1.csv")
dim(dat3)

dat$sample_id
dat2$sample_id <- paste(dat2$Sample_ID, dat2$Location, dat2$Horizon, sep = "_")
dat2$sample_id <- gsub("-", "_", dat2$sample_id )
dat2$sample_id <- gsub(",", "_", dat2$sample_id )

dat3$sample_id <- paste(dat3$Sample_ID, dat3$Location, dat3$Horizon, sep = "_")
dat3$sample_id <- gsub("-", "_", dat3$sample_id )
dat3$sample_id <- gsub(",", "_", dat3$sample_id )

dat_1_2 <- merge(dat, dat2, by.x = "sample_id", by.y = "sample_id", all = T)
dim(dat_1_2)
dat_1_2_3 <- merge(dat_1_2, dat3, by.x = "sample_id", by.y = "sample_id", all = T)
dim(dat_1_2_3)

write.csv(dat_1_2_3, file = "~/Desktop/soilDataMerged.csv")
