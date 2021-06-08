rm(list=ls())
library(tidyr)
library(tibble)
library(dplyr)

dat <- read.csv("./raw_data_notIncluding_sequenceData//rawODKdata/alldata_ODK_combined_edited.csv",
                stringsAsFactors = F, 
                fill = T)
#checking that the number of sites present is as expected
length(unique(dat$siteLabel))

#Goal is a data.frame with site, Daub plot, taxon and abundance class.

#Order of operations:
#1. fill empty rows (currently many blanks)
#2. Daubemire plots need unique id (we get this from lat long data)
#3. Convert cover classes to values between 0-100
#Cover classes are: 0 (0%) 1 (0.1–5%), 2 (5–25%), 3 (25–50%), 4 (50–75%), 5 (75–95%), and 6 (>95%)
#4. Calculate Shannon's entropy for the veg. cover data for each site. This may be a bit tricky.
#6. Export the Shannon's entropy to processed data, then bring it in to the metadata
#during the metadata wrangling scripts.

#1-2. Fill rows (#2 is accomplished here as well)
dat <- as_tibble(dat)

#have to change empty strings to NAs prior to fill
dat[dat == ""] <- NA

dat <- fill(dat,c(2:37, 39:55), .direction="down")

#3. 
dat$daubPlots__VegCompositionBlock__coverClass <- recode(as.character(dat$daubPlots__VegCompositionBlock__coverClass), 
       "1"="0.1",
       "2"="5",
       "3"="25",
       "4"="50",
       "5"="75",
       "6"="95",
       )

dat$daubPlots__daubplotbare <- recode(as.character(dat$daubPlots__daubplotbare), 
                                                         "1"="0.1",
                                                         "2"="5",
                                                         "3"="25",
                                                         "4"="50",
                                                         "5"="75",
                                                         "6"="95",
)
dat$daubPlots__daubplotltter <- recode(as.character(dat$daubPlots__daubplotltter), 
                                      "1"="0.1",
                                      "2"="5",
                                      "3"="25",
                                      "4"="50",
                                      "5"="75",
                                      "6"="95",
)
dat$daubPlots__daubplotgraminoid <- recode(as.character(dat$daubPlots__daubplotgraminoid), 
                                       "1"="0.1",
                                       "2"="5",
                                       "3"="25",
                                       "4"="50",
                                       "5"="75",
                                       "6"="95",
)
dat$daubPlots__daubplotforb <- recode(as.character(dat$daubPlots__daubplotforb), 
                                       "1"="0.1",
                                       "2"="5",
                                       "3"="25",
                                       "4"="50",
                                       "5"="75",
                                       "6"="95",
)
dat$daubPlots__daubplotshrub <- recode(as.character(dat$daubPlots__daubplotshrub), 
                                       "1"="0.1",
                                       "2"="5",
                                       "3"="25",
                                       "4"="50",
                                       "5"="75",
                                       "6"="95",
)
dat$daubPlots__daubplottree <- recode(as.character(dat$daubPlots__daubplottree), 
                                       "1"="0.1",
                                       "2"="5",
                                       "3"="25",
                                       "4"="50",
                                       "5"="75",
                                       "6"="95",
)

#4. Need to convert cover measurements to total percentage for the whole site
#Could combine data for all identical taxa within a site..
#Could convert cover measurements to rows so, four rows for each site, with x columns, for taxa
#then calc shannons for each row and take average shannons...or
#go through by site and do Shannons of the veg composition column for whole site

#remove NAs in the species code
dat <- dat[!is.na(dat$daubPlots__VegCompositionBlock__speciesCode),]

newdat <- as_tibble(data.frame(dat$siteLabel,
                     dat$gpsLocation.Latitude,
                     dat$daubPlots__VegCompositionBlock__speciesCode,
                     dat$daubPlots__VegCompositionBlock__coverClass))
newdat[newdat$dat.siteLabel=="1_1",]

shannons_site_veg <- aggregate(as.numeric(newdat$dat.daubPlots__VegCompositionBlock__coverClass) ~ 
                  newdat$dat.siteLabel, FUN = function(x){
                          vegan::diversity(x, index = "shannon")
                  })

write.csv(shannons_site_veg, 
          file = "./processedData/shannons_siteLevel_veg.csv", row.names = F)
