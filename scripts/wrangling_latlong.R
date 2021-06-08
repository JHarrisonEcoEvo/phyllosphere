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

latlong <- data.frame(unique(dat$siteLabel),
  unique(dat$gpsLocation.Latitude),
unique(dat$gpsLocation.Altitude),
unique(dat$gpsLocation.Longitude)
)

write.csv(latlong, 
          file = "./processedData/latlong_elevation_by_site.csv", row.names = F)
