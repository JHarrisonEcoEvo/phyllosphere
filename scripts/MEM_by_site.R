############################################
# Load data, packages, and specify options #
############################################
rm(list=ls())

set.seed(666)
options(scipen = 999)
library(vegan)
library(adespatial)
library(ade4)
library(adegraphics)
library(raster)
library(geosphere)

dat <- read.csv("./processedData/treatments_metadata.csv")

simpledat <- data.frame(matrix(nrow = 21, ncol = 11))
k <- 1
for(i in unique(dat$unique.dat.siteLabel.)){
  simpledat[k,] <- dat[which(dat$unique.dat.siteLabel. == i)[1], ]
  k <- k + 1
}

Not providing expected values. start here
#haversine distances
dist_mat <- distHaversine(p1 = cbind( 
      simpledat$X4,
      simpledat$X2))


#Run dbmem using defaults for threshold, which is min length in spanning tree
#This computes only MEMs with positive eigenvalues (the function's default)

dbmem_out <-  dbmem(dist_mat, 
                    silent = F) #prints some useful info, including trunc. level

#str(dbmem_out)
length(attributes(dbmem_out)$values) #number of dbMEM eigenvectors

# Plot the associated spatial weighting matrix
s.label(cbind(dat$POINT_X, 
              dat$POINT_Y), 
        nb = attr(dbmem_out, "listw"))

#Plot a specific MEM
#The built in plotting functions are a pain, so use this.

plot( cbind(dat$POINT_X, 
            dat$POINT_Y),
      ylab = "y coord",
      xlab = "x coord",
      cex = abs(dbmem_out$MEM4), 
      col = ifelse(dbmem_out$MEM4 > 0, "red", "gray")
)

# Compute and test associated Moran's I values
# Eigenvalues are proportional to Moran's I

test <- moran.randtest(dbmem_out,
                       nrepet = 99, 
                       alter="two-sided")

plot(test$obs, 
     attr(dbmem_out, "values"), 
     xlab = "Moran's I", 
     ylab = "Eigenvalues")

# Decreasing values of Moran's I for the successive MEM.
# The red line is the expected value of Moran's I under H0.

plot(test$obs, 
     xlab="MEM rank", 
     ylab="Moran's I")
abline(h = -1 / (length(dat$POINT_X) - 1), 
       col="red") 

#Extract significant eigenvectors from Moran's I analysis
dbmem_out[which(test$pvalue < 0.05)]
