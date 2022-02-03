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

simpledat <- data.frame(matrix(nrow = 20, ncol = 13))
k <- 1
for(i in unique(dat$label)){
  simpledat[k,] <- dat[which(dat$label == i)[1], ]
  k <- k + 1
}

#haversine distances
dist_mat <- distm(x = cbind( 
      simpledat$X4,
      simpledat$X2))

#Run dbmem using defaults for threshold, which is min length in spanning tree
#This computes only MEMs with positive eigenvalues (the function's default)

dbmem_out <-  dbmem(as.dist(dist_mat), 
                    silent = F) #prints some useful info, including trunc. level

#str(dbmem_out)
length(attributes(dbmem_out)$values) #number of dbMEM eigenvectors

# Plot the associated spatial weighting matrix
s.label(cbind( 
  simpledat$X4,
  simpledat$X2), 
        nb = attr(dbmem_out, "listw"))

#Plot a specific MEM
#The built in plotting functions are a pain, so use this.

plot( cbind( 
  simpledat$X4,
  simpledat$X2),
      ylab = "y coord",
      xlab = "x coord",
      cex = abs(dbmem_out$MEM1), 
      col = ifelse(dbmem_out$MEM1 > 0, "red", "gray")
)

# Compute and test associated Moran's I values
# Eigenvalues are proportional to Moran's I

test <- moran.randtest(dbmem_out,
                       nrepet = 99, 
                       alter="two-sided")

#Extract significant eigenvectors from Moran's I analysis
dbmem_out[which(test$pvalue < 0.05)]


write.csv(data.frame(unique(dat$label), 
                     dbmem_out), 
          file = "./processedData/moransEigenVectors.csv",
          row.names = F)
