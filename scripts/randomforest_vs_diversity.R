#install.packages("randomForest")
library(randomForest)
library(Metrics)
options(scipen = 99)


meta1 <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", 
                 stringsAsFactors = F, header = T)

predictors<-c(
  "area_cm2",
  "mass_g",
  "mass_extracted_g",
  #  "leaves_extracted",
  "circumStem",
  "width",
  "plantMeasurements__width2",
  "plantMeasurements__height",
  "height_sample",
  "Ambient_Humidity",
  "Ambient_Temperature",
  "Leaf_Angle",
  "Leaf_Temp_Differential",
  "LEF",
  "Light_Intensity..PAR.",
  "NPQt",
  "Phi2",
  "PhiNO",
  "PhiNPQ",
  "Relative_Chlorophyll",
  "thickness",
  "absorbance_420",
  # "absorbance_530",
  # "absorbance_605",
  # "absorbance_650",
  # "absorbance_730",
  # "absorbance_850",
  # "absorbance_880",
  # "absorbance_940",
  "angle",
  "contactless_temp",
  "ecs_initial",
  "ecs_max",
  "ecs_r_squared",
  "flatten_slope",
  "FmPrime",
  "FoPrime",
  "Fs",
  "FvP.FmP",
  "G",
  "gH.",
  "humidity",
  "humidity2",
  "light2_raw1",
  # "light2_raw2",
  # "light2_raw3",
  # "light6_raw1",
  # "light6_raw2",
  # "light6_raw3",
  "light_intensity",
  "light_intensity_raw",
  "MPF_rsquared",
  "MPF_slope",
  "NPQt_MPF",
  # "NPQt_noMPF",
  #"Phi2_MPF",
  "PhiNPQ_MPF",
  # "PhiNPQ_noMPF",
  "pitch",
  "pressure",
  "pressure2",
  "qL",
  "R",
  "ratio_MPF.noMPF_Phi2",
  "ratio_MPF.noMPF_PhiNO",
  "ratio_MPF.noMPF_PhiNPQ",
  "Rel_Chl_intensity",
  "RFd",
  "roll",
  "SPAD_420", #stand ins for other spad measurements
  "SPAD_420_intensity",
  # "SPAD_530",
  # "SPAD_530_intensity",
  # "SPAD_605",
  # "SPAD_605_intensity",
  # "SPAD_650",
  # "SPAD_650_intensity",
  # "SPAD_730",
  # "SPAD_730_intensity",
  # "SPAD_850",
  # "SPAD_850_intensity",
  # "SPAD_880",
  # "SPAD_880_intensity",
  # "spad_raw1",
  # "spad_raw2",
  # "spad_raw3",
  "temperature",
  "temperature2",
  "TimeofDay",
  "Latitude",
  "Longitude",
  "waterRetention",
  "toughness",
  "elev_m",
  "slope_perc",
  #"Startday0to365", #not calculated
  #"canopyCover",
  "treeRich",
  "shrubRich",
  "deadDown",
  # "numSamples", #not correlated with much of anything
 # "taxon.x",
  "shannons",
  "simpsons" ,                
 "shannonsISD" , #Note that the info content is very similar regardless of whether it was normalized or not
 "simpsonsISD"   
)
meta1$waterRetention <- as.numeric(meta1$waterRetention)

meta <- meta1[,names(meta1) %in% predictors]
meta <- data.frame(scale(meta))

#impute data to fill the matrix
table(is.na(meta))

meta_impute <- rfImpute(x = meta[,1:67], y = meta$shannons,
         iter=5, ntree=300)

samps <- seq(1,dim(meta_impute)[1],1)

#to repeat
set.seed(666)
dummy_taxon <- fastDummies::dummy_cols(meta1$taxon_final)

rf_out <- list()
for(i in 1:100){
  indices <- sample(samps, replace = F, size = 0.15*length(samps))
  train <- meta_impute[indices,
                       2:length(meta_impute)]
  
  rf <- randomForest(x = data.frame(train, 
                                    dummy_taxon[indices, 2:65]),
                     y = meta$shannons[indices],
                     importance = T,						#save variable importance metrics
                     ntree=4000,							#number of trees to make
                     nodesize = 5,						#this is the number of rows of data in the terminal node, defaults at 5 for regression
                     mtry = length(meta_impute)/3				#this is m, the number of predictors to consider at each split. 
                     #The default for regression is p/3 (which I just spelled out here for clarity. p is the number of predictors by convention)
  )
  #Extract % var explained in this cumbersome way. 
  rf$var_explained <- 100*rf$rsq[length(rf$rsq)]
  
  out <- predict(object = rf, 
                 newdata = data.frame(
                   meta_impute[-indices,
                               2:length(meta_impute)], 
                   dummy_taxon[-indices, 2:65]),
                 type = "response")
  
  rf$rmse <- Metrics::rmse(out, meta$shannons[-indices])
  rf$cor <- cor.test(out, meta$shannons[-indices])
  #plot(out, meta$shannons[-indices])
  
  rf_out[[i]] <- rf
}

# plot(out, meta$shannons) 
# 
# partialPlot(rf,
#             pred.data = train,					#training data
#             x.var = "height_sample")	#variable name for which u want to see the partial plot, 
# #has to be the character name NOT the actual object (eg., train$area will not work)!
# 
# cor.test(meta$shannons, meta$height_sample)
# plot(meta$shannons, meta$height_sample)

#CONSIDERATIONS
#1. scaling and centering features
#2. train/test split

reg <- lm(meta1$shannons ~ . , data = data.frame(train[,1:(length(train)-1)], meta1$taxon.y))
summary(reg)
cor.test(meta$shannons, meta$Light_Intensity..PAR.)
dim(train)

######
# 16S #
#####

rm(list=ls())
meta1 <- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv", 
                  stringsAsFactors = F, header = T)

predictors<-c(
  "area_cm2",
  "mass_g",
  "mass_extracted_g",
  #  "leaves_extracted",
  "circumStem",
  "width",
  "plantMeasurements__width2",
  "plantMeasurements__height",
  "height_sample",
  "Ambient_Humidity",
  "Ambient_Temperature",
  "Leaf_Angle",
  "Leaf_Temp_Differential",
  "LEF",
  "Light_Intensity..PAR.",
  "NPQt",
  "Phi2",
  "PhiNO",
  "PhiNPQ",
  "Relative_Chlorophyll",
  "thickness",
  "absorbance_420",
  # "absorbance_530",
  # "absorbance_605",
  # "absorbance_650",
  # "absorbance_730",
  # "absorbance_850",
  # "absorbance_880",
  # "absorbance_940",
  "angle",
  "contactless_temp",
  "ecs_initial",
  "ecs_max",
  "ecs_r_squared",
  "flatten_slope",
  "FmPrime",
  "FoPrime",
  "Fs",
  "FvP.FmP",
  "G",
  "gH.",
  "humidity",
  "humidity2",
  "light2_raw1",
  # "light2_raw2",
  # "light2_raw3",
  # "light6_raw1",
  # "light6_raw2",
  # "light6_raw3",
  "light_intensity",
  "light_intensity_raw",
  "MPF_rsquared",
  "MPF_slope",
  "NPQt_MPF",
  # "NPQt_noMPF",
  #"Phi2_MPF",
  "PhiNPQ_MPF",
  # "PhiNPQ_noMPF",
  "pitch",
  "pressure",
  "pressure2",
  "qL",
  "R",
  "ratio_MPF.noMPF_Phi2",
  "ratio_MPF.noMPF_PhiNO",
  "ratio_MPF.noMPF_PhiNPQ",
  "Rel_Chl_intensity",
  "RFd",
  "roll",
  "SPAD_420", #stand ins for other spad measurements
  "SPAD_420_intensity",
  # "SPAD_530",
  # "SPAD_530_intensity",
  # "SPAD_605",
  # "SPAD_605_intensity",
  # "SPAD_650",
  # "SPAD_650_intensity",
  # "SPAD_730",
  # "SPAD_730_intensity",
  # "SPAD_850",
  # "SPAD_850_intensity",
  # "SPAD_880",
  # "SPAD_880_intensity",
  # "spad_raw1",
  # "spad_raw2",
  # "spad_raw3",
  "temperature",
  "temperature2",
  "TimeofDay",
  "Latitude",
  "Longitude",
  "waterRetention",
  "toughness",
  "elev_m",
  "slope_perc",
  #"Startday0to365", #not calculated
  #"canopyCover",
  "treeRich",
  "shrubRich",
  "deadDown",
  # "numSamples", #not correlated with much of anything
  # "taxon.x",
  "shannons",
  "simpsons" ,                
  "shannonsISD" , #Note that the info content is very similar regardless of whether it was normalized or not
  "simpsonsISD"   
)
meta1$waterRetention <- as.numeric(meta1$waterRetention)

meta <- meta1[,names(meta1) %in% predictors]

#impute data to fill the matrix
table(is.na(meta))

library(randomForest)
meta_impute <- rfImpute(x = meta[,1:67], y = meta$shannons,
                        iter=5, ntree=300)

train <- meta_impute[,2:length(meta_impute)]

set.seed(666)
rf <- randomForest(x = data.frame(train, 
                                  fastDummiesmeta1$taxon_final),
                   y = meta$shannons,
                   importance = T,						#save variable importance metrics
                   ntree=4000,							#number of trees to make
                   nodesize = 5,						#this is the number of rows of data in the terminal node, defaults at 5 for regression
                   mtry = length(meta_impute)/3				#this is m, the number of predictors to consider at each split. 
                   #The default for regression is p/3 (which I just spelled out here for clarity. p is the number of predictors by convention)
)

varImpPlot(rf, sort=T)	

reg <- lm(meta1$shannons ~ . , data = data.frame(train[,1:(length(train)-1)], meta1$taxon_final))
summary(reg)
dim(train)
