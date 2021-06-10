#install.packages("corrplot")
library(corrplot)

meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", 
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
  "absorbance_530",
  "absorbance_605",
  "absorbance_650",
  "absorbance_730",
  "absorbance_850",
  "absorbance_880",
  "absorbance_940",
 # "angle",
  "contactless_temp",
  "ecs_initial",
  "ecs_max",
  "ecs_r_squared",
  #"flatten_slope",
  "FmPrime",
  "FoPrime",
  "Fs",
  "FvP.FmP",
  "G",
  "gH.",
  "humidity",
  # "humidity2",
  # "light2_raw1",
  # "light2_raw2",
  # "light2_raw3",
  # "light6_raw1",
  # "light6_raw2",
  # "light6_raw3",
  # "light_intensity",
  # "light_intensity_raw",
  # "MPF_rsquared",
  # "MPF_slope",
  "NPQt_MPF",
 # "NPQt_noMPF",
  #"Phi2_MPF",
  "PhiNPQ",
 # "PhiNPQ_noMPF",
  #"pitch",
  "pressure",
  #"pressure2", #not sure what this is, but dont need it as seems very similar to pressure. 
  "qL",
  "R",
  # "ratio_MPF.noMPF_Phi2", #Not really sure what the value of these are so am omitting
  # "ratio_MPF.noMPF_PhiNO",
  # "ratio_MPF.noMPF_PhiNPQ",
  "Rel_Chl_intensity",
  "RFd",
  #"roll",
  "SPAD_420", #stand ins for other spad measurements
   "SPAD_420_intensity",
  "SPAD_530",
  "SPAD_530_intensity",
  "SPAD_605",
  "SPAD_605_intensity",
  "SPAD_650",
  "SPAD_650_intensity",
  "SPAD_730",
  "SPAD_730_intensity",
  "SPAD_850",
  "SPAD_850_intensity",
  "SPAD_880",
  "SPAD_880_intensity",
  "spad_raw1",
  "spad_raw2",
  "spad_raw3",
  "temperature",
  # "temperature2",
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
"latitude"    ,                    
"altitude"    ,                    
"longitude"                       
,"densitometer.y"   
,"shannons_flora.y"                
, "MEM1.y"                          
, "MEM2.y"                          
, "julianDate"                      
# , "shannons"                        
# , "simpsons"                        
, "shannonsISD"
,"simpsonsISD"
)
meta$waterRetention <- as.numeric(meta$waterRetention)
cormat <- cor(meta[,names(meta) %in% predictors], use = "pairwise.complete.obs")
pdf(width = 12, height = 12, file = "./visuals/correlogram_kitchensink.pdf")
  corrplot(cormat, method="circle")
dev.off()



