rm(list=ls())
#install.packages("corrplot")
library(corrplot)

meta <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv", 
                 stringsAsFactors = F, header = T)

names(meta)

meta <- meta[,names(meta) %in% c("area_cm2"                                                 
                             , "mass_extracted_g"                                         
                             , "leaves_extracted"                                         
                             , "circumStem"                                               
                             , "height_sample"                                            
                             , "Ambient_Humidity"                                         
                             , "Ambient_Temperature"                                      
                             , "Leaf_Temp_Differential"                                   
                             , "LEF"                                                      
                             , "Light_Intensity..PAR."                                    
                             , "Phi2"                                                     
                             , "PhiNO"                                                    
                             , "PhiNPQ"                                                   
                             , "Relative_Chlorophyll"                                     
                             , "thickness"                                                
                             , "absorbance_420"                                           
                             , "absorbance_940"                                           
                             , "B"                                                        
                             , "contactless_temp"                                         
                             , "ecs_initial"                                              
                             , "ecs_max"                                                  
                             , "FmPrime"                                                  
                             , "FoPrime"                                                  
                             , "Fs"                                                       
                             , "FvP.FmP"                                                  
                             , "G"                                                        
                             , "gH."                                                      
                             #, "NPQt_MPF"                                                 
                             , "pressure"                                                 
                             , "qL"                                                       
                             , "R"                                                        
                             , "Rel_Chl_intensity"                                        
                             , "RFd"                                                      
                             , "SPAD_420"                                                 
                             , "SPAD_420_intensity"                                       
                             , "TimeofDay"                                                
                             , "lat"                                                 
                             , "long"                                                
                             , "waterRetention"                                           
                             , "toughness"                                                
                             , "elev_m"                                                   
                             , "slope_perc"                                               
                             , "treeRich"                                                 
                             , "shrubRich"                                                
                             , "deadDown"                                                 
                             , "precip_april_in.x"                                        
                             , "densitometer"                                             
                             , "shannons_flora"
                             # , "shannonsISD"
                             , "julianDate"                                               
                             , "mean_temp_april.y"                                        
                             , "plant_vol"                                                
                             , "sla"                                                      
                           )]


meta$leaves_extracted[meta$leaves_extracted == "partial"] <- 0.5

meta$waterRetention <- as.numeric(meta$waterRetention)
meta$leaves_extracted <- as.numeric(meta$leaves_extracted)

for(i in names(meta)){
  if(length(unique(meta[, names(meta) == i])) < 2){
    print(i)
  }
}


names(meta) <- c(
  "Area cm^2",
  "Mass extracted g",
  "Leaves extracted",
  "Circumference (trees only)",
  "Height sample taken",
  "Ambient humidity",
  "Ambient temperature",
  "Leaf temp. differential",
  "Linear electron flow",
  "Light intensity (PAR)",
  "Phi2",
  "PhiNO",
  "PhiNPQ",
  "Rel. chlorophyll",
  "Leaf thickness",
  "Absorbance 420 nm",
  "Absorbance 940 nm",
  "Blueness",
  "Leaf temperature",
  "ECS initial",
  "ECS max.",
  "Fm prime",
  "Fo prime",
  "Fs",
  "FvP.FmP",
  "Greenness",
  "gH",
  "Pressure",
  "qL",
  "Redness",
  "Rel. chlorophyll intensity",
  "RFd",
  "SPAD 420 nm",
  "SPAD 420 intensity",
  "Time of the day",
  "Water retention",
  "Toughness",
  "Latitude",
  "Longitude",
  "Elevation",
  "Slope %",
  "Tree richness",
  "Shrub richness",
  "Dead and down trees",
  "Precip. in April",
  "Densitometer",
  "Shannon's flora",
  "Julian date",
  "Mean temperature in April",
  "Plant volume")


cormat <- cor(meta,
              use = "pairwise.complete.obs")

pdf(width = 12, height = 12, file = "./visuals/correlogram_kitchensink.pdf")
  corrplot(cormat, method="circle")
dev.off()



