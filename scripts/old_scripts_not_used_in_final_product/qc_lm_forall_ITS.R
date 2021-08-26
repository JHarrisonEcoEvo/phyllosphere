rm(list=ls())
taxa <- read.csv("./processedData/otuTables/its97_cluster_otus_algorithm_otutableCLEAN",
                  stringsAsFactors = F)
taxa <- read.csv("./processedData//ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv")

zotus <- taxa$OTUID

#Need to get metadata into the same order as the modeled data
sample <- gsub("_ITS_[1,2]_\\w+$", "",names(taxa))
sample <- gsub("ITS_[ATCG]*_[ATCG]*_", "",sample)
sample <- toupper(sample)
sample <- gsub("(\\d+)(E)", "\\1_\\2",sample)

taxa <- data.frame(t(taxa[2:length(taxa[,1]),2:length(taxa)]))
taxa[,1:3]
names(taxa) <- zotus[2:length(zotus)]
taxa$sample <- sample[2:length(sample)]

X <- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv",
                    stringsAsFactors = F)

#######################
# Feature engineering #
######################

X <- X[X$substrate == "plant",]

X <- X[,names(X) != "shannonsISD"] #remove microbial diversity

#Make a leaf density variable
X$sla = X$mass_extracted_g / X$area_cm2

X$sla[is.infinite(X$sla)] <- NA

#Select only the features that I care about
#only a single one of a group of highly correlated features was chosen
#e.g., absorbance measurements were redundant in some cases.
#See the correlagram.R script
X <- X[,names(X) %in% c(
  "absorbance_420",
  "absorbance_940", 
  "area_cm2",
  "Ambient_Humidity",
  "Ambient_Temperature",
  "B",
  #"canopyCover" , #Not sure where my measurement of this went. Maybe I didn't mak them?
  "circumStem",
  "compartment", #check this out
  "contactless_temp",
  "deadDown",
  "densitometer",
  "ecs_initial",
  "ecs_max",
  "elev_m",
  "FmPrime",
  "FoPrime",
  "Fs",
  "FvP.FmP",
  "G",
  "gH.",
  "habit", 
  "height_sample",
  "julianDate",
  "lat",
  "Leaf_Temp_Differential",
  "LEF",
  "leaves_extracted",
  "Light_Intensity..PAR.",
  "long",
  "mass_extracted_g",
  "mean_temp_april.y" ,          
  "MEM1.y",
  "MEM2.y",
  "NPQt_MPF",
  "phenology",
  "Phi2",
  "PhiNO",
  "PhiNPQ",
  "plant_vol",
  "precip_april_in.x",
  "pressure",
  "qL",
  "R",
  "Rel_Chl_intensity",
  "Relative_Chlorophyll",
  "RFd",
  "sample",
  "shannons_flora",
  "shrubRich",
  #"sidePlantSampled",
  "sla",
  "slope_perc",
  "SPAD_420_intensity",
  "SPAD_420",
  "thickness",
  "TimeofDay",
  "taxon_final",
  "toughness",
  "treeRich",
  "waterRetention"
)]

#Convert partial leaf to 0.5
X$leaves_extracted <- (as.character(X$leaves_extracted))
X$leaves_extracted[X$leaves_extracted == "partial"] <- 0.5
X$leaves_extracted <- (as.numeric(X$leaves_extracted))

X$waterRetention <- (as.character(X$waterRetention))
indices <- which(X$waterRetention == "90+")
X$waterRetention[indices] <- 90
X$waterRetention <- (as.numeric(X$waterRetention))

#Do dummy encoding of all host taxa and other categoricals
#Writing out all the categorical features for clarity
categoricals <- c(
  "habit",
  "compartment",
  "taxon_final",
  "phenology"
)
X <- fastDummies::dummy_cols(.data = X, select_columns = categoricals)

table(taxa$sample %in%  X$sample )

merged_dat <- merge(X, taxa, by.x = "sample", by.y = "sample", all.y = T)
dim(merged_dat)


###########

#normalize
merged_dat[,130:length(merged_dat)] <- merged_dat[,130:length(merged_dat)] / merged_dat$ISD 

#do regression
features <- merged_dat[,names(merged_dat) %in%
                         c(#"response_taxon"                                           
                           #, "sample"                                                   
                            "area_cm2"                                                 
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
                           , "NPQt_MPF"                                                 
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
                           , "julianDate"                                               
                           , "mean_temp_april.y"                                        
                           , "plant_vol"                                                
                           , "sla"                                                      
                           , "habit_forb"                                               
                           , "habit_graminoid"                                          
                           , "habit_shrub"                                              
                           , "habit_tree"                                               
                           , "compartment_EN"                                           
                           , "compartment_EP"                                           
                           , "taxon_final_Abies.concolor"                               
                           , "taxon_final_Abies.grandis"                                
                           , "taxon_final_Antennaria.media"                             
                           , "taxon_final_Aquilegia.caerula"                            
                           , "taxon_final_Arnica.cordifolia"                            
                           , "taxon_final_Artemisia.tridentata"                         
                           , "taxon_final_Astragalus.alpinus"                           
                           , "taxon_final_Astragalus.kentrophyta"                       
                           , "taxon_final_Astragalus.miser"                             
                           , "taxon_final_Bistorta.bistortoides"                        
                           , "taxon_final_Carex.paysonis"                               
                           , "taxon_final_Carex.scopulorum.var..bracteosa"              
                           , "taxon_final_Chamaenerion.angustifolium.var..angustifolium"
                           , "taxon_final_Delphinium.occidentale"                       
                           , "taxon_final_Erigeron.glacialis"                           
                           , "taxon_final_Eriogonum.umbellatum"                         
                           , "taxon_final_Eucephalus.elegans"                           
                           , "taxon_final_Eucephalus.engelmannii"                       
                           , "taxon_final_Fragaria.virginiana"                          
                           , "taxon_final_Frasera.speciosa"                             
                           , "taxon_final_Geranium.viscossimum.var..viscosissimum"      
                           , "taxon_final_Helianthella.uniflora"                        
                           , "taxon_final_Heracleum.maximum"                            
                           , "taxon_final_Juncus.balticus."                             
                           , "taxon_final_Juncus.parryi"                                
                           , "taxon_final_Juncus.sp."                                   
                           , "taxon_final_Juniperus.communis"                           
                           , "taxon_final_Ligusticum.filicinum"                         
                           , "taxon_final_Lupinus.argenteus"                            
                           , "taxon_final_Mahonia.repens"                               
                           , "taxon_final_Mertensiana.ciliata.var..ciliata"             
                           , "taxon_final_Minuartia.obtusiloba"                         
                           , "taxon_final_Osmorhiza.depauperata"                        
                           , "taxon_final_Oxyria.digyna"                                
                           , "taxon_final_Paxistima.myrsinites"                         
                           , "taxon_final_Picea.engelmannii"                            
                           , "taxon_final_Pinus.albicaulis"                             
                           , "taxon_final_Pinus.contorta"                               
                           , "taxon_final_Poa.pratensis"                                
                           , "taxon_final_Poa.wheeleri"                                 
                           , "taxon_final_Polemonium.viscosum"                          
                           , "taxon_final_Populus.angustifolium"                        
                           , "taxon_final_Populus.tremulloides"                         
                           , "taxon_final_Potentilla.diversifolia"                      
                           , "taxon_final_Potentilla.fruticosa"                         
                           , "taxon_final_Potentilla.pulcherrima"                       
                           , "taxon_final_Primula.parryi"                               
                           , "taxon_final_Pseudotsuga.menziesii"                        
                           , "taxon_final_Ribes.montigenum"                             
                           , "taxon_final_Salix.glauca.var..villosa"                    
                           , "taxon_final_Salix.reticulata.var..nana"                   
                           , "taxon_final_Sedum.lanceolatum"                            
                           , "taxon_final_Symphoricarpos.albus"                         
                           ,"taxon_final_Symphyotrichum.foliaceum.var..apricum"        
                           , "taxon_final_Thalictrum.occidentale"                       
                           , "taxon_final_Trifolium.parryi.var..montanense"             
                           , "taxon_final_Unknown.fir"                                  
                           , "taxon_final_Unknown.pine"                                 
                           , "taxon_final_Unknown.Spruce"                               
                           , "taxon_final_Vaccinium.membranaceum"                       
                           , "taxon_final_Vaccinium.scoparium"                          
                           , "taxon_final_Wyethia.amplexicaulis"                        
                           , "phenology_flowering"                                      
                           , "phenology_fruiting"                                       
                           , "phenology_vegetative"                                     
                         )]

#median imputation
features <- data.frame(lapply(features,function(x) {
  if(is.numeric(x)) ifelse(is.na(x),median(x,na.rm=T),x) else x}))

merged_dat <- data.frame(lapply(merged_dat,function(x) {
  if(is.numeric(x)) ifelse(is.infinite(x),median(x,na.rm=T),x) else x}))

taxon <- NA
rsq <- NA
for(i in 130:(length(merged_dat))){
  tormv <- which(is.na(merged_dat[,i]))
  response <- na.omit(merged_dat[,i])

  reg <- lm(response~ ., data = features[-tormv, ])
  taxon[i] <- names(merged_dat)[i]
  rsq[i] <- summary(reg)$r.squared
}
summary(rsq)




