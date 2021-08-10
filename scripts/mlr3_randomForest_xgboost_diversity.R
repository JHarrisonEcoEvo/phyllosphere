rm(list=ls())
library(mlr3) 
library(xgboost)
library(mlr3learners)
library(mlr3verse)
library(mlr3fselect)
library("mlr3hyperband")
library(fastDummies)
library(randomForest)
library(caret)


set.seed(666)

X<- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv")

X <- X[X$substrate == "plant",]

#Make a leaf density variable
X$sla = X$mass_extracted_g / X$area_cm2

X$sla[is.infinite(X$sla)] <- NA

#Convert partial leaf to 0.5
X$leaves_extracted <- as.character(X$leaves_extracted)
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
X <- dummy_cols(.data = X, select_columns = categoricals)

#remove our bookkeeping data and various unneeded fields
X <- X[,names(X) %in%
                           c("response_taxon"                                           
                             , "sample"                                                   
                             , "area_cm2"                                                 
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
                             , "shannonsISD"
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


#impute
merged_dat[,
           c(1,3:(length(merged_dat)-2))] <- randomForest::rfImpute(response_taxon ~ .,
                                                                    iter=5, 
                                                                    ntree=300,
                                                                    data = merged_dat[,
                                                                                      c(1,3:(length(merged_dat)-2))])

#scale and center
merged_dat[,3:54] <- scale(merged_dat[,3:54], scale = T, center = T)


#########################################
#Use mlr to define a task and a learner#
########################################

#Tasks are a way to organize data and the learner is th emodel 
#see: https://mlr3book.mlr-org.com/tasks.html
phyllo_task = as_task_regr(merged_dat
                           , target = "response_taxon", id = "phyllo_data")

#Define roles for all features
phyllo_task$col_roles$feature <- names(merged_dat)[which(names(merged_dat)=="area_cm2") : which(names(merged_dat)=="phenology_vegetative")]

#QC only
#phyllo_task$col_roles$feature <- c(names(merged_dat)[which(names(merged_dat)=="area_cm2") : which(names(merged_dat)=="phenology_vegetative")], "simmed")

phyllo_task$col_roles$stratum <- "stratify"
phyllo_task$col_roles$name <- "sample" 


#XGBoost
# learner <- lrn("regr.xgboost", 
#                objective = "reg:squarederror")

#random forest
# learner <- lrn("regr.ranger", importance = "permutation") # should search for available CPUs

#lasso regression. I guess this does cv automatically
learner <- lrn("regr.cv_glmnet", alpha = 1, family = 'gaussian')

##########
# tuning #
##########

#For xgboost
#define hyperparameters to check via CV
# params <- ps(
#   eta = p_dbl(lower = 0.1, upper = 0.3),
#   nrounds = p_int(lower = 1, upper = 16, tags = "budget"),
#   booster = p_fct(levels = c("gbtree", "gblinear", "dart")),
#   gamma = p_int(lower = 0, upper = 15),
#   max_depth = p_int(lower = 4, upper = 15),
#   min_child_weight = p_int(lower = 1, upper = 10),
#   subsample = p_dbl(lower = 0.5, upper = 0.8),
#   colsample_bytree = p_dbl(lower = 0.5, upper = 1)
# )

###################################################
# Nested resampling to estimate model performance #
###################################################

at = AutoTuner$new(
  learner = learner,
  resampling = rsmp("cv", folds = 3),
  measure = msr("regr.rsq"),
  terminator = trm("none"),
  tuner = tnr("hyperband", eta = 3),
  search_space = params)

#pass back to resampling for nested resampling
outer_resampling = rsmp("cv", folds = 10)
rr = resample(task = phyllo_task, 
              learner = at, 
              outer_resampling, 
              store_models = TRUE)

outer_resampling$instantiate(phyllo_task)

extract_inner_tuning_results(rr)

# #outer tuning results
rr$score()

#This is the unbiased estimate of model performance
rsq <- rr$aggregate( measure = msr("regr.rsq")) #can swap out different measures
mse <- rr$aggregate(measure = msr("regr.mse") ) #can swap out different measures


##################################################################
# Extract important features
##################################################################
