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
library(ranger)

set.seed(666)

#CHOOSE INPUT 
#X<- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv")
X<- read.csv("./processedData/16smetadat_wrangled_for_post_modeling_analysis.csv")

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
         c("shannonsISD",
           "simpsonsISD",
           "sample"
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
           , "taxon_final_Abies concolor"                               
           , "taxon_final_Abies grandis"                                
           , "taxon_final_Antennaria media"                             
           , "taxon_final_Aquilegia caerula"                            
           , "taxon_final_Arnica cordifolia"                            
           , "taxon_final_Artemisia tridentata"                         
           , "taxon_final_Astragalus alpinus"                           
           , "taxon_final_Astragalus kentrophyta"                       
           , "taxon_final_Astragalus miser"                             
           , "taxon_final_Bistorta bistortoides"                        
           , "taxon_final_Carex paysonis"                               
           , "taxon_final_Carex scopulorum var. bracteosa"              
           , "taxon_final_Chamaenerion angustifolium var. angustifolium"
           , "taxon_final_Delphinium occidentale"                       
           , "taxon_final_Erigeron glacialis"                           
           , "taxon_final_Eriogonum umbellatum"                         
           , "taxon_final_Eucephalus elegans"                           
           , "taxon_final_Eucephalus engelmannii"                       
           , "taxon_final_Fragaria virginiana"                          
           , "taxon_final_Frasera speciosa"                             
           , "taxon_final_Geranium viscossimum var. viscosissimum"      
           , "taxon_final_Helianthella uniflora"                        
           , "taxon_final_Heracleum maximum"                            
           , "taxon_final_Juncus balticus."                             
           , "taxon_final_Juncus parryi"                                
           , "taxon_final_Juncus sp."                                   
           , "taxon_final_Juniperus communis"                           
           , "taxon_final_Ligusticum filicinum"                         
           , "taxon_final_Lupinus argenteus"                            
           , "taxon_final_Mahonia repens"                               
           , "taxon_final_Mertensiana ciliata var. ciliata"             
           , "taxon_final_Minuartia obtusiloba"                         
           , "taxon_final_Osmorhiza depauperata"                        
           , "taxon_final_Oxyria digyna"                                
           , "taxon_final_Paxistima myrsinites"                         
           , "taxon_final_Picea engelmannii"                            
           , "taxon_final_Pinus albicaulis"                             
           , "taxon_final_Pinus contorta"                               
           , "taxon_final_Poa pratensis"                                
           , "taxon_final_Poa wheeleri"                                 
           , "taxon_final_Polemonium viscosum"                          
           , "taxon_final_Populus angustifolium"                        
           , "taxon_final_Populus tremulloides"                         
           , "taxon_final_Potentilla diversifolia"                      
           , "taxon_final_Potentilla fruticosa"                         
           , "taxon_final_Potentilla pulcherrima"                       
           , "taxon_final_Primula parryi"                               
           , "taxon_final_Pseudotsuga menziesii"                        
           , "taxon_final_Ribes montigenum"                             
           , "taxon_final_Salix glauca var. villosa"                    
           , "taxon_final_Salix reticulata var. nana"                   
           , "taxon_final_Sedum lanceolatum"                            
           , "taxon_final_Symphoricarpos albus"                         
           ,"taxon_final_Symphyotrichum foliaceum var. apricum"        
           , "taxon_final_Thalictrum occidentale"                       
           , "taxon_final_Trifolium parryi var. montanense"             
           , "taxon_final_Unknown fir"                                  
           , "taxon_final_Unknown pine"                                 
           , "taxon_final_Unknown Spruce"                               
           , "taxon_final_Vaccinium membranaceum"                       
           , "taxon_final_Vaccinium scoparium"                          
           , "taxon_final_Wyethia amplexicaulis"                        
           , "phenology_flowering"                                      
           , "phenology_fruiting"                                       
           , "phenology_vegetative"                                     
         )]

#scale and center (now doing this in pipeline to avoid leakage)
#X[,c(1:46,48:52, 55)] <- scale(X[,c(1:46,48:52, 55)], scale = T, center = T)

X <- X[, -which(names(X) == "simpsonsISD")]
X <- X[, -which(names(X) == "sample")]
#X$sample <- as.character(X$sample)

#Convert to numeric (makes it easier when doing imputing)
for(i in 1:length(X)){
  X[,i] <- as.numeric(X[,i])  
}

#########################################
#Use mlr to define a pipeline, task and a learner#
########################################
names(X) <- gsub("\\s", "_", names(X))

#Tasks are a way to organize data and the learner is the model 
#see: https://mlr3book.mlr-org.com/tasks.html
phyllo_task = as_task_regr(X
                           , target = "shannonsISD", id = "phyllo_data")

#Define roles for all features
#phyllo_task$col_roles$name <- "sample" 

phyllo_task$col_roles$features <-  names(X)[!(names(X) %in% c("sample", "shannonsISD"))]


#Build pipeline for scaling and imputation
imp_missind <- po("missind")

#Note that imputation of numeric will NOT work with integer class features. 
#Vice versa doesn't work either
imp_num <- po("imputehist", param_vals = list(affect_columns = selector_type("numeric")))

#simple pipe
#graph <-  po("imputehist", param_vals = list(affect_columns = selector_type("numeric"))) %>>% 
# po("scale") #pass to learner
#graph$plot()


###########
# XGBoost #
###########

#make a graph learner
graph <-  po("imputehist", param_vals = list(affect_columns = selector_type("numeric"))) %>>% 
  po("scale") %>>%
  po( lrn("regr.xgboost", objective = "reg:squarederror"))

g1 <- GraphLearner$new(graph)

##########
# tuning #
##########
#see https://mlr3gallery.mlr-org.com/posts/2021-03-10-practical-tuning-series-tune-a-preprocessing-pipeline/

#For xgboost, 
#define hyperparameters to check
#Note the prefix of regr.xgboost which gets added since I am passing in a complicated learner
#that has parameters for the imputer too. 
#suse as.data.table(g1$param_set) to see what parameters to use. 
params <- ps(
  regr.xgboost.eta = p_dbl(lower = 0.1, upper = 0.3),
  regr.xgboost.nrounds = p_int(lower = 1, upper = 16, tags = "budget"),
  regr.xgboost.booster = p_fct(levels = c("gbtree", "gblinear", "dart")),
  regr.xgboost.gamma = p_int(lower = 0, upper = 15),
  regr.xgboost.max_depth = p_int(lower = 4, upper = 15),
  regr.xgboost.min_child_weight = p_int(lower = 1, upper = 10),
  regr.xgboost.subsample = p_dbl(lower = 0.5, upper = 0.8),
  regr.xgboost.colsample_bytree = p_dbl(lower = 0.5, upper = 1)
)

###################################################
# Nested resampling to estimate model performance #
###################################################

at = AutoTuner$new(
  learner = g1,
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
#rr$score()

#This is the unbiased estimate of model performance
rr$aggregate( measure = msr("regr.rsq")) #can swap out different measures
rr$aggregate(measure = msr("regr.mse") ) #can swap out different measures

##################################################################
# try with random forest
##################################################################

graph <-  po("imputehist", param_vals = list(affect_columns = selector_type("numeric"))) %>>% 
  po("scale") %>>%
  po( lrn("regr.ranger", importance = "permutation"))

g1 <- GraphLearner$new(graph)

##########
# tuning #
##########

#For ranger (randomforest)
params <- ps(
  regr.ranger.num.trees = p_int(lower = 50, upper = 400, tags = "budget"), #number of trees in the dark forest
  regr.ranger.sample.fraction = p_dbl(lower = 0.01, upper = 0.3), #The sample.fraction parameter specifies the fraction of observations to be used in each tree. Smaller fractions lead to greater diversity, and thus less correlated trees which often is desirable.
  regr.ranger.splitrule = p_fct(levels = c("extratrees", "variance")),
  regr.ranger.min.node.size = p_int(lower = 1, upper = 25)
)

###################################################
# Nested resampling to estimate model performance #
###################################################

at <- AutoTuner$new(
  learner = g1,
  resampling = rsmp("cv", folds = 3),
  measure = msr("regr.rsq"),
  terminator = trm("none"),
  tuner = tnr("hyperband", eta = 3),
  search_space = params)

#pass back to resampling for nested resampling
outer_resampling <- rsmp("cv", folds = 10)
rr <- resample(task = phyllo_task, 
               learner = at, 
               outer_resampling, 
               store_models = TRUE)

outer_resampling$instantiate(phyllo_task)

extract_inner_tuning_results(rr)

# #outer tuning results
rr$score()

#This is the unbiased estimate of model performance
rr$aggregate(measure = msr("regr.rsq")) #can swap out different measures
rr$aggregate(measure = msr("regr.mse") ) #can swap out different measures
