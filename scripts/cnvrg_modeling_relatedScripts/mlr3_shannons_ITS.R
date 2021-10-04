rm(list=ls())
library(mlr3) #cant get mlr3 to install on work computer
library(xgboost)
#install.packages("mlr3learners")
library(mlr3learners)
library(mlr3verse)
library(mlr3fselect)
library("mlr3hyperband")
library(fastDummies)
#library(mlr3pipelines)
library(ranger)
library(caret)

#library(doParallel)
#all_cores <- parallel::detectCores(logical = FALSE)
#registerDoParallel(cores = all_cores)

set.seed(666)

###########################
#Input and minor wrangling#
###########################

X<- read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv")

#######################
# Feature engineering #
######################

X <- X[X$substrate == "plant",]

#Make a leaf density variable
X$sla = X$mass_extracted_g / X$area_cm2

X$sla[is.infinite(X$sla)] <- NA

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
X$taxon_final <- gsub(" ", "", X$taxon_final)

X <- dummy_cols(.data = X, select_columns = categoricals)

#Get rid of stuff we don't need 
merged_dat <- X[,names(X) %in%
                           c(
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
                             ,"MEM1"
                             , "MEM2"
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
,"taxon_final_Abiesconcolor"
,"taxon_final_Abiesgrandis"
,"taxon_final_Antennariamedia"
,"taxon_final_Aquilegiacaerula"
,"taxon_final_Arnicacordifolia"
,"taxon_final_Artemisiatridentata"
,"taxon_final_Astragalusalpinus"
,"taxon_final_Astragaluskentrophyta"
,"taxon_final_Astragalusmiser"
,"taxon_final_Bistortabistortoides"
,"taxon_final_Carexpaysonis"
,"taxon_final_Carexscopulorumvar.bracteosa"
,"taxon_final_Chamaenerionangustifoliumvar.angustifolium"
,"taxon_final_Delphiniumoccidentale"
,"taxon_final_Erigeronglacialis"
,"taxon_final_Eriogonumumbellatum"
,"taxon_final_Eucephaluselegans"
,"taxon_final_Eucephalusengelmannii"
,"taxon_final_Fragariavirginiana"
,"taxon_final_Fraseraspeciosa"
,"taxon_final_Geraniumviscossimumvar.viscosissimum"
,"taxon_final_Helianthellauniflora"
,"taxon_final_Heracleummaximum"
,"taxon_final_Juncusbalticus."
,"taxon_final_Juncusparryi"
,"taxon_final_Juncussp."
,"taxon_final_Juniperuscommunis"
,"taxon_final_Ligusticumfilicinum"
,"taxon_final_Lupinusargenteus"
,"taxon_final_Mahoniarepens"
,"taxon_final_Mertensianaciliatavar.ciliata"
,"taxon_final_Minuartiaobtusiloba"
,"taxon_final_Osmorhizadepauperata"
,"taxon_final_Oxyriadigyna"
,"taxon_final_Paxistimamyrsinites"
,"taxon_final_Piceaengelmannii"
,"taxon_final_Pinusalbicaulis"
,"taxon_final_Pinuscontorta"
,"taxon_final_Poapratensis"
,"taxon_final_Poawheeleri"
,"taxon_final_Polemoniumviscosum"
,"taxon_final_Populusangustifolium"
,"taxon_final_Populustremulloides"
,"taxon_final_Potentilladiversifolia"
,"taxon_final_Potentillafruticosa"
,"taxon_final_Potentillapulcherrima"
,"taxon_final_Primulaparryi"
,"taxon_final_Pseudotsugamenziesii"
,"taxon_final_Ribesmontigenum"
,"taxon_final_Salixglaucavar.villosa"
,"taxon_final_Salixreticulatavar.nana"
,"taxon_final_Sedumlanceolatum"
,"taxon_final_Symphoricarposalbus"
,"taxon_final_Symphyotrichumfoliaceumvar.apricum"
,"taxon_final_Thalictrumoccidentale"
,"taxon_final_Trifoliumparryivar.montanense"
,"taxon_final_Unknownfir"
,"taxon_final_Unknownpine"
,"taxon_final_UnknownSpruce"
,"taxon_final_Vacciniummembranaceum"
,"taxon_final_Vacciniumscoparium"
,"taxon_final_Wyethiaamplexicaulis"
                             , "phenology_flowering"                                      
                             , "phenology_fruiting"                                       
                             , "phenology_vegetative"                                     
                           )]

#Convert to numeric (makes it easier when doing imputing)
for(i in 2:length(merged_dat)){
  merged_dat[,i] <- as.numeric(merged_dat[,i])  
}

#########################################
#Use mlr to define a pipeline, task and a learner#
########################################

#Make a stratum column that is just a dupe of the response
merged_dat$shannonsISDstratum <- ifelse(merged_dat$shannonsISD > 
                                          mean(merged_dat$shannonsISD), 1, 0)

names(merged_dat) <- gsub("\\s", "_", names(merged_dat))

#Tasks are a way to organize data and the learner is the model 
#see: https://mlr3book.mlr-org.com/tasks.html
phyllo_task = TaskRegr$new(backend = merged_dat
                           , target = "shannonsISD",
                           id = "phyllo_data")

#Define roles for all features
phyllo_task$set_col_roles("sample", roles = "name")
phyllo_task$set_col_roles("shannonsISDstratum", roles = "stratum")

phyllo_task$col_roles$feature <-  names(merged_dat)[!(names(merged_dat) %in% c("sample", "shannonsISD","shannonsISDstratum"))]


#Build pipeline for scaling and imputation
imp_missind <- po("missind")

#Note that imputation of numeric will NOT work with integer class features. 
#Vice versa doesn't work either
imp_num <- po("imputehist", param_vals = list(affect_columns = selector_type("numeric")))

######
# QC #
######
# Simulate some data to see if the model works properly with a better feature. 
# 
# merged_dat$simmed <- merged_dat$shannonsISD * rnorm(mean = 0, sd = 100, n = length(merged_dat$response_taxon))
# merged_dat$simmed <- merged_dat$shannonsISD * rnorm(mean = 0, sd = 10, n = length(merged_dat$response_taxon))

##########
# RF #
##########

graph <-  po("imputehist", param_vals = list(affect_columns = selector_type("numeric"))) %>>% 
  po("scale", param_vals = list(scale = T, center = T)) %>>%
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
  resampling = rsmp("cv", folds = 4),
  measure = msr("regr.rsq"),
  terminator = trm("none"),
  tuner = tnr("hyperband", eta = 3),
  search_space = params, 
  store_models = TRUE)

#pass back to resampling for nested resampling
outer_resampling <- rsmp("cv", folds = 3)
rr <- resample(task = phyllo_task, 
               learner = at, 
               outer_resampling, 
               store_models = TRUE)

outer_resampling$instantiate(phyllo_task)

extract_inner_tuning_results(rr)

# #outer tuning results
rr$score()

#This is the unbiased estimate of model performance
rsq <- rr$aggregate(measure = msr("regr.rsq")) #can swap out different measures
mse <- rr$aggregate(measure = msr("regr.mse") ) #can swap out different measures


##################################################################
# Extract important features
##################################################################

#Extremely clunky to get the variable importance. This took me a lot of digging to even find
#since you can't query the object effectively with str()

tained_at <- at$train(phyllo_task)

var.imp <- data.frame(tained_at$model$learner$model$regr.ranger$model$variable.importance)

out <- data.frame(matrix(nrow = 1, ncol = 1))
out$rsq_nested_resampling <- rsq
out$mse_nested_resampling <- mse

#rank by importance
viewdf <- data.frame(
  row.names(var.imp)[order(var.imp$tained_at.model.learner.model.regr.ranger.model.variable.importance)],
  var.imp[order(var.imp$tained_at.model.learner.model.regr.ranger.model.variable.importance),1] 
)
names(viewdf) <- ""
viewdf

write.csv(var.imp, file = paste("variableImportanceshannonsISD_ITS.csv"), row.names = T)

##########################################
# Do over AFTER reducing the feature set#
##########################################

#####################
# Feature selection #
#####################
#Note, I am going to do this prior to tuning, otherwise I would want to do 
#this after tuning, then retune, which seems somewhat in the realm of
#diminishing returns and would increase runtime.

#https://mlr3book.mlr-org.com/fs.html

feat_selector = FSelectInstanceSingleCrit$new(
  task = phyllo_task,
  learner = g1,
  resampling = rsmp("cv", folds = 3),
  measure = msr("regr.rsq"),
  terminator = trm("evals", n_evals = 200)
)
#fselector = fs("rfe")
fselector = fs("sequential")

# reduce logging output
lgr::get_logger("bbotk")$set_threshold("warn")

fselector$optimize(feat_selector)
#See what we are left with
feat_selector$result_feature_set

#Select the best features
if(length(feat_selector$result_feature_set) > 1){
  phyllo_task$select(feat_selector$result_feature_set)
}


###################################################
# Nested resampling to estimate model performance #
###################################################

at_reduced <- AutoTuner$new(
  learner = g1,
  resampling = rsmp("cv", folds = 4),
  measure = msr("regr.rsq"),
  terminator = trm("none"),
  tuner = tnr("hyperband", eta = 3),
  search_space = params)

#pass back to resampling for nested resampling
outer_resampling <- rsmp("cv", folds = 3)
rr <- resample(task = phyllo_task, 
               learner = at_reduced, 
               outer_resampling, 
               store_models = TRUE)

outer_resampling$instantiate(phyllo_task)

#This is the unbiased estimate of model performance
rsq_reduced <- rr$aggregate(measure = msr("regr.rsq")) #can swap out different measures
mse_reduced <- rr$aggregate(measure = msr("regr.mse") ) #can swap out different measures

out$rsq_reduced<- rsq_reduced
out$mse_reduced <- mse_reduced

tained_at_reduced <- at$train(phyllo_task)

var.impR <- data.frame(tained_at_reduced$model$learner$model$regr.ranger$model$variable.importance)

write.csv(var.impR, file = paste("variableImportanceReducedshannonsISD_ITS.csv"), row.names = T)

write.csv(out, file = paste("results_shannonosISD_ITS.csv"), row.names = F)

