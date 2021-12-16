#extract best predicted taxa. Figure out who they are, then see if they are predicted by sort of the same things. 

#Then extract their partial dependence plots to see if shape of responses are the same
#to certain influential predictors. This last part may be tricky. 

#Goal have a plot that shows the top best predicted taxa and their responses to certain features.
#overplotted lines, thus showing the idiosyncratic or similar response, as the case may be. 
#Do at landsacpe level for relative abundance and then do for the same taxa using absolute abundance
#Two panel plot probably....or maybe four panel if bacteria is in there too. OR maybe have an example
#from within a host, to show that the response is different when across host data are not considered...?
#The story may become that model structure is critical and leads to radically different inferences,
#Thus compounding the challenge of prediction. Or maybe responses will be sort of similar across models

rm(list=ls())
options(scipen = 99)
dat <- read.csv("./modelingResults/results_hella_noduds/all_hella_its_noduds.csv", stringsAsFactors = F)
head(dat)
dat <- dat[dat$taxon != "taxon",]
dim(dat)
dat$rsq_nested_resampling <- as.numeric(dat$rsq_nested_resampling)
summary(dat$rsq_nested_resampling)
length(dat[,1])
table(dat$rsq_nested_resampling > 0.1)
dat[which.max(dat$rsq_nested_resampling),]
#get the best predicted taxa
dat <- dat[dat$rsq_nested_resampling > 0.25,]
#See what they are
taxa <- read.csv("processedData/taxonomy_modeled_fungi_landscape.csv", stringsAsFactors = F)
taxa <- taxa[taxa$V1 %in% dat$taxon,]

#See what these blast as...doing outside of R
#THey all show up as fungi of various sorts
taxa[taxa$V4 == "d:Fungi",]


#Bring in variable importance metrics for only those taxa that were predicted
varimpdata <- list.files("modelingResults/varImp_landscape_hellinger/")
varimpdata <- varimpdata[grep("*ITS*", varimpdata)]

#find those files for the stuff that was predicted
varimpdata <- grep(paste(taxa$V1 , "FITS*", sep = "",collapse="|"), 
            varimpdata, value=TRUE)
#QC
length(varimpdata) == length(taxa$V1)

varimps <- lapply(paste("modelingResults/results_hella_noduds//noduds_", varimpdata, sep = ""), read.csv)

#summarize counts for important features across taxa
important <- vector()

for(i in 1:length(varimps)){
  important <- c(important,
                 as.character(varimps[[i]][head(rev(order(varimps[[i]][,2])), n = 10),1]))
}
table(important)
xtable::xtable(sort(table(important)))
#RESULT #1 Not a single feature was important for all taxa. 

length(varimps)


###################################################
# Extract Partial dependence plots for each taxon #
###################################################
#Strategy is to model each taxon, save the PDP to disk as a pdf, using identical plot dimensions. 
#Then overlay the lines in Affinity Designer. Will consider trying to save the equation of the line here and do
#plotting in this script, but it may not be worth the trouble. 

rm(list=ls()[ls() != "dat"])
library(mlr3)
library(mlr3learners)
library(mlr3verse)
library(mlr3fselect)
library("mlr3hyperband")
library(fastDummies)
#library(mlr3pipelines)
library(ranger)
library(caret)

set.seed(666)

features2examine <- 
  c("deadDown",
"MEM2",
"precip_april_in.x", 
"lat",
"mean_temp_april.y", 
"julianDate",
#"elev_m",
"shrubRich")


focal_taxon_line <- list()
k <- 1

for(z in 1:length(features2examine)){
  for (m in 1:length(dat$taxon)) {
  # # #Debugging stuff or running for a focal taxon
  taxa <-
    read.csv(
      "./processedData//otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG"
    )
  possibles <-
    read.csv("./processedData/ITS_taxa_to_model_via_randomforest.csv")
  focal_taxa <- dat$taxon
  X <-
    read.csv("./processedData/ITSmetadat_wrangled_for_post_modeling_analysis.csv")
  isdNorm_tf <- F
  focal_taxon <- dat$taxon[m]
  
  #######################
  # Feature engineering #
  ######################
  
  taxa$sample <- gsub("X", "", taxa$sample)
  
  taxa[, 3:length(taxa)] <- taxa[, 3:length(taxa)] - 1
  
  if (isdNorm_tf == "T") {
    print("dividing by the isd")
    taxa[, 3:length(taxa)] <- taxa[, 3:length(taxa)] / taxa$ISD
  } else{
    print("not dividing by the isd")
    taxa[, 3:length(taxa)]  <-
      vegan::decostand(taxa[, 3:(length(taxa) - 2)],
                       method = "hellinger")
  }
  
  X <- X[X$substrate == "plant", ]
  
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
  categoricals <- c("habit",
                    "compartment",
                    "taxon_final",
                    "phenology")
  #remove spaces in taxon_final for ease of output file manipulation.
  X$taxon_final <- gsub(" ", "", X$taxon_final)
  
  X <- dummy_cols(.data = X, select_columns = categoricals)
  
  #the missing stuff were things I removed due to poor sequencing
  table(X$sample %in% taxa$sample)
  table(taxa$sample %in%  X$sample)
  
  merged_dat <-
    merge(X,
          taxa,
          by.x = "sample",
          by.y = "sample",
          all.y = T)
  
  #For convenience make response variable
  response_taxon <- merged_dat[, names(merged_dat) == focal_taxon]
  
  #Get rid of stuff we don't need in merged_dat
  merged_dat <- merged_dat[, names(merged_dat) %in%
                             c(
                               as.character(focal_taxon),
                               "sample",
                               "area_cm2"
                               ,
                               "mass_extracted_g"
                               ,
                               "leaves_extracted"
                               ,
                               "circumStem"
                               ,
                               "height_sample"
                               ,
                               "Ambient_Humidity"
                               ,
                               "Ambient_Temperature"
                               ,
                               "Leaf_Temp_Differential"
                               ,
                               "LEF"
                               ,
                               "Light_Intensity..PAR."
                               # , "Phi2"
                               ,
                               "PhiNO"
                               ,
                               "PhiNPQ"
                               ,
                               "Relative_Chlorophyll"
                               ,
                               "thickness"
                               #, "absorbance_420"
                               ,
                               "absorbance_940"
                               # , "B"
                               #, "contactless_temp"
                               ,
                               "ecs_initial"
                               #, "ecs_max"
                               #, "FmPrime"
                               #, "FoPrime"
                               ,
                               "Fs"
                               #, "FvP.FmP"
                               #, "G"
                               ,
                               "gH."
                               #, "NPQt_MPF"
                               #, "pressure"
                               ,
                               "qL"
                               #, "R"
                               ,
                               "Rel_Chl_intensity"
                               # , "RFd"
                               #, "SPAD_420"
                               ,
                               "SPAD_420_intensity"
                               ,
                               "TimeofDay"
                               ,
                               "lat"
                               ,
                               "long"
                               ,
                               "waterRetention"
                               ,
                               "toughness"
                               ,
                               "elev_m"
                               ,
                               "slope_perc"
                               ,
                               "treeRich"
                               ,
                               "shrubRich"
                               ,
                               "deadDown"
                               ,
                               "precip_april_in.x"
                               ,
                               "densitometer"
                               ,
                               "shannons_flora"
                               # , "shannonsISD"
                               ,
                               "julianDate"
                               ,
                               "mean_temp_april.y"
                               ,
                               "plant_vol"
                               #,"div_raw"
                               ,
                               "sla"
                               ,
                               "compartment_EN"
                               ,
                               "habit_forb"    #
                               ,
                               "habit_graminoid"  #
                               ,
                               "habit_shrub"  #
                               ,
                               "habit_tree" #
                               ,
                               "taxon_final_Abiesconcolor"
                               ,
                               "taxon_final_Abiesgrandis"
                               ,
                               "taxon_final_Antennariamedia"
                               ,
                               "taxon_final_Aquilegiacaerula"
                               ,
                               "taxon_final_Arnicacordifolia"
                               ,
                               "taxon_final_Artemisiatridentata"
                               ,
                               "taxon_final_Astragalusalpinus"
                               ,
                               "taxon_final_Astragaluskentrophyta"
                               ,
                               "taxon_final_Astragalusmiser"
                               ,
                               "taxon_final_Bistortabistortoides"
                               ,
                               "taxon_final_Carexpaysonis"
                               ,
                               "taxon_final_Carexscopulorumvar.bracteosa"
                               ,
                               "taxon_final_Chamaenerionangustifoliumvar.angustifolium"
                               ,
                               "taxon_final_Delphiniumoccidentale"
                               ,
                               "taxon_final_Erigeronglacialis"
                               ,
                               "taxon_final_Eriogonumumbellatum"
                               ,
                               "taxon_final_Eucephaluselegans"
                               ,
                               "taxon_final_Eucephalusengelmannii"
                               ,
                               "taxon_final_Fragariavirginiana"
                               ,
                               "taxon_final_Fraseraspeciosa"
                               ,
                               "taxon_final_Geraniumviscossimumvar.viscosissimum"
                               ,
                               "taxon_final_Helianthellauniflora"
                               ,
                               "taxon_final_Heracleummaximum"
                               ,
                               "taxon_final_Juncusbalticus."
                               ,
                               "taxon_final_Juncusparryi"
                               ,
                               "taxon_final_Juncussp."
                               ,
                               "taxon_final_Juniperuscommunis"
                               ,
                               "taxon_final_Ligusticumfilicinum"
                               ,
                               "taxon_final_Lupinusargenteus"
                               ,
                               "taxon_final_Mahoniarepens"
                               ,
                               "taxon_final_Mertensianaciliatavar.ciliata"
                               ,
                               "taxon_final_Minuartiaobtusiloba"
                               ,
                               "taxon_final_Osmorhizadepauperata"
                               ,
                               "taxon_final_Oxyriadigyna"
                               ,
                               "taxon_final_Paxistimamyrsinites"
                               ,
                               "taxon_final_Piceaengelmannii"
                               ,
                               "taxon_final_Pinusalbicaulis"
                               ,
                               "taxon_final_Pinuscontorta"
                               ,
                               "taxon_final_Poapratensis"
                               ,
                               "taxon_final_Poawheeleri"
                               ,
                               "taxon_final_Polemoniumviscosum"
                               ,
                               "taxon_final_Populusangustifolium"
                               ,
                               "taxon_final_Populustremulloides"
                               ,
                               "taxon_final_Potentilladiversifolia"
                               ,
                               "taxon_final_Potentillafruticosa"
                               ,
                               "taxon_final_Potentillapulcherrima"
                               ,
                               "taxon_final_Primulaparryi"
                               ,
                               "taxon_final_Pseudotsugamenziesii"
                               ,
                               "taxon_final_Ribesmontigenum"
                               ,
                               "taxon_final_Salixglaucavar.villosa"
                               ,
                               "taxon_final_Salixreticulatavar.nana"
                               ,
                               "taxon_final_Sedumlanceolatum"
                               ,
                               "taxon_final_Symphoricarposalbus"
                               ,
                               "taxon_final_Symphyotrichumfoliaceumvar.apricum"
                               ,
                               "taxon_final_Thalictrumoccidentale"
                               ,
                               "taxon_final_Trifoliumparryivar.montanense"
                               ,
                               "taxon_final_Unknownfir"
                               ,
                               "taxon_final_Unknownpine"
                               ,
                               "taxon_final_UnknownSpruce"
                               ,
                               "taxon_final_Vacciniummembranaceum"
                               ,
                               "taxon_final_Vacciniumscoparium"
                               ,
                               "taxon_final_Wyethiaamplexicaulis"
                               ,
                               "phenology_flowering"
                               ,
                               "phenology_fruiting"
                               ,
                               "phenology_vegetative"
                               ,
                               "MEM1"
                               ,
                               "MEM2"
                             )]
  
  #Convert to numeric (makes it easier when doing imputing)
  for (i in 2:length(merged_dat)) {
    merged_dat[, i] <- as.numeric(merged_dat[, i])
  }
  
  #######################
  # Making variable for stratification #
  ######################
  
  #Make a one hot variable that is a 1 if the focal taxon (the response variable)
  #is present above its mean abund. and a zero otherwise.
  #This lets us stratify during splitting so we don't end up with a
  #train/test split that has just high or low values of the response.
  
  if (any(is.na(response_taxon))) {
    merged_dat <- merged_dat[!is.na(response_taxon), ]
    response_taxon <- response_taxon[!is.na(response_taxon)]
  }
  
  if (any(is.infinite(response_taxon))) {
    merged_dat <- merged_dat[!is.infinite(response_taxon), ]
    response_taxon <- response_taxon[!is.infinite(response_taxon)]
  }
  
  merged_dat$focal_one_hot <-
    ifelse(response_taxon > mean(response_taxon), 1, 0)
  table(merged_dat$focal_one_hot)
  
  merged_dat$stratify <-
    paste(merged_dat$compartment, merged_dat$focal_one_hot)
  
  #########################################
  #Use mlr to define a pipeline, task and a learner#
  ########################################
  
  names(merged_dat) <- gsub("\\s", "_", names(merged_dat))
  
  #Tasks are a way to organize data and the learner is the model
  #see: https://mlr3book.mlr-org.com/tasks.html
  phyllo_task = TaskRegr$new(backend = merged_dat
                             ,
                             target = as.character(focal_taxon),
                             id = "phyllo_data")
  
  #Define roles for all features
  phyllo_task$set_col_roles("sample", roles = "name")
  phyllo_task$set_col_roles("stratify", roles = "stratum")
  
  phyllo_task$col_roles$feature <-
    names(merged_dat)[!(names(merged_dat) %in% c(
      "sample",
      "focal_one_hot",
      "stratify",
      as.character(focal_taxon)
    ))]
  
  
  #Build pipeline for scaling and imputation
  imp_missind <- po("missind")
  
  #Note that imputation of numeric will NOT work with integer class features.
  #Vice versa doesn't work either
  imp_num <-
    po("imputehist", param_vals = list(affect_columns = selector_type("numeric")))
  
  
  ######
  # QC #
  ######
  # Simulate some data to see if the model works properly with a better feature.
  #
  # merged_dat$simmed <- merged_dat$response_taxon * rnorm(mean = 0, sd = 100, n = length(merged_dat$response_taxon))
  # merged_dat$simmed <- merged_dat$response_taxon * rnorm(mean = 0, sd = 10, n = length(merged_dat$response_taxon))
  
  ##########
  # RF #
  ##########
  
  graph <-
    po("imputehist", param_vals = list(affect_columns = selector_type("numeric"))) %>>%
    po("scale", param_vals = list(scale = T, center = T)) %>>%
    po(lrn("regr.ranger", importance = "permutation"))
  
  g1 <- GraphLearner$new(graph)
  
  ##########
  # tuning #
  ##########
  
  #For ranger (randomforest)
  params <- ps(
    regr.ranger.num.trees = p_int(
      lower = 50,
      upper = 400,
      tags = "budget"
    ),
    #number of trees in the dark forest
    regr.ranger.sample.fraction = p_dbl(lower = 0.01, upper = 0.3),
    #The sample.fraction parameter specifies the fraction of observations to be used in each tree. Smaller fractions lead to greater diversity, and thus less correlated trees which often is desirable.
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
    search_space = params,
    store_models = TRUE
  )
  
  #pass back to resampling for nested resampling
  outer_resampling <- rsmp("cv", folds = 2)
  
  rr <- resample(
    task = phyllo_task,
    learner = at,
    outer_resampling,
    store_models = TRUE
  )
  
  outer_resampling$instantiate(phyllo_task)
  
  extract_inner_tuning_results(rr)
  
  # #outer tuning results
  #This is the unbiased estimate of model performance
  rr$score()
  rr$aggregate(measure = msr("regr.rsq")) #can swap out different measures
  
  # Extract important features
  tained_at <- at$train(phyllo_task)
  var.imp <-
    data.frame(tained_at$model$learner$model$regr.ranger$model$variable.importance)
  
  
  #Had to build a different learner to get the backend to pass to iml later.
  #WHAT A PAIN
  graph2 <-
    po("imputehist", param_vals = list(affect_columns = selector_type("numeric"))) %>>%
    po("scale", param_vals = list(scale = T, center = T))
  
  graph2$train(phyllo_task)
  test <- graph2$predict(phyllo_task)
  forgraphing <- test$scale.output$backend
  
  engineered_df <- data.frame(forgraphing$data(
    rows = 1:dim(merged_dat)[1],
    cols = colnames(merged_dat)
  ))
  #
  # #sanity check to make sure what came through the pipeop is as expected
  # table(engineered_df$shannons_flora ==
  #         scale(merged_dat$shannons_flora, center = T, scale = T))
  
  
  library(pdp)  # for partial dependence plots
  library(ggplot2)  # for autoplot() generic
  
  p1 <- tained_at$model$learner$model$regr.ranger$model %>%
    partial(pred.var = features2examine[z],
            train = engineered_df[, !(
              names(engineered_df) %in% c("sample", "merged_dat$rich", "rich", "richISDstratum")
            )]) %>%
    autoplot(
      smooth = TRUE,
      train = engineered_df[, !(
        names(engineered_df) %in% c("sample", "merged_dat$rich", "rich", "richISDstratum")
      )],
      rug = TRUE,
      xlab = "",
      ylab = "Partial dependence"
    ) +
    theme_light()
  
  #Clunky way to get the line out of ggplot so we can use a more sensible plotting tool
  components_ggplotsux <- ggplot_build(p1)
  
  focal_taxon_line[[k]] <- components_ggplotsux$data[[3]]
  k <- k + 1
  }
  save(focal_taxon_line, file = paste("modelingResults/linesFromGGPLOT_PDPs_of_best_models", features2examine[z],".Rdata", sep = ""))
  
}
# Qucik paste of the focal variables
# features2examine <- 
#   c("deadDown",
#     "MEM2",
#     "precip_april_in.x", 
#     "lat",
#     "mean_temp_april.y", 
#     "julianDate",
#     #"elev_m",
#     "shrubRich")
# 

rm(list=ls())
load("modelingResults/linesFromGGPLOT_PDPs_of_best_modelsdeadDown.Rdata")
length(focal_taxon_line)
plot(NULL, ylim = c(0, 0.06), xlim = c(-1,3), las = 2)
for(i in 1:length(focal_taxon_line)){
  lines(focal_taxon_line[[i]])
}

#THIS LAST RDATA FILE HAS THE PRECEEDING LINES IN IT, in intervals of 18
#E.g., if we look at the first 18 lines they match the plot above. 
rm(list=ls())
load("modelingResults/linesFromGGPLOT_PDPs_of_best_modelsshrubRich.Rdata")
length(focal_taxon_line)
plot(NULL, ylim = c(0, 0.06), xlim = c(-2,3), las = 2)
for(i in 1:18){
  lines(focal_taxon_line[[i]])
}

#PLot em all
rm(list=ls())
load("modelingResults/linesFromGGPLOT_PDPs_of_best_modelsshrubRich.Rdata")
pdf(width = 8, height = 5, file = "visuals/pdp_bestModeled_taxa.pdf")
par(mfrow = c(2,4))

plotnames <-c("# of dead/down trees","MEM #2","Precip. in April",
      "Latitude",
      "Mean temp. in Apriil",
      "Sampling date",
      "Shrub richness")
for(i in 1:7){
  plot(NULL, ylim = c(0, 0.06), xlim = c(-2,3), las = 2, ylab = "Partial dependence", xlab = plotnames[i])
  
  for(j in (18*i-17):(18*i)){
    lines(focal_taxon_line[[j]])
  }
}
#tack on elevation
load("modelingResults/linesFromGGPLOT_PDPs_of_best_modelsElevation.Rdata")
length(focal_taxon_line)
plot(NULL, ylim = c(0, 0.06), xlim = c(-2,3), las = 2, ylab = "Partial dependence", xlab = "Elevation")
for(i in 1:18){
  lines(focal_taxon_line[[i]])
}
dev.off()


