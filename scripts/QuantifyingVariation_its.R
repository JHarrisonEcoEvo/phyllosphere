rm(list=ls())

its <- read.csv("./processedData/ITSp_estimates_wrangled_for_post_modeling_analysis_divided_by_ISD.csv")
itsCounts <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG")
itsCounts[,grep("Zotu", names(itsCounts))] <- itsCounts[,grep("Zotu", names(itsCounts))] -1
itsCounts[,grep("Zotu", names(itsCounts))] <- itsCounts[,grep("Zotu", names(itsCounts))] / itsCounts$ISD

#Calculate variation in each of these data sets
vars <- vector()
for(i in grep("Zotu", names(its))){
  vars <- c(vars, var(its[,i]))
}
summary(vars)

varsC <- vector()
for(i in grep("Zotu", names(itsCounts))){
  varsC <- c(varsC, var(itsCounts[,i]))
}
summary(varsC)

#See what the summed read counts look like across all samples. 
#Again, the minimum and max are higher for raw count data compared to multinomial
summary(colSums(its[, grep("Zotu", names(its))]))
summary(colSums(itsCounts[, grep("Zotu", names(itsCounts))]))

#Do the same thing, but just for the very rare things.
#Curious if there is more variation in the rarer stuff for count data.
lowCounts <- which(colSums(its[, grep("Zotu", names(its))]) < 50)

summary(colSums(its[, grep("Zotu", names(its))])[lowCounts])
summary(colSums(itsCounts[, grep("Zotu", names(itsCounts))])[lowCounts])

#calculate variation
var(colSums(its[, grep("Zotu", names(its))])[lowCounts])
var(colSums(itsCounts[, grep("Zotu", names(itsCounts))])[lowCounts])


#TRY FOR ABUNDANT STUFF

highCounts <- which(colSums(its[, grep("Zotu", names(its))]) > 200)

summary(colSums(its[, grep("Zotu", names(its))])[highCounts])
summary(colSums(itsCounts[, grep("Zotu", names(itsCounts))])[highCounts])

#calculate variation
var(colSums(its[, grep("Zotu", names(its))])[highCounts])
var(colSums(itsCounts[, grep("Zotu", names(itsCounts))])[highCounts])

#This analysis is crude, but it seems that the lows are not as low and the highs are not as high for the 
#multinomial data compared to the Hellinger data.