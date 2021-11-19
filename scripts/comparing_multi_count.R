rm(list=ls())

counts <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                   stringsAsFactors = F)
multi <- read.csv("./processedData/16sp_estimates_wrangled_for_post_modeling_analysis.csv",
                  stringsAsFactors = F)
#QC IMPORTANT to ensure order matches
table(names(counts) == names(multi))
counts$sample <- gsub("X","",counts$sample)
table(counts$sample == multi$sample)

#Convert counts to proportions. 
counts[,3:length(counts)]  <- lapply( counts[,3:length(counts)] , 
                                     function(x) x/sum(x, na.rm=TRUE) )


plot(as.numeric(counts[1,3:length(multi)]),
  as.numeric(multi[1,3:length(counts)]),
  ylim = c(0,1),
  xlim =c(0,0.1))

abline(lm(as.numeric(multi[1,3:length(counts)])~
         as.numeric(counts[1,3:length(multi)])))

# for(i in 1:length(counts[,1])){
#   points(as.numeric(counts[i,3:length(multi)]),
#      as.numeric(multi[i,3:length(counts)]))
# 
#   abline(lm(as.numeric(multi[i,3:length(counts)])~
#             as.numeric(counts[i,3:length(multi)])))
# }

cor.test(as.numeric(multi[4,3:length(counts)]),
         as.numeric(counts[4,3:length(multi)]))
