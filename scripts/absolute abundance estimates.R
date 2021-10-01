rm(list=ls())
its <- read.table("processedData/smallmem97_ITS.sintax", fill = T)
head(its)

its_otus <- read.csv("./processedData/otuTables/smallmem97_ITS_for_modeling_rearranged_for_CNVRG",
                     stringsAsFactors = F)
length(names(its_otus)) -4

its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] - 1
its_otus[,3:length(its_otus)] <- its_otus[,3:length(its_otus)] / its_otus$ISD

table(rowSums(its_otus[grep("EN",its_otus$sample),3:(length(its_otus)-2)]) <= 1)
268/(935+268) #22%
table(rowSums(its_otus[grep("EP",its_otus$sample),3:(length(its_otus)-2)]) <= 1)
73/(73+1141) #6%

#How many molecules does this translate to
#0.03pg x 6 is our ISD = 0.18pg
#our ITS ISD weighs this much
#51656.44 Da #https://www.bioinformatics.org/sms2/dna_mw.html
#8.577752914e-8 pg is how much this many Da weighs  
#so if we divide 0.18 by this we will get how many molecules we have
#.18/8.577752914e-8
#2,098,452

#The 16s ISD weighs 50293.50 Da 
#which is 8.35143142e-8 pg
#.18 / 8.35143142e-8
#2155319

#Lets check out bacteria
bacts <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                  stringsAsFactors = F)

tail(names(bacts))
names(bacts)[length(names(bacts))-4]


bacts[,3:length(bacts)] <- bacts[,3:length(bacts)] - 1
bacts[,3:length(bacts)] <- bacts[,3:length(bacts)] / bacts$ISD
table(rowSums(bacts[grep("EN",bacts$sample),3:(length(bacts)-4)]) <= 1)
150/(150+1003) #13%
table(rowSums(bacts[grep("EP",bacts$sample),3:(length(bacts)-4)]) <= 1)
124/(124+1073) #10%

table(bacts$plant <= 1)
1107/(sum(table(bacts$plant <= 1)))

summary(rowSums(bacts[grep("EN",bacts$sample),3:(length(bacts)-4)])) 
summary(rowSums(its_otus[grep("EN",its_otus$sample),3:(length(its_otus)-2)])) 

summary(rowSums(bacts[grep("EP",bacts$sample),3:(length(bacts)-4)])) 
summary(rowSums(its_otus[grep("EP",its_otus$sample),3:(length(its_otus)-2)])) 

bacts <- read.csv("processedData/otuTables/smallmem97_16S_for_modeling_rearranged_for_CNVRG",
                  stringsAsFactors = F)
sum(bacts$plant)
