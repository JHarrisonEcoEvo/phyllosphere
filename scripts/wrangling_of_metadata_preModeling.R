#This code is terrifically confusing. Wrangling for this project was a nightmare.
#Horrible

rm(list=ls())
#load the two demux keys, edit and combine
demux <- read.csv("raw_data_notIncluding_sequenceData/demultiplexingKeys/novaseq3_demux.csv", 
                  header = T, stringsAsFactors = F)
demux$forward_barcode <-toupper(demux$forward_barcode)
demux$reverse_barcode <-toupper(demux$reverse_barcode)
demux$combo <- paste(demux$locus, 
                     demux$forward_barcode, 
                     demux$reverse_barcode, 
                     demux$samplename, sep = "_")
demux2 <- read.csv("raw_data_notIncluding_sequenceData/demultiplexingKeys/NovaSeq2_DemuxJH.csv", 
                   header = T, stringsAsFactors = F)
demux2$forward_barcode <-toupper(demux2$forward_barcode)
demux2$reverse_barcode <-toupper(demux2$reverse_barcode)
demux2$combo <- paste(demux2$locus, 
                     demux2$forward_barcode, 
                     demux2$reverse_barcode, 
                     demux2$samplename, sep = "_")

# demux2[intersect(grep("GATACGCCAA", demux2$forward_barcode), 
# grep("GACGAGAC", demux2$reverse_barcode)),]
# 
# demux2[intersect(grep("TCCATATCAA", demux2$forward_barcode), 
#                  grep("GACCGCGA", demux2$reverse_barcode)),]
# 
# demux2[intersect(grep("GATGCGCA", demux2$forward_barcode), 
#                  grep("GACCGCGA", demux2$reverse_barcode)),]
# 
# demux2[intersect(grep("TTAACTCCAA", demux2$forward_barcode), 
#                  grep("GACGAGAC", demux2$reverse_barcode)),]

#remove extraneous projects
demux2 <- demux2[demux2$project == "Harrison",]
demux <- demux[demux$project == "sm18biog",]
demux2 <- demux2[,-(which(names(demux2) == "key"))]
demux_all <- rbind(demux2, demux)
table(duplicated(demux_all$combo)) #all false, as expected.
table(demux_all$project)
table(demux_all$locus)

demux_all$plant <- gsub("[A-Za-z_]*(\\d+_\\d+_\\d+_\\d+)\\w+", "\\1", demux_all$samplename)

#table(demux_all$plant)[which(table(demux_all$plant) > 8)]

#merge with leaf area and size information
dim(demux_all)
leafarea <- read.csv("raw_data_notIncluding_sequenceData/area_mg_leafCount.csv", 
                     header=T, stringsAsFactors =F)
dim(leafarea)
leafarea <- leafarea[leafarea$notes != "dupe",]

demux_all_leaf <- merge(demux_all, leafarea, by.x = "plant", by.y = "label", all.x = T)
dim(demux_all_leaf)

locus <- read.csv("raw_data_notIncluding_sequenceData/LocationSize.csv", header=T, stringsAsFactors =F)
dim(locus)
demux_all_leaf_locus <- merge(demux_all_leaf, locus, by.x = "plant", by.y = "label",  all.x = T)
dim(demux_all_leaf_locus)

#bring in the multispeq
multi <- read.csv("raw_data_notIncluding_sequenceData/multispeq_data.csv", header=T, stringsAsFactors =F)
multi <- multi[-grep("[a-z]", multi$label),]
dim(multi)
multi <- multi[-grep("dupe",multi$note),]

demux_all_leaf_locus_multi <- merge(demux_all_leaf_locus, multi, 
                                    by.x = "plant", by.y = "label", 
                                    all.x = T)
#multi[!(multi$label %in% demux_all_leaf_locus$plant),]
#demux_all_leaf_locus$plant[!(demux_all_leaf_locus$plant%in% multi$label)]

dim(demux_all_leaf_locus)

#bring in leaf toughness and water retention
tough <- read.csv("raw_data_notIncluding_sequenceData/Toughness_waterRetention.csv",
                  header=T, stringsAsFactors =F)
dim(tough)

demux_all_leaf_locus_multi_tough <- merge(demux_all_leaf_locus_multi, tough,
                                          by.x = "plant", by.y = "sample", 
                                    all.x = T)

dim(demux_all_leaf_locus_multi_tough)


#bring in site data
site <- read.csv("raw_data_notIncluding_sequenceData/siteData.csv", 
                 header=T, stringsAsFactors =F)
#make region and site labels for data
demux_all_leaf_locus_multi_tough$region_site <- gsub("(\\d+_\\d+)_\\d+_\\d+","\\1",
                                                     demux_all_leaf_locus_multi_tough$plant)
#table(demux_all_leaf_locus_multi_tough$region_site )
demux_all_leaf_locus_multi_tough_site <- merge(demux_all_leaf_locus_multi_tough, site, 
                                               by.x = "region_site", by.y = "label", 
                                          all.x = T)
dim(demux_all_leaf_locus_multi_tough_site)

#bring in plant taxon data
taxa <- read.csv("raw_data_notIncluding_sequenceData/TaxaSampled.csv", header=T, 
                 stringsAsFactors =F)
#make region and site labels for data
demux_all_leaf_locus_multi_tough_site$region_site_plant <-
  gsub("(\\d+_\\d+_\\d+)_\\d+","\\1",demux_all_leaf_locus_multi_tough_site$plant)

demux_all_leaf_locus_multi_tough_site_taxa <- merge(demux_all_leaf_locus_multi_tough_site,
                                                    taxa, by.x = "region_site_plant",
                                                    by.y = "label", 
                                               all.x = T)
dim(demux_all_leaf_locus_multi_tough_site_taxa)

demux_all_leaf_locus_multi_tough_site_taxa <- 
  demux_all_leaf_locus_multi_tough_site_taxa[,
    !names(demux_all_leaf_locus_multi_tough_site_taxa) %in%
  c("region", "site", "species", "indiv.x", "indiv.y" )]


#POSTHOC ADDITION:
#Used this code to figure out that I want to remove the 2_1_1_2ep sample from consideration
# dat <- read.csv("./processedData/otuTables/its90_otuTableCLEAN", 
#                 fill = T, header = T, stringsAsFactors = F)
# metadat <- read.csv("./processedData/metadata.csv", stringsAsFactors = F)
# metadat <- metadat[metadat$locus == "ITS",]
# metadat$combo[grep("2_1_1_2", metadat$combo)]
# colSums(dat[,grep("2_1_1_2",names(dat))])
# dat[which(dat[,grep("ITS_CGAACGCA_ACTCCGCCA_2_1_1_2ep", names(dat))] > 0),
#     c(1,grep("ITS_CGAACGCA_ACTCCGCCA_2_1_1_2ep", names(dat)))]
# #unclear if this is a weak EN or a strong EP, so removing it. Same goes for other of these samples
#best to just not include them in analysis.

#Addding this here in case the metadata  needs to be remade. As of Feb 8, 2021 I removed this line in the
#excel file. But if metadat needs edited, then I wanted this removal to be automated.
demux_all_leaf_locus_multi_tough_site_taxa <- 
  demux_all_leaf_locus_multi_tough_site_taxa[-grep("2_1_1_2",
                                                   demux_all_leaf_locus_multi_tough_site_taxa$combo),]


#assign compartment to all samples, either EN or EP
demux_all_leaf_locus_multi_tough_site_taxa$compartment <- NA 
demux_all_leaf_locus_multi_tough_site_taxa$compartment[grep("ep|EP",demux_all_leaf_locus_multi_tough_site_taxa$samplename)] <- "EP"
demux_all_leaf_locus_multi_tough_site_taxa$compartment[grep("en|EN",demux_all_leaf_locus_multi_tough_site_taxa$samplename)] <- "EN"

dim(demux_all_leaf_locus_multi_tough_site_taxa)

#demux_all_leaf_locus_multi_tough_site_taxa$samplename[grep("[MO]$",
#                                                         demux_all_leaf_locus_multi_tough_site_taxa$samplename)]

#fix misassignment of soil samples and blanks to plants
demux_all_leaf_locus_multi_tough_site_taxa$substrate[grep("[MO]$",
      demux_all_leaf_locus_multi_tough_site_taxa$samplename)] <- "soil"

demux_all_leaf_locus_multi_tough_site_taxa$substrate[grep("Harrison_soil",
                                                          demux_all_leaf_locus_multi_tough_site_taxa$samplename)] <- "soil"

demux_all_leaf_locus_multi_tough_site_taxa$substrate[grep("lank",
                                                          demux_all_leaf_locus_multi_tough_site_taxa$samplename)] <- "blank"

#check that compartment assignment worked, it didnt. Will check again later. 

demux_all_leaf_locus_multi_tough_site_taxa[is.na(demux_all_leaf_locus_multi_tough_site_taxa$compartment) & 
                                             demux_all_leaf_locus_multi_tough_site_taxa$substrate == "plant",1:10]


#Remove hyphens. They are only present in the names of a few soil samples
demux_all_leaf_locus_multi_tough_site_taxa$samplename <- 
  gsub("-","_",  demux_all_leaf_locus_multi_tough_site_taxa$samplename)

#make lower case en and ep uppercase and get the underscore in there
demux_all_leaf_locus_multi_tough_site_taxa$samplename <- 
  gsub("(\\d+)en","\\1_EN",  demux_all_leaf_locus_multi_tough_site_taxa$samplename)
demux_all_leaf_locus_multi_tough_site_taxa$samplename <- 
  gsub("(\\d+)ep","\\1_EP",  demux_all_leaf_locus_multi_tough_site_taxa$samplename)

#REname the dataset for ease of typing
metadata <- demux_all_leaf_locus_multi_tough_site_taxa

#Get rid of stuff that wasn't sequenced
dim(metadata)
metadata <- metadata[!is.na(metadata$forward_barcode),]

#Extract duplicate from sample name and put that info in a new field
metadata$dupe <- "no"
metadata$dupe[grep("dupe", metadata$samplename)] <- "yes"

#clean up rewashed field
metadata$rewashed[is.na(metadata$rewashed)] <- "no"
metadata$rewashed[metadata$rewashed == ""] <- "no"

#Get rid of samples from other project that have 19 as site
metadata <- metadata[metadata$region_site != 19,]
dim(metadata)

#posthoc QC: all have locus, but still some don't have compartment. Will check again later.
demux_all_leaf_locus_multi_tough_site_taxa$locus[is.na(demux_all_leaf_locus_multi_tough_site_taxa$compartment) & 
                                             demux_all_leaf_locus_multi_tough_site_taxa$substrate == "plant"]


#for samples that don't have locus in the name add that back in there
metadata$samplename[grep("_[ENP]*$", metadata$samplename)] <-
      paste(metadata$samplename[grep("_[ENP]*$", metadata$samplename)], "_",
  metadata$locus[grep("_[ENP]*$", metadata$samplename)], sep = "")

#clean up specific sample names that were incorrect

if(dim(metadata)[1] != 10678){
  print("THIS IS VERY FRAGILE...! THE DATASET CHANGED IN DIMENSIONS, so this is probably broken")
}

metadata$samplename[grep("1_3_2_2", metadata$samplename)]
if(metadata$samplename[1184] == "1_3_2_2"){
  metadata$samplename[1184]  <- "1_3_2_2_EN_16S"
}
if(metadata$samplename[1186] == "1_3_2_2"){
  metadata$samplename[1186]  <- "1_3_2_2_EN_ITS"
}
if(metadata$samplename[1189] == "1_3_2_2"){
  metadata$samplename[1189]  <- "1_3_2_2_EN_ITS"
}
if(metadata$samplename[1195] == "1_3_2_2"){
  metadata$samplename[1195]  <- "1_3_2_2_EN_ITS"
}

#switching to more robust recoding method. Will fix above if I need to redo this heinous script
metadata[grep("1_1_3_9", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "ITS_TGAGAAGAA_CTCAGGTAA_1_1_3_9"] <- "1_1_3_9_EN_ITS"
metadata$samplename[metadata$combo == "16S_CATGCAGAA_CTCAGGTAA_1_1_3_9"] <- "1_1_3_9_EN_16S"
metadata$samplename[metadata$combo == "16S_TGAGAAGAA_CTCAGGTAA_1_1_3_9"] <- "1_1_3_9_EN_16S"
metadata$samplename[metadata$combo == "ITS_CATGCAGAA_CTCAGGTAA_1_1_3_9"] <- "1_1_3_9_EN_ITS"

metadata[grep("1_2_1_4", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "ITS_ACGAGTCAA_CATTACGCCA_1_2_1_4"] <- "1_2_1_4_EN_ITS"
metadata$samplename[metadata$combo == "16S_GGTAATGCAA_CATTACGCCA_1_2_1_4"] <- "1_2_1_4_EN_16S"
metadata$samplename[metadata$combo == "16S_ACGAGTCAA_CATTACGCCA_1_2_1_4"] <- "1_2_1_4_EN_16S"
metadata$samplename[metadata$combo == "ITS_GGTAATGCAA_CATTACGCCA_1_2_1_4"] <- "1_2_1_4_EN_ITS"

metadata[grep("1_2_7_6", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_GCATCAGAA_AGCAATGAA_1_2_7_6"] <- "1_2_7_6_EN_16S"
metadata$samplename[metadata$combo == "16S_AGGCCAGAA_AGCAATGAA_1_2_7_6"] <- "1_2_7_6_EN_16S"
metadata$samplename[metadata$combo == "ITS_GCATCAGAA_AGCAATGAA_1_2_7_6"] <- "1_2_7_6_EN_ITS"
metadata$samplename[metadata$combo == "ITS_AGGCCAGAA_AGCAATGAA_1_2_7_6"] <- "1_2_7_6_EN_ITS"

metadata[grep("1_3_1_7", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_CGCTGAGAA_GGAATAGA_1_3_1_7"] <- "1_3_1_7_EN_16S"
metadata$samplename[metadata$combo == "16S_ACGCCGCA_GGAATAGA_1_3_1_7"] <- "1_3_1_7_EN_16S"
metadata$samplename[metadata$combo == "ITS_CGCTGAGAA_GGAATAGA_1_3_1_7"] <- "1_3_1_7_EN_ITS"
metadata$samplename[metadata$combo == "ITS_ACGCCGCA_GGAATAGA_1_3_1_7"] <- "1_3_1_7_EN_ITS"

metadata[grep("1_3_5_7", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_GTCGACCA_CAGTCCTAA_1_3_5_7"] <- "1_3_5_7_EN_16S"
metadata$samplename[metadata$combo == "16S_TCTCGCATAA_CAGTCCTAA_1_3_5_7"] <- "1_3_5_7_EN_16S"
metadata$samplename[metadata$combo == "ITS_TCTCGCATAA_CAGTCCTAA_1_3_5_7"] <- "1_3_5_7_EN_ITS"
metadata$samplename[metadata$combo == "ITS_GTCGACCA_CAGTCCTAA_1_3_5_7"] <- "1_3_5_7_EN_ITS"

metadata[grep("2_2_3_10", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_CGAACGCA_CCTAGAGA_2_2_3_10"] <- "2_2_3_10_EN_16S"
metadata$samplename[metadata$combo == "16S_ATAAGAGAA_CCTAGAGA_2_2_3_10"] <- "2_2_3_10_EN_16S"
metadata$samplename[metadata$combo == "ITS_ATAAGAGAA_CCTAGAGA_2_2_3_10"] <- "2_2_3_10_EN_ITS"
metadata$samplename[metadata$combo == "ITS_CGAACGCA_CCTAGAGA_2_2_3_10"] <- "2_2_3_10_EN_ITS"

metadata[grep("2_2_5_1", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_ATAATCCA_TGCTGCTA_2_2_5_1"] <- "2_2_5_1_EN_16S"
metadata$samplename[metadata$combo == "16S_CGAACGCA_TGCTGCTA_2_2_5_1"] <- "2_2_5_1_EN_16S"
metadata$samplename[metadata$combo == "ITS_ATAATCCA_TGCTGCTA_2_2_5_1"] <- "2_2_5_1_EN_ITS"
metadata$samplename[metadata$combo == "ITS_CGAACGCA_TGCTGCTA_2_2_5_1"] <- "2_2_5_1_EN_ITS"

metadata[grep("3_1_2_9", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_TCCATATCAA_TGAGCTGCCA_3_1_2_9"] <- "3_1_2_9_EN_16S"
metadata$samplename[metadata$combo == "16S_GGTAATGCAA_TGAGCTGCCA_3_1_2_9"] <- "3_1_2_9_EN_16S"
metadata$samplename[metadata$combo == "ITS_GGTAATGCAA_TGAGCTGCCA_3_1_2_9"] <- "3_1_2_9_EN_ITS"
metadata$samplename[metadata$combo == "ITS_TCCATATCAA_TGAGCTGCCA_3_1_2_9"] <- "3_1_2_9_EN_ITS"

metadata[grep("3_3_5_8", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_GGCGCCGAA_GCGCTATAA_3_3_5_8"] <- "3_3_5_8_EN_16S"
metadata$samplename[metadata$combo == "16S_TTATTCGAA_GCGCTATAA_3_3_5_8"] <- "3_3_5_8_EN_16S"
metadata$samplename[metadata$combo == "ITS_GGCGCCGAA_GCGCTATAA_3_3_5_8"] <- "3_3_5_8_EN_ITS"
metadata$samplename[metadata$combo == "ITS_TTATTCGAA_GCGCTATAA_3_3_5_8"] <- "3_3_5_8_EN_ITS"

metadata[grep("3_2_2_5", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_GGAATGGAA_TGGACTGA_3_2_2_5"] <- "3_2_2_5_EN_16S"
metadata$samplename[metadata$combo == "16S_GAGCTTCA_TGGACTGA_3_2_2_5"] <- "3_2_2_5_EN_16S"
metadata$samplename[metadata$combo == "ITS_GAGCTTCA_TGGACTGA_3_2_2_5"] <- "3_2_2_5_EN_ITS"
metadata$samplename[metadata$combo == "ITS_GGAATGGAA_TGGACTGA_3_2_2_5"] <- "3_2_2_5_EN_ITS"

metadata[grep("3_2_7_2", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_CGCCTCCA_AATCCAGA_3_2_7_2"] <- "3_2_7_2_EN_16S"
metadata$samplename[metadata$combo == "16S_ACGCCGCA_AATCCAGA_3_2_7_2"] <- "3_2_7_2_EN_16S"
metadata$samplename[metadata$combo == "ITS_CGCCTCCA_AATCCAGA_3_2_7_2"] <- "3_2_7_2_EN_ITS"
metadata$samplename[metadata$combo == "ITS_ACGCCGCA_AATCCAGA_3_2_7_2"] <- "3_2_7_2_EN_ITS"

metadata[grep("4_3_1_7", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_GGCGCCGAA_GATAGATAA_4_3_1_7"] <- "4_3_1_7_EN_16S"
metadata$samplename[metadata$combo == "16S_TAGAACGAA_GATAGATAA_4_3_1_7"] <- "4_3_1_7_EN_16S"
metadata$samplename[metadata$combo == "ITS_TAGAACGAA_GATAGATAA_4_3_1_7"] <- "4_3_1_7_EN_ITS"
metadata$samplename[metadata$combo == "ITS_GGCGCCGAA_GATAGATAA_4_3_1_7"] <- "4_3_1_7_EN_ITS"

metadata[grep("6_1_4_6", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_ATAAGAGAA_TTCGTAGA_6_1_4_6"] <- "6_1_4_6_EN_16S"
metadata$samplename[metadata$combo == "16S_TCCATAGAA_TTCGTAGA_6_1_4_6"] <- "6_1_4_6_EN_16S"
metadata$samplename[metadata$combo == "ITS_TCCATAGAA_TTCGTAGA_6_1_4_6"] <- "6_1_4_6_EN_ITS"
metadata$samplename[metadata$combo == "ITS_ATAAGAGAA_TTCGTAGA_6_1_4_6"] <- "6_1_4_6_EN_ITS"

metadata[grep("7_1_1_4", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_CTTGCAGCAA_ACCTAGTAA_7_1_1_4"] <- "7_1_1_4_EN_16S"
metadata$samplename[metadata$combo == "16S_AATTACCA_ACCTAGTAA_7_1_1_4"] <- "7_1_1_4_EN_16S"
metadata$samplename[metadata$combo == "ITS_AATTACCA_ACCTAGTAA_7_1_1_4"] <- "7_1_1_4_EN_ITS"
metadata$samplename[metadata$combo == "ITS_CTTGCAGCAA_ACCTAGTAA_7_1_1_4"] <- "7_1_1_4_EN_ITS"

metadata[grep("6_2_5_2", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_ACGAGTCAA_ACCGTTACCA_6_2_5_2"] <- "6_2_5_2_EN_16S"
metadata$samplename[metadata$combo == "16S_TCCTCTTGAA_ACCGTTACCA_6_2_5_2"] <- "6_2_5_2_EN_16S"
metadata$samplename[metadata$combo == "ITS_ACGAGTCAA_ACCGTTACCA_6_2_5_2"] <- "6_2_5_2_EN_ITS"
metadata$samplename[metadata$combo == "ITS_TCCTCTTGAA_ACCGTTACCA_6_2_5_2"] <- "6_2_5_2_EN_ITS"

metadata[grep("6_2_6_5", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_CTTCATCA_AAGCCAACCA_6_2_6_5"] <- "6_2_6_5_EN_16S"
metadata$samplename[metadata$combo == "16S_GATGCGCA_AAGCCAACCA_6_2_6_5"] <- "6_2_6_5_EN_16S"
metadata$samplename[metadata$combo == "ITS_GATGCGCA_AAGCCAACCA_6_2_6_5"] <- "6_2_6_5_EN_ITS"
metadata$samplename[metadata$combo == "ITS_CTTCATCA_AAGCCAACCA_6_2_6_5"] <- "6_2_6_5_EN_ITS"

metadata[grep("6_2_5_8", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_GGAATGGAA_TAGTAATA_6_2_5_8"] <- "6_2_5_8_EN_16S"
metadata$samplename[metadata$combo == "16S_AGAACCGCAA_TAGTAATA_6_2_5_8"] <- "6_2_5_8_EN_16S"
metadata$samplename[metadata$combo == "ITS_GGAATGGAA_TAGTAATA_6_2_5_8"] <- "6_2_5_8_EN_ITS"
metadata$samplename[metadata$combo == "ITS_AGAACCGCAA_TAGTAATA_6_2_5_8"] <- "6_2_5_8_EN_ITS"

metadata[grep("6_3_1_10", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_CCGTCAAGAA_CAACTTGAA_6_3_1_10"] <- "6_3_1_10_EN_16S"
metadata$samplename[metadata$combo == "16S_GCATCAGAA_CAACTTGAA_6_3_1_10"] <- "6_3_1_10_EN_16S"
metadata$samplename[metadata$combo == "ITS_GCATCAGAA_CAACTTGAA_6_3_1_10"] <- "6_3_1_10_EN_ITS"
metadata$samplename[metadata$combo == "ITS_CCGTCAAGAA_CAACTTGAA_6_3_1_10"] <- "6_3_1_10_EN_ITS"

metadata[grep("6_3_2_7", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_TTGGTAAGAA_CCTCGAACCA_6_3_2_7"] <- "6_3_2_7_EN_16S"
metadata$samplename[metadata$combo == "16S_CTCCGAAGAA_CCTCGAACCA_6_3_2_7"] <- "6_3_2_7_EN_16S"
metadata$samplename[metadata$combo == "ITS_CTCCGAAGAA_CCTCGAACCA_6_3_2_7"] <- "6_3_2_7_EN_ITS"
metadata$samplename[metadata$combo == "ITS_TTGGTAAGAA_CCTCGAACCA_6_3_2_7"] <- "6_3_2_7_EN_ITS"

metadata[grep("6_3_1_4", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_ACTGACGAA_CCAGGATAA_6_3_1_4"] <- "6_3_1_4_EN_16S"
metadata$samplename[metadata$combo == "16S_CCGCGCGAA_CCAGGATAA_6_3_1_4"] <- "6_3_1_4_EN_16S"
metadata$samplename[metadata$combo == "ITS_CCGCGCGAA_CCAGGATAA_6_3_1_4"] <- "6_3_1_4_EN_ITS"
metadata$samplename[metadata$combo == "ITS_ACTGACGAA_CCAGGATAA_6_3_1_4"] <- "6_3_1_4_EN_ITS"

metadata[grep("6_3_5_10", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_ATAATCCA_ATCCGATA_6_3_5_10"] <- "6_3_5_10_EN_16S"
metadata$samplename[metadata$combo == "16S_AGAACCGCAA_ATCCGATA_6_3_5_10"] <- "6_3_5_10_EN_16S"
metadata$samplename[metadata$combo == "ITS_AGAACCGCAA_ATCCGATA_6_3_5_10"] <- "6_3_5_10_EN_ITS"
metadata$samplename[metadata$combo == "ITS_ATAATCCA_ATCCGATA_6_3_5_10"] <- "6_3_5_10_EN_ITS"

metadata[grep("7_2_2_3", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_CTCCGAAGAA_CGCTCAACCA_7_2_2_3"] <- "7_2_2_3_EN_16S"
metadata$samplename[metadata$combo == "16S_CCGTCAAGAA_CGCTCAACCA_7_2_2_3"] <- "7_2_2_3_EN_16S"
metadata$samplename[metadata$combo == "ITS_CTCCGAAGAA_CGCTCAACCA_7_2_2_3"] <- "7_2_2_3_EN_ITS"
metadata$samplename[metadata$combo == "ITS_CCGTCAAGAA_CGCTCAACCA_7_2_2_3"] <- "7_2_2_3_EN_ITS"

metadata[grep("7_2_2_7", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_TTATTCGAA_GTAACCTAA_7_2_2_7"] <- "7_2_2_7_EN_16S"
metadata$samplename[metadata$combo == "16S_CTTGCAGCAA_GTAACCTAA_7_2_2_7"] <- "7_2_2_7_EN_16S"
metadata$samplename[metadata$combo == "ITS_CTTGCAGCAA_GTAACCTAA_7_2_2_7"] <- "7_2_2_7_EN_ITS"
metadata$samplename[metadata$combo == "ITS_TTATTCGAA_GTAACCTAA_7_2_2_7"] <- "7_2_2_7_EN_ITS"

metadata[grep("7_2_1_5", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_TAGAACGAA_AAGGAATAA_7_2_1_5"] <- "7_2_1_5_EN_16S"
metadata$samplename[metadata$combo == "16S_CGCTGAGAA_AAGGAATAA_7_2_1_5"] <- "7_2_1_5_EN_16S"
metadata$samplename[metadata$combo == "ITS_CGCTGAGAA_AAGGAATAA_7_2_1_5"] <- "7_2_1_5_EN_ITS"
metadata$samplename[metadata$combo == "ITS_TAGAACGAA_AAGGAATAA_7_2_1_5"] <- "7_2_1_5_EN_ITS"

#Fix some of the stuff that had the "check post seq" note. 
#Turns out these samples are mostly EPs, but are highly correlated with ENs
#Which is interesting. 
metadata$samplename <- gsub("(1_2_6_._)EN(_[IT16]*S_[12]_primary)_CHECK_post_seq",
     "\\1EP\\2", metadata$samplename)

metadata$samplename <- gsub("(1_2_7_\\d+_)EN(_[IT16]*S_[12]_primary)_CHECK_post_seq",
                            "\\1EP\\2", metadata$samplename)

metadata$samplename <- gsub("(1_2_7_5_)EN(_[16SIT]*).*primary",
                            "\\1EP\\2", metadata$samplename)

metadata$samplename[metadata$plate=="Ha8" & 
                      metadata$samplename=="1_2_7_4_EN_ITS"] <- "1_2_7_4_EP_ITS"

metadata$samplename <- gsub("(1_2_7_._)EN(_[IT16]*S_[12]_primary)_CHECK_post_seq",
                            "\\1EP\\2", metadata$samplename)

metadata <- metadata[-grep("p2", metadata$samplename),]
metadata <- metadata[-grep("p3", metadata$samplename),]

#grep("p2", metadata$samplename, value = T)

#2_1_1_5 Decided was ok that there was extra dupes of EP bc the whole plate was EP

#After I get the matching script in "combine_pcr..." to function for 16s then do this:
#Need to check 1_2_7_4 16S and make sure Ha8 should be EP

# 
# metadata$samplename[grep("7_2_1_5$", metadata$samplename)]
# metadata[10304,]

#cehck out dupes
table(metadata$samplename[metadata$dupe == "no"])[
  which(table(metadata$samplename[metadata$dupe == "no"]) > 2)]

#Spot checked a bunch of these and looks like primarily the extras are duplicates
#that are properly labeled.
table(demux_all$plant)[which(table(demux_all$plant) > 8)]

metadata[metadata$samplename == "1_3_2_2_EN_ITS",1:10]

# 1_3_2_2_EN_ITS #seems good
# 2_1_1_4_EP_rewash
# 2_1_1_7_EP_rewash
# 2_1_2_1_EP_rewash 
# 2_1_2_2_EP_rewash 
# 2_1_3_3_EP_rewash 
# 2_1_3_5_EP_rewash 
# 2_1_4_2_EP_rewash
# 2_1_4_6 
# 2_1_6_2_EP_rewash 
# 2_1_6_5_EP_rewash 
# 2_1_6_7_EP_rewash 
# 2_1_6_9_EP_rewash
# 2_2_5_1_EN_16S
#metadata[grep( "2_2_5_1",metadata$samplename),1:10]

# 2_2_5_1_EN_ITS 
# 3_1_2_1O
# 3_2_7_4 *
# 3_3_1_8_EP_rewash 
# 3_3_1_9 
# 3_3_2_10_EP_rewash
# 3_3_2_5_EP_rewash 
# 3_3_3_10
# 4_2_5_4_EP_rewash
# 5_2_7_6 *
# 6_3_1_3

#COme back to. All have extras, because whole plates were sequenced twice
#2_1_1_1 #fine
#2_1_1_10 #fine
#2_1_1_3 #fine
#2_1_4_6 as I need to look at reads to see what is EN v EP


metadata[grep("2_1_4_6", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_TCCGCTCA_GACCGCGA_2_1_4_6"] <- "2_1_4_6_EN_16S"
metadata$samplename[metadata$combo == "16S_CTTCATCA_GACCGCGA_2_1_4_6"] <- "2_1_4_6_EN_16S"
metadata$samplename[metadata$combo == "ITS_TCCGCTCA_GACCGCGA_2_1_4_6"] <- "2_1_4_6_EN_ITS"
metadata$samplename[metadata$combo == "ITS_CTTCATCA_GACCGCGA_2_1_4_6"] <- "2_1_4_6_EN_ITS"

#Checked the otu table names all the way through to 2_1_6_9
#and all of them seemed fine.

metadata[grep("6_3_1_3", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_GAGCTTCA_TATAAGGA_6_3_1_3"] <- "6_3_1_3_EN_16S"
metadata$samplename[metadata$combo == "16S_TCCGCTCA_TATAAGGA_6_3_1_3"] <- "6_3_1_3_EN_16S"
metadata$samplename[metadata$combo == "ITS_GAGCTTCA_TATAAGGA_6_3_1_3"] <- "6_3_1_3_EN_ITS"
metadata$samplename[metadata$combo == "ITS_TCCGCTCA_TATAAGGA_6_3_1_3"] <- "6_3_1_3_EN_ITS"

metadata[grep("3_3_1_9", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_TCCATATCAA_ATAGTCGTAA_3_3_1_9"] <- "3_3_1_9_EN_16S"
metadata$samplename[metadata$combo == "16S_GATGCGCA_ATAGTCGTAA_3_3_1_9"] <- "3_3_1_9_EN_16S"
metadata$samplename[metadata$combo == "ITS_TCCATATCAA_ATAGTCGTAA_3_3_1_9"] <- "3_3_1_9_EN_ITS"
metadata$samplename[metadata$combo == "ITS_GATGCGCA_ATAGTCGTAA_3_3_1_9"] <- "3_3_1_9_EN_ITS"

metadata[grep("3_3_3_10", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_TCCTCTTGAA_CTGCAGACCA_3_3_3_10"] <- "3_3_3_10_EN_16S"
metadata$samplename[metadata$combo == "16S_TTGGTAAGAA_CTGCAGACCA_3_3_3_10"] <- "3_3_3_10_EN_16S"
metadata$samplename[metadata$combo == "ITS_TTGGTAAGAA_CTGCAGACCA_3_3_3_10"] <- "3_3_3_10_EN_ITS"
metadata$samplename[metadata$combo == "ITS_TCCTCTTGAA_CTGCAGACCA_3_3_3_10"] <- "3_3_3_10_EN_ITS"

metadata[grep("3_2_7_4", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_AATTACCA_GCTGCGTAA_3_2_7_4"] <- "3_2_7_4_EN_16S"
metadata$samplename[metadata$combo == "16S_AGGCCAGAA_GCTGCGTAA_3_2_7_4"] <- "3_2_7_4_EN_16S"
metadata$samplename[metadata$combo == "ITS_AATTACCA_GCTGCGTAA_3_2_7_4"] <- "3_2_7_4_EN_ITS"
metadata$samplename[metadata$combo == "ITS_AGGCCAGAA_GCTGCGTAA_3_2_7_4"] <- "3_2_7_4_EN_ITS"

metadata[grep("5_2_7_6", metadata$samplename),c(14,1, 4:10,171,172)]

metadata$samplename[metadata$combo == "16S_TGAGAAGAA_TGGACGTAA_5_2_7_6"] <- "5_2_7_6_EN_16S"
metadata$samplename[metadata$combo == "16S_GTCGACCA_TGGACGTAA_5_2_7_6"] <- "5_2_7_6_EN_16S"
metadata$samplename[metadata$combo == "ITS_TGAGAAGAA_TGGACGTAA_5_2_7_6"] <- "5_2_7_6_EN_ITS"
metadata$samplename[metadata$combo == "ITS_GTCGACCA_TGGACGTAA_5_2_7_6"] <- "5_2_7_6_EN_ITS"

metadata$substrate[grep("48_6_3_24_37", metadata$samplename)] <- "soil"

metadata$substrate[grep("lank", metadata$samplename)]<- "blank" 

metadata$substrate[grep("soil", metadata$samplename)] <- "soil" 

#bring in growthhabit info
host <- read.csv("./raw_data_notIncluding_sequenceData/Hosttaxon_lifehistory.csv",
                 stringsAsFactors = F, header = T)

dim(metadata)
metadata <- merge(metadata, host, by.x = "taxon.y", by.y = "taxon", all.x = T)
dim(metadata)
metadata$taxon_final <- metadata$taxon.y

metadata$taxon_final[metadata$taxon.y == "Artemisia tridentate var. vaseyana"] <- "Artemisia tridentata"

metadata$taxon_final[metadata$taxon.y == "Artemisia tridentate"] <- "Artemisia tridentata"
metadata$taxon_final[metadata$taxon.y == "Astragalus kentrophyta cf."] <- 
  "Astragalus kentrophyta"
metadata$taxon_final[metadata$taxon.y == "Astragalus kentrophyta var. tegetarius"] <- 
  "Astragalus kentrophyta"
metadata$taxon_final[metadata$taxon.y == "Eriogonum umbellatum var. majus"] <- 
  "Eriogonum umbellatum"
metadata$taxon_final[metadata$taxon.y == "Lupinus argenteus var. argenteus"] <- 
  "Lupinus argenteus"
metadata$taxon_final[metadata$taxon.y == "Potentilla diversifolia var. diversifolia"] <- 
  "Potentilla diversifolia"
metadata$taxon_final[metadata$taxon.y == "Unknown spruce"] <- 
  "Unknown Spruce"

#load treatment vectors that were used during modeling

treats16s <- read.csv("./processedData/treatments_for_modeling_16s.csv", 
                      stringsAsFactors = F)

treats16s$x <- gsub("(\\d+_\\d+_\\d+E[NP])_.*", "\\1", treats16s$x)

treatsITS <- read.csv("./processedData/treatments_for_modeling_ITS.csv", 
                      stringsAsFactors = F)
treatsITS$x <- gsub("(\\d+_\\d+_\\d+E[NP])_.*", "\\1", treatsITS$x)

table(treats16s == treatsITS) #yay


write.csv(metadata, file = "./processedData/metadata_2.csv")

metadata$samplename[metadata$plant == "1_1_1_1"]
#checking to see if Ha26_ITSoC6 was PCRd twice
# dat <- read.csv("./processedData/otuTables/its90_otuTableCLEAN", 
#                 fill = T, header = T, stringsAsFactors = F)
# metadat <- read.csv("./processedData/metadata.csv", stringsAsFactors = F)
# metadat <- metadat[metadat$locus == "ITS",]
# #extract the samples that might have got PCRd twice
# a <- metadat[grep("26", metadat$plate),]
# a <-  a[grep("ITS0C6", a$midplate),]
# a$combo
# 
# #now find the PCR dupes of these to compare, they end with "ITS_2_primary" instead of ITS_1_primary
# #metadat$combo[grep("2_1_1_1_EN_ITS", metadat$combo)]
# 
# #compare the ITS_1 and ITS_2 reps of these samples to see if the read counts are roughly similar,
# #as they should be if PCR was done the same number of times. 
# 
# #make ITS 2 version of the names, so we can search data for the other reps of the offending samples
# b <- gsub("ITS_1_primary", "ITS_2_primary", a$combo)
# #get rid of the MID sequences
# searches <- gsub("ITS_[ATCG]*_[ATCG]*","",b)
# #extract the b reps from the data using this handy way to grep for multiple things
# bdf <- dat[,grep(paste(searches, collapse="|"), names(dat))]
# 
# adf <- dat[,which(names(dat) %in% a$combo)]
# 
# #not the greatest correlation, but clearly not twice as many reads in one than the other. So
# #it seems unlikley that one of these plates got amplified twice.
# cbind(colSums(adf),
#       colSums(bdf))
# sum(colSums(adf))
# sum(colSums(bdf))



