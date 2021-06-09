
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
dim(demux_all)
#table(demux_all$plant)[which(table(demux_all$plant) > 8)]

#merge with leaf area and size information
leafarea <- read.csv("raw_data_notIncluding_sequenceData/area_mg_leafCount.csv", 
                     header=T, stringsAsFactors =F)
dim(leafarea)
leafarea <- leafarea[leafarea$notes != "dupe",]
demux_all_leaf <- merge(demux_all, leafarea, by.x = "plant", by.y = "label", all = T)
dim(demux_all_leaf)

locus <- read.csv("raw_data_notIncluding_sequenceData/LocationSize.csv", header=T, stringsAsFactors =F)
dim(locus)
demux_all_leaf_locus <- merge(demux_all_leaf, locus, by.x = "plant", by.y = "label", all = T)
dim(demux_all_leaf_locus)

#bring in the multispeq
multi <- read.csv("raw_data_notIncluding_sequenceData/multispeq_data.csv", header=T, stringsAsFactors =F)
multi <- multi[-grep("[a-z]", multi$label),]
dim(multi)
multi <- multi[-grep("dupe",multi$note),]

demux_all_leaf_locus_multi <- merge(demux_all_leaf_locus, multi, 
                                    by.x = "plant", by.y = "label", 
                                    all = T)
#multi[!(multi$label %in% demux_all_leaf_locus$plant),]

dim(demux_all_leaf_locus_multi)

#bring in leaf toughness and water retention
tough <- read.csv("raw_data_notIncluding_sequenceData/Toughness_waterRetention.csv",
                  header=T, stringsAsFactors =F)
dim(tough)

demux_all_leaf_locus_multi_tough <- merge(demux_all_leaf_locus_multi, tough,
                                          by.x = "plant", by.y = "sample", 
                                    all = T)

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
                                          all = T)
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
                                               all = T)
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

#fix misassignment of soil samples to plants
demux_all_leaf_locus_multi_tough_site_taxa$substrate[grep("[MO]$",
      demux_all_leaf_locus_multi_tough_site_taxa$samplename)] <- "soil"

#"blank" substrate doublecheck. Seems good
#demux_all_leaf_locus_multi_tough_site_taxa$samplename[grep("blank",
#                                                           demux_all_leaf_locus_multi_tough_site_taxa$samplename)]

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

#get rid of Ast_Alta samples
metadata <- metadata[-grep("Ast_Alta",metadata$region_site),]

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

#for samples that don't have locus in the name add that back in there
metadata$samplename[grep("_[ENP]*$", metadata$samplename)] <-
      paste(metadata$samplename[grep("_[ENP]*$", metadata$samplename)], "_",
  metadata$locus[grep("_[ENP]*$", metadata$samplename)], sep = "")

#clean up
metadata$samplename[grep("1_3_2_2", metadata$samplename)]
if(metadata$samplename[1210] == "1_3_2_2"){
  metadata$samplename[1210]  <- "1_3_2_2_EN_16S"
}
if(metadata$samplename[1219] == "1_3_2_2"){
  metadata$samplename[1219]  <- "1_3_2_2_EN_ITS"
}
if(metadata$samplename[1221] == "1_3_2_2"){
  metadata$samplename[1221]  <- "1_3_2_2_EN_ITS"
}
if(metadata$samplename[1246] == "1_3_2_2"){
  metadata$samplename[1246]  <- "1_3_2_2_EN_ITS"
}

metadata$samplename[151] <- "1_1_3_9_EN_ITS"
metadata$samplename[162] <- "1_1_3_9_EN_16S"
metadata$samplename[165] <- "1_1_3_9_EN_16S"
metadata$samplename[202] <- "1_1_3_9_EN_ITS"

metadata$samplename[562] <- "1_2_1_4_EN_ITS"
metadata$samplename[568] <- "1_2_1_4_EN_16S"
metadata$samplename[580] <- "1_2_1_4_EN_16S"
metadata$samplename[620] <- "1_2_1_4_EN_ITS"

metadata$samplename[1024] <- "1_2_7_6_EN_ITS"
metadata$samplename[1046] <- "1_2_7_6_EN_16S"
metadata$samplename[1069] <- "1_2_7_6_EN_ITS"
metadata$samplename[1085] <- "1_2_7_6_EN_16S"

metadata$samplename[1106] <- "1_3_1_7_EN_ITS"
metadata$samplename[1129] <- "1_3_1_7_EN_16S"
metadata$samplename[1143] <- "1_3_1_7_EN_ITS"
metadata$samplename[1174] <- "1_3_1_7_EN_16S"

metadata$samplename[1425] <- "1_3_5_7_EN_ITS"
metadata$samplename[1427] <- "1_3_5_7_EN_16S"
metadata$samplename[1468] <- "1_3_5_7_EN_16S"
metadata$samplename[1497] <- "1_3_5_7_EN_ITS"

metadata$samplename[2452] <- "2_2_3_10_EN_16S"
metadata$samplename[2467] <- "2_2_3_10_EN_ITS"
metadata$samplename[2470] <- "2_2_3_10_EN_ITS"
metadata$samplename[2510] <- "2_2_3_10_EN_16S"

metadata$samplename[2597] <- "2_2_5_1_EN_ITS"
metadata$samplename[2604] <- "2_2_5_1_EN_ITS"
metadata$samplename[2607] <- "2_2_5_1_EN_16S"
metadata$samplename[2651] <- "2_2_5_1_EN_16S"

metadata$samplename[2526] <- "2_2_5_1_EN_ITS"
metadata$samplename[2558] <- "2_2_5_1_EN_16S"
metadata$samplename[2574] <- "2_2_5_1_EN_ITS"
metadata$samplename[2587] <- "2_2_5_1_EN_16S"

metadata$samplename[3494] <- "3_1_2_9_EN_ITS"
metadata$samplename[3540] <- "3_1_2_9_EN_16S"
metadata$samplename[3543] <- "3_1_2_9_EN_16S"
metadata$samplename[3564] <- "3_1_2_9_EN_ITS"

metadata$samplename[4959] <- "3_3_5_8_EN_ITS"
metadata$samplename[4958] <- "3_3_5_8_EN_16S"
metadata$samplename[4970] <- "3_3_5_8_EN_ITS"
metadata$samplename[5011] <- "3_3_5_8_EN_16S"

metadata$samplename[4102] <- "3_2_2_5_EN_16S"
metadata$samplename[4130] <- "3_2_2_5_EN_16S"
metadata$samplename[4144] <- "3_2_2_5_EN_ITS"
metadata$samplename[4151] <- "3_2_2_5_EN_ITS"

metadata$samplename[4481] <- "3_2_7_2_EN_16S"
metadata$samplename[4486] <- "3_2_7_2_EN_ITS"
metadata$samplename[4527] <- "3_2_7_2_EN_16S"
metadata$samplename[4546] <- "3_2_7_2_EN_ITS"

metadata$samplename[6454] <- "4_3_1_7_EN_16S"
metadata$samplename[6494] <- "4_3_1_7_EN_ITS"
metadata$samplename[6514] <- "4_3_1_7_EN_16S"
metadata$samplename[6515] <- "4_3_1_7_EN_ITS"

metadata$samplename[8742] <- "6_1_4_6_EN_ITS"
metadata$samplename[8747] <- "6_1_4_6_EN_ITS"
metadata$samplename[8750] <- "6_1_4_6_EN_16S"
metadata$samplename[8753] <- "6_1_4_6_EN_16S"

metadata$samplename[9903] <- "7_1_1_4_EN_ITS"
metadata$samplename[9904] <- "7_1_1_4_EN_ITS"
metadata$samplename[9905] <- "7_1_1_4_EN_16S"
metadata$samplename[9908] <- "7_1_1_4_EN_16S"

metadata$samplename[9157] <- "6_2_5_2_EN_16S"
metadata$samplename[9158] <- "6_2_5_2_EN_16S"
metadata$samplename[9175] <- "6_2_5_2_EN_ITS"
metadata$samplename[9177] <- "6_2_5_2_EN_ITS"

metadata$samplename[9227] <- "6_2_6_5_EN_ITS"
metadata$samplename[9296] <- "6_2_6_5_EN_16S"
metadata$samplename[9297] <- "6_2_6_5_EN_ITS"
metadata$samplename[9303] <- "6_2_6_5_EN_16S"

metadata$samplename[9202] <- "6_2_5_8_EN_ITS"
metadata$samplename[9203] <- "6_2_5_8_EN_ITS"
metadata$samplename[9205] <- "6_2_5_8_EN_16S"
metadata$samplename[9208] <- "6_2_5_8_EN_16S"

metadata$samplename[9385] <- "6_3_1_10_EN_ITS"
metadata$samplename[9395] <- "6_3_1_10_EN_ITS"
metadata$samplename[9396] <- "6_3_1_10_EN_16S"
metadata$samplename[9399] <- "6_3_1_10_EN_16S"

metadata$samplename[9517] <- "6_3_2_7_EN_ITS"
metadata$samplename[9521] <- "6_3_2_7_EN_16S"
metadata$samplename[9523] <- "6_3_2_7_EN_ITS"
metadata$samplename[9524] <- "6_3_2_7_EN_16S"

metadata$samplename[9418] <- "6_3_1_4_EN_ITS"
metadata$samplename[9451] <- "6_3_1_4_EN_16S"
metadata$samplename[9453] <- "6_3_1_4_EN_ITS"
metadata$samplename[9454] <- "6_3_1_4_EN_16S"

metadata$samplename[9710] <- "6_3_5_10_EN_16S"
metadata$plant[9710] <- "6_3_5_10"

metadata$samplename[9719] <- "6_3_5_10_EN_ITS"
metadata$plant[9719] <- "6_3_5_10"

metadata$samplename[9720] <- "6_3_5_10_EN_16S"
metadata$plant[9720] <- "6_3_5_10"

metadata$samplename[9721] <- "6_3_5_10_EN_ITS"
metadata$plant[9721] <- "6_3_5_10"

metadata$samplename[10407] <- "7_2_2_3_EN_ITS"
metadata$samplename[10411] <- "7_2_2_3_EN_16S"
metadata$samplename[10412] <- "7_2_2_3_EN_ITS"
metadata$samplename[10414] <- "7_2_2_3_EN_16S"

metadata$samplename[10359] <- "7_2_2_7_EN_ITS"
metadata$samplename[10360] <- "7_2_2_7_EN_16S"
metadata$samplename[10362] <- "7_2_2_7_EN_16S"
metadata$samplename[10364] <- "7_2_2_7_EN_ITS"

metadata$samplename[10298] <- "7_2_1_5_EN_16S"
metadata$samplename[10299] <- "7_2_1_5_EN_ITS"
metadata$samplename[10302] <- "7_2_1_5_EN_16S"
metadata$samplename[10304] <- "7_2_1_5_EN_ITS"

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
metadata$compartment[metadata$samplename =="2_1_4_6"] <- "EN"
metadata$samplename[metadata$samplename =="2_1_4_6"] <-
  paste(metadata$samplename[metadata$samplename =="2_1_4_6"],
       "EN", metadata$locus[metadata$samplename =="2_1_4_6"], sep ="_" )

#Checked the otu table names all the way through to 2_1_6_9
#and all of them seemed fine.

#metadata$samplename[grep("6_3_1_3",metadata$samplename)]
metadata$samplename[metadata$samplename =="6_3_1_3"] <-
  paste(metadata$samplename[metadata$samplename =="6_3_1_3"],
        "EN", metadata$locus[metadata$samplename =="6_3_1_3"], sep ="_" )

#metadata$samplename[grep("3_3_1_9",metadata$samplename)]
metadata$samplename[metadata$samplename =="3_3_1_9"] <-
  paste(metadata$samplename[metadata$samplename =="3_3_1_9"],
        "EN", metadata$locus[metadata$samplename =="3_3_1_9"], sep ="_" )

#metadata$samplename[grep("3_3_3_10",metadata$samplename)]
metadata$samplename[metadata$samplename =="3_3_3_10"] <-
  paste(metadata$samplename[metadata$samplename =="3_3_3_10"],
        "EN", metadata$locus[metadata$samplename =="3_3_3_10"], sep ="_" )

metadata$samplename[metadata$samplename =="3_2_7_4"] <-
  paste(metadata$samplename[metadata$samplename =="3_2_7_4"],
        "EN", metadata$locus[metadata$samplename =="3_2_7_4"], sep ="_" )

metadata$samplename[metadata$samplename =="5_2_7_6"] <-
  paste(metadata$samplename[metadata$samplename =="5_2_7_6"],
        "EN", metadata$locus[metadata$samplename =="5_2_7_6"], sep ="_" )

metadata$substrate[grep("48_6_3_24_37", metadata$samplename)] <- "soil"

metadata$substrate[grep("lank", metadata$samplename)]<- "blank" 

metadata$substrate[grep("soil", metadata$samplename)] <- "soil" 

#bring in growthhabit info
host <- read.csv("./raw_data_notIncluding_sequenceData/Hosttaxon_lifehistory.csv",
                 stringsAsFactors = F, header = T)

dim(metadata)
metadata <- merge(metadata, host, by.x = "taxon.x", by.y = "taxon", all.x = T)
dim(metadata)
metadata$taxon_final <- metadata$taxon.x
metadata$taxon_final[metadata$taxon.x == "Artemisia tridentate var. vaseyana"] <- "Artemisia tridentata"

metadata$taxon_final[metadata$taxon.x == "Artemisia tridentate"] <- "Artemisia tridentate"
metadata$taxon_final[metadata$taxon.x == "Astragalus kentrophyta cf."] <- 
  "Astragalus kentrophyta"
metadata$taxon_final[metadata$taxon.x == "Astragalus kentrophyta var. tegetarius"] <- 
  "Astragalus kentrophyta"
metadata$taxon_final[metadata$taxon.x == "Eriogonum umbellatum var. majus"] <- 
  "Eriogonum umbellatum"
metadata$taxon_final[metadata$taxon.x == "Lupinus argenteus var. argenteus"] <- 
  "Lupinus argenteus"
metadata$taxon_final[metadata$taxon.x == "Potentilla diversifolia var. diversifolia"] <- 
  "Potentilla diversifolia"

#load treatment vectors that were used during modeling

treats16s <- read.csv("./processedData/treatments_for_modeling_16s.csv", 
                      stringsAsFactors = F)
treats16s$x <- gsub("(\\d+_\\d+_\\d+E[NP])_.*", "\\1", treats16s$x)

treatsITS <- read.csv("./processedData/treatments_for_modeling_ITS.csv", 
                      stringsAsFactors = F)
treatsITS$x <- gsub("(\\d+_\\d+_\\d+E[NP])_.*", "\\1", treatsITS$x)

table(treats16s == treatsITS) #yay

flora <- read.csv(file = "./processedData/shannons_siteLevel_veg.csv",
                  stringsAsFactors = F, header = T)
flora$newdat.dat.siteLabel <- gsub("(\\d+_\\d+).*", "\\1", flora$newdat.dat.siteLabel)
flora$newdat.dat.siteLabel <- gsub(".*(\\d+_\\d+)", "\\1", flora$newdat.dat.siteLabel)

#dat <- read.csv("processedData/treatments_metadata.csv", stringsAsFactors = F)

newdat <- merge(metadata, flora, by.y = "newdat.dat.siteLabel",
                by.x = "region_site", all.x = T)
names(newdat)
dim(newdat)
dim(metadata)
names(newdat)[length(newdat)] <- "shannons_flora"

#Densitometer
denso <- read.csv(file = "./processedData/densitometer_mean_vs_site.csv",
                  stringsAsFactors = F, header = T)
denso$dat.siteLabel <- gsub("(\\d+_\\d+).*", "\\1", denso$dat.siteLabel)
denso$dat.siteLabel <- gsub(".*(\\d+_\\d+)", "\\1", denso$dat.siteLabel)

newdat <- merge(y = denso, x = newdat, 
                by.x = "region_site",
                by.y = "dat.siteLabel", all.x = T)
names(newdat)[length(newdat)] <- "densitometer"
dim(newdat)

#latlong
latlong <- read.csv(file = "./processedData/latlong_elevation_by_site.csv",
                    stringsAsFactors = F, header = T)
latlong$unique.dat.siteLabel. <- gsub("(\\d+_\\d+).*", "\\1", latlong$unique.dat.siteLabel.)
latlong$unique.dat.siteLabel. <- gsub(".*(\\d+_\\d+)", "\\1", latlong$unique.dat.siteLabel.)


newdat <- merge(y = latlong, x = newdat, 
                by.x = "region_site",
                by.y = "unique.dat.siteLabel.", all.x = T)
names(newdat)[length(newdat)] <- "latlong"
dim(newdat)

mems <- read.csv("processedData/moransEigenVectors.csv", stringsAsFactors = F)
newdat <- merge(y = mems, x = newdat, 
                by.x = "region_site",
                by.y = "unique.dat.label.", all.x = T)

write.csv(newdat, file = "./processedData/metadata_2.csv")

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



