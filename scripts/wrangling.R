
rm(list=ls())
#load the two demux keys, edit and combine
demux <- read.csv("raw_data_notIncluding_sequenceData/demultiplexingKeys/novaseq3_demux.csv", header = T, stringsAsFactors = F)
demux$forward_barcode <-toupper(demux$forward_barcode)
demux$reverse_barcode <-toupper(demux$reverse_barcode)
demux$combo <- paste(demux$locus, 
                     demux$forward_barcode, 
                     demux$reverse_barcode, 
                     demux$samplename, sep = "_")
demux2 <- read.csv("raw_data_notIncluding_sequenceData/demultiplexingKeys/NovaSeq2_DemuxJH.csv", header = T, stringsAsFactors = F)
demux2$forward_barcode <-toupper(demux2$forward_barcode)
demux2$reverse_barcode <-toupper(demux2$reverse_barcode)
demux2$combo <- paste(demux2$locus, 
                     demux2$forward_barcode, 
                     demux2$reverse_barcode, 
                     demux2$samplename, sep = "_")
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

#merge with leaf area and size information
leafarea <- read.csv("raw_data_notIncluding_sequenceData/area_mg_leafCount.csv", header=T, stringsAsFactors =F)
dim(leafarea)
leafarea <- leafarea[leafarea$notes != "dupe",]
demux_all_leaf <- merge(demux_all, leafarea, by.x = "plant", by.y = "label", all.x = T, all.y = F)
dim(demux_all_leaf)

locus <- read.csv("raw_data_notIncluding_sequenceData/LocationSize.csv", header=T, stringsAsFactors =F)
dim(locus)
demux_all_leaf_locus <- merge(demux_all_leaf, locus, by.x = "plant", by.y = "label", all.x = T, all.y = F)
dim(demux_all_leaf_locus)

#bring in the multispeq
multi <- read.csv("raw_data_notIncluding_sequenceData/multispeq_data.csv", header=T, stringsAsFactors =F)
multi <- multi[-grep("[a-z]", multi$label),]
dim(multi)
multi <- multi[-grep("dupe",multi$note),]

demux_all_leaf_locus_multi <- merge(demux_all_leaf_locus, multi, by.x = "plant", by.y = "label", 
                                    all.x = T, all.y = F)
#multi[!(multi$label %in% demux_all_leaf_locus$plant),]

dim(demux_all_leaf_locus_multi)

#bring in leaf toughness and water retention
tough <- read.csv("raw_data_notIncluding_sequenceData/Toughness_waterRetention.csv", header=T, stringsAsFactors =F)
dim(tough)

demux_all_leaf_locus_multi_tough <- merge(demux_all_leaf_locus_multi, tough, by.x = "plant", by.y = "sample", 
                                    all.x = T, all.y = F)

dim(demux_all_leaf_locus_multi_tough)


#bring in site data
site <- read.csv("raw_data_notIncluding_sequenceData/siteData.csv", header=T, stringsAsFactors =F)
#make region and site labels for data
demux_all_leaf_locus_multi_tough$region_site <- gsub("(\\d+_\\d+)_\\d+_\\d+","\\1",demux_all_leaf_locus_multi_tough$plant)
#table(demux_all_leaf_locus_multi_tough$region_site )
demux_all_leaf_locus_multi_tough_site <- merge(demux_all_leaf_locus_multi_tough, site, by.x = "region_site", by.y = "label", 
                                          all.x = T, all.y = F)
dim(demux_all_leaf_locus_multi_tough_site)

#bring in plant taxon data
taxa <- read.csv("raw_data_notIncluding_sequenceData/TaxaSampled.csv", header=T, stringsAsFactors =F)
#make region and site labels for data
demux_all_leaf_locus_multi_tough_site$region_site_plant <- gsub("(\\d+_\\d+_\\d+)_\\d+","\\1",demux_all_leaf_locus_multi_tough_site$plant)
demux_all_leaf_locus_multi_tough_site_taxa <- merge(demux_all_leaf_locus_multi_tough_site, taxa, by.x = "region_site_plant", by.y = "label", 
                                               all.x = T, all.y = F)
dim(demux_all_leaf_locus_multi_tough_site_taxa)

demux_all_leaf_locus_multi_tough_site_taxa <- demux_all_leaf_locus_multi_tough_site_taxa[,!names(demux_all_leaf_locus_multi_tough_site_taxa) %in%
  c("region", "site", "species", "indiv.x", "indiv.y" )]


write.csv(demux_all_leaf_locus_multi_tough_site_taxa,
          file = "metadata.csv", row.names = F)


#minor manipulation in Excel
