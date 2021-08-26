done <- read.table("~/Desktop/teton/coplete", stringsAsFactors = F, fill = T)
head(done)
dat <- read.csv("./processedData/sixteenS_otus_to_analyze.csv")
head(dat)

done$V1 <- gsub(".*(Zotu\\d+_\\w+)", "\\1", done$V1)
donep<- paste(done$V1, done$V2, sep = "")

todo <- paste(dat$X1, dat$X2,"16s.csv", sep = "_")      
todo <- gsub(" ","", todo)
todo <- gsub("_16s","16s", todo)

donep <- gsub("var.\\w+_","", donep)
todo <- gsub("var.\\w+1","1", todo)

donep[-grep("16s.csv", donep)] <- paste(donep[-grep("16s.csv", donep)], "16s.csv", sep = "")
donep[!(donep %in% todo)]

still <- todo[!(todo %in% donep)]

