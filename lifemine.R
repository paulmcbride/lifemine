rm(list=ls())
if(!require(stringr)) { message("ERROR: Requires the 'stringr' package") } 
if(!require(taxize)) { message("ERROR: Requires the 'taxize' package") }

#Not a function yet 

#not yet used
#genes <- c("COX1","cytb","ND2")
cox1.syn <- c("COX1","COI","CO1","COXI","cytochrome c oxidase subunit I ","cytochrome c oxidase subunit 1","cytochrome oxidase subunit 1","cytochrome oxidase subunit I ")
cox2.syn <- c("COX2","COII","COXII","cytochrome c oxidase subunit II","cytochrome c oxidase subunit 2","cytochrome oxidase subunit 2","cytochrome oxidase subunit II ")
nd1.syn <- c("ND1","NAD1", "NADH1","NADH dehydrogenase subunit 1","NADH dehydrogenase subunit I ")
nd2.syn <- c("ND2","NAD2","NADH2","NADH dehydrogenase subunit 2","NADH dehydrogenase subunit II ")
cytb.syn <- c("cytb","COB","cytochrome b","cyt b")

#My arbitrary range of gene lengths
lth <- "600:20000"

#Parulidae is a test taxon here because I know what I'm expecting to find if all goes well
taxon_tsn <- as.numeric(get_tsn("Columbidae"))
all_spp <- itis_downstream(taxon_tsn,downto="Species")
all_gen <- data.frame(genus=unique(all_spp$parentName))
all_gen$spp.count <- lapply(all_gen$genus, function(x) length(all_spp$parentName[all_spp$parentName==x]))
all_spp$cox1 <- 0
all_spp$cox1.length <- 0
all_spp$cox2 <- 0
all_spp$cox2.length <- 0
all_spp$nd1 <- 0
all_spp$nd1.length <- 0
all_spp$nd2 <- 0
all_spp$nd2.length <- 0
all_spp$cytb <- 0
all_spp$cytb.length <- 0

gene.lengths<- lapply(all_spp$taxonName, function(x) get_genes_avail(x, lth, getrelated = FALSE))

for (i in 1:nrow(all_spp)) {
  df.current <- data.frame(gene.lengths[i], stringsAsFactors=F)
  list.cox1 <- unlist(lapply(cox1.syn, function(x) c(as.character(df.current$access_num[str_detect(df.current$genesavail, x)]),df.current$length[str_detect(df.current$genesavail, x)])))
  list.cox1 <- data.frame(accession=list.cox1[nchar(list.cox1)>6],length=list.cox1[nchar(list.cox1)<6],stringsAsFactors=F)
  list.cox2 <- unlist(lapply(cox2.syn, function(x) c(as.character(df.current$access_num[str_detect(df.current$genesavail, x)]),df.current$length[str_detect(df.current$genesavail, x)])))
  list.cox2 <- data.frame(accession=list.cox2[nchar(list.cox2)>6],length=list.cox2[nchar(list.cox2)<6],stringsAsFactors=F)
  list.nd1 <- unlist(lapply(nd1.syn, function(x) c(as.character(df.current$access_num[str_detect(df.current$genesavail, x)]),df.current$length[str_detect(df.current$genesavail, x)])))
  list.nd1 <- data.frame(accession=list.nd1[nchar(list.nd1)>6],length=list.nd1[nchar(list.nd1)<6],stringsAsFactors=F)
  list.nd2 <- unlist(lapply(nd2.syn, function(x) c(as.character(df.current$access_num[str_detect(df.current$genesavail, x)]),df.current$length[str_detect(df.current$genesavail, x)])))
  list.nd2 <- data.frame(accession=list.nd2[nchar(list.nd2)>6],length=list.nd2[nchar(list.nd2)<6],stringsAsFactors=F)
  list.cytb <- unlist(lapply(cytb.syn, function(x) c(as.character(df.current$access_num[str_detect(df.current$genesavail, x)]),df.current$length[str_detect(df.current$genesavail, x)])))
  list.cytb <- data.frame(accession=list.cytb[nchar(list.cytb)>6],length=list.cytb[nchar(list.cytb)<6],stringsAsFactors=F)
  all_spp$cox1[i] <- list.cox1$accession[list.cox1$length==max(as.numeric(list.cox1$length))][1]
  all_spp$cox2[i] <- list.cox2$accession[list.cox2$length==max(as.numeric(list.cox2$length))][1]
  all_spp$nd1[i] <- list.nd1$accession[list.nd1$length==max(as.numeric(list.nd1$length))][1]
  all_spp$nd2[i] <- list.nd2$accession[list.nd2$length==max(as.numeric(list.nd2$length))][1]
  all_spp$cytb[i] <- list.cytb$accession[list.cytb$length==max(as.numeric(list.cytb$length))][1]
  all_spp$cox1.length[i] <- list.cox1$length[list.cox1$length==max(as.numeric(list.cox1$length))][1]
  all_spp$cox2.length[i] <- list.cox2$length[list.cox2$length==max(as.numeric(list.cox2$length))][1]
  all_spp$nd1.length[i] <- list.nd1$length[list.nd1$length==max(as.numeric(list.nd1$length))][1]
  all_spp$nd2.length[i] <- list.nd2$length[list.nd2$length==max(as.numeric(list.nd2$length))][1]
  all_spp$cytb.length[i] <- list.cytb$length[list.cytb$length==max(as.numeric(list.cytb$length))][1]
}
all_gen$cox1 <- 0
all_gen$cox2 <- 0
all_gen$nd1 <- 0
all_gen$nd2 <- 0
all_gen$cytb <- 0

for(j in 1:nrow(all_gen)) {
  all_gen$cox1[j] <-  100*(length(na.omit(all_spp$cox1[all_spp$parentName==all_gen$genus[j]])))/as.numeric(all_gen$spp.count[j])
  all_gen$cox2[j] <-  100*length(na.omit(all_spp$cox2[all_spp$parentName==all_gen$genus[j]]))/as.numeric(all_gen$spp.count[j])
  all_gen$nd1[j] <-  100*length(na.omit(all_spp$nd1[all_spp$parentName==all_gen$genus[j]]))/as.numeric(all_gen$spp.count[j])
  all_gen$nd2[j] <-  100*length(na.omit(all_spp$nd2[all_spp$parentName==all_gen$genus[j]]))/as.numeric(all_gen$spp.count[j])
  all_gen$cytb[j] <-  100*length(na.omit(all_spp$cytb[all_spp$parentName==all_gen$genus[j]]))/as.numeric(all_gen$spp.count[j])
}