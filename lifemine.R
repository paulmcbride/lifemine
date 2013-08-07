rm(list=ls())
require(stringr)
require(taxize)

#Not a function yet 

#not yet used
#genes <- c("COX1","cytb")

#Parulidae is a test taxon here because I know what I'm expecting to find if all goes well
taxon_tsn <- as.numeric(get_tsn("Parulidae"))
all_spp <- itis_downstream(taxon_tsn,downto="Species")
all_gen <- data.frame(genus=unique(all_spp$parentName))
all_gen$spp.count <- lapply(all_gen$genus, function(x) length(all_spp$parentName[all_spp$parentName==x]))

gene.lengths<- lapply(all_spp$taxonName, function(x) get_genes_avail(x, "700:20000", getrelated = FALSE))

df1 <- data.frame(gene.lengths[1])
