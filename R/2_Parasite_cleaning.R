source("R/1_filtering.R")
source("R/Correlation_net.R")

library(DECIPHER)
library(phangorn)
library(ShortRead)
library(dada2)
library(dplyr)
library("Hmisc", lib.loc="/usr/local/lib/R/site-library/")

#### In this script we clean up parasite ASVs by individual acessing if each parasite genus likely has several species
# we do this by looking at co-occurrence correlation networks per genus and see if the ASVs (that come from different amplicons) make a hairball structure, suggesting ASVs are from the same taxa or a modular network. If it's modular, then likely ASV's here represent different species and we need to investigate further with phylogenetic analysis. 

# there's evidence of 2 different taxa annotated as oxyurida. Probably Apiculuris and Syphacia
############# We need to distinguis between Oxyurida ASVs
## create reference sequences for Syhacia and Apiculuris and align them to ASV
Oxy <- subset_taxa(PS.TSS, Genus%in%"Oxyurida")
seqs <- DNAStringSet(getSequences(colnames(Oxy@otu_table)))
parasite <- "Oxyurida"
amp <-  amp_func(PS.lTSS, parasite="Oxyurida")
names(seqs) <- amp_func(PS.lTSS, parasite="Oxyurida")
seqs18 <- seqs[-grep("D3A_", names(seqs))]
seqs28 <- seqs[grep("D3A_", names(seqs))]
#writeFasta(seqs, "tmp/Oxyurida.fa")
## get 18S sequences from Aspiculuris, Syphacia and...
access.aspi18 <- c("EF464551.1", "MH215350.1", "KP338606.1", "KY462827.1", "MT755640.1", "KY462828.1", "MT613322.1")
access.syph18 <- c("EF464553.1", "EF464554.1", "AB629697.1", "KY462829.1", "KY462826.1","OK138900.1", "OK138907.1", "OK138885.1", "OK138886.1", "OK138897.1", "OK138899.1", "OK138904.1", "OK138908.1","EU263105.2", "MT135057.1", "MT135058.1")
# retrieve sequences from NBCI and convert formats and label
aspi18.db <- read.GenBank(access.aspi18)
aspi18.db.ID <- paste(names(aspi18.db), attr(aspi18.db, "species"), sep="_")
aspi18.DB <- aspi18.db %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
syph18.db <- read.GenBank(access.syph18)
syph18.db.ID <- paste(names(syph18.db), attr(syph18.db, "species"), sep="_")
syph18.DB <- syph18.db %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(aspi18.DB) <- aspi18.db.ID
names(syph18.DB) <- syph18.db.ID
outg <- read.GenBank(c("JF934731.1", "FR687850.1", "AB626598.1"))# don't forget an outgroup
outg.ID <- paste(names(outg), attr(outg, "species"), sep="_")
#convert to DNAStringset
outg.db <- outg %>% as.character %>% lapply(.,paste0,collapse="") %>% unlist %>% DNAStringSet
names(outg.db) <- outg.ID
# uncomment if rerunning is necessary
#align18 <-  AlignSeqs(c(outg.db, syph18.DB, aspi18.DB, seqs18), anchor=NA, iterations=20, refinements=20, processors=90)
#writeFasta(align18, "tmp/Oxyurida_tree/oxyurida18S.fa")
#run tree outside R
#~/iqtree-2.2.0-Linux/bin/iqtree2 -s tmp/Oxyurida_tree/oxyurida18S.fa -m MFP -B 5000 -T AUTO
# read tree into R
library(ggtree)
phyloOxy <- ggtree::read.tree("tmp/Oxyurida_tree/oxyurida18S.fa.contree")
phyloO <- root(phyloOxy, outg.ID)
clade.nr <- ggtree::MRCA(phyloO, syph18.db.ID)
clade.nr.aspi <- ggtree::MRCA(phyloO, aspi18.db.ID)
syph <- extract.clade(phyloO, clade.nr)
aspi <- extract.clade(phyloO, clade.nr.aspi)
syphASV <- which(names(seqs)%in%syph$tip.label)
aspiASV <- which(names(seqs)%in%aspi$tip.label)
ampdf <-as.data.frame(amp)
ampdf$species <- "name"
ampdf$species[syphASV] <- "Syphacia"
ampdf$species[aspiASV] <- "Aspiculuris"
ampdf$species[ampdf$species=="name"] <- "28S"
fPS@tax_table[which(fPS@tax_table[,6]=="Oxyurida"),7] <- ampdf$species
PS.TSS@tax_table[which(PS.TSS@tax_table[,6]=="Oxyurida"),7] <- ampdf$species
# now the correlation networks
######## Now we do a correlation network to see the likelihood of ASV within the same genus being from the same species.
## I decided to manually evaluate each parasite genus
parasite <- "Oxyurida"
Correlation_net(PS.lTSS, PS.l, PS.TSS, "Oxyurida") # 2 species
parasite <- "Tritrichomonas"
Correlation_net(PS.lTSS, PS.l, PS.TSS, "Tritrichomonas")# likely 1 species: merge
parasite <- "Cryptosporidium"
Correlation_net(PS.lTSS, PS.l, PS.TSS, "Cryptosporidium")# too sparce but likely 1 sp, merge?
parasite <- "Ascaridida"
Correlation_net(PS.lTSS, PS.l, PS.TSS, "Ascaridida")# likely 1 sp, merge
parasite <- "Spirurida"
Correlation_net(PS.lTSS, PS.l, PS.TSS, "Spirurida")# too sparce, likely 1 sp, merge?
parasite <- "Cyclophyllidea"
Correlation_net(PS.lTSS, PS.l, PS.TSS, "Cyclophyllidea")# 1 species, merge
parasite <- "Trichocephalida"
Correlation_net(PS.lTSS, PS.l, PS.TSS, "Trichocephalida")# 1 species, merge

## merging ASVs based on correlations and for oxyuridae, based on phylogenetic information
Clean_parasites <- function(PS.T){
# Oxyurida
PS.T@tax_table[which(PS.T@tax_table[,7]=="28S"),7]<- "Aspiculuris"
PS.T@tax_table[which(PS.T@tax_table[,7]=="Aspiculuris"),6]<- "Aspiculuris"
PS.T@tax_table[which(PS.T@tax_table[,7]=="Syphacia"),6]<- "Syphacia"
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Syphacia")])
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Aspiculuris")])
sample_data(PS.T)$Syphacia_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Syphacia"))
sample_data(PS.T)$Aspiculuris_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Aspiculuris"))
#cor.test(log(1+sample_data(PS.T)$Aspiculuris_asv), log(1+sample_data(PS.T)$Aspiculuris_sp))
#cor.test(log(1+sample_data(PS.T)$Syphacia_asv), log(1+sample_data(PS.T)$Syphacia_sp))
#
# Tritrichomonas
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Tritrichomonas")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Tritrichomonas"),7] <- "Tritrichomonas_sp"
sample_data(PS.T)$Tritrichomonas_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Tritrichomonas"))
#
# Cryptosporidum
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Cryptosporidium")])
sample_data(PS.T)$Crypto_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Cryptosporidium"))
#
#Ascaridida # heterakis?
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Ascaridida")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Ascaridida"),7] <- "Ascaridida"
sample_data(PS.T)$Ascaridida_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Ascaridida"))
#cor.test(log(1+sample_data(PS.T)$Ascaridida_asv), log(1+sample_data(PS.T)$Heterakis))
#plot(log(1+sample_data(PS.T)$Ascaridida_asv), log(1+sample_data(PS.T)$Heterakis))
#
# Spirurida, likely Mastophorus (high correlation with worm counts)
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Spirurida")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Spirurida"),7] <- "Mastophorus_muris"
PS.T@tax_table[which(PS.T@tax_table[,6]=="Spirurida"),6] <- "Mastophorus"
sample_data(PS.T)$Mastophorus_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Mastophorus"))
#cor.test(log(1+sample_data(PS.T)$Mastophorus_asv), log(1+sample_data(PS.T)$Mastophorus_muris))
#plot(sample_data(PS.T)$Mastophorus_asv, sample_data(PS.T)$Mastophorus_muris)
#
#Trichocephalida, likely Trichuris muris
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Trichocephalida")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Trichocephalida"),7] <- "Trichuris_muri"
PS.T@tax_table[which(PS.T@tax_table[,6]=="Trichocephalida"),6] <- "Trichuris"
sample_data(PS.T)$Trichuris_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Trichuris"))
#cor.test(log(1+sample_data(PS.T)$Trichuris_asv), log(1+sample_data(PS.T)$Trichuris_muri))
#plot(log(1+sample_data(PS.T)$Trichuris_asv), log(1+sample_data(PS.T)$Trichuris_muri))
#
# Cyclophyllidae, likely Hymenolepis, not Catenotaenia_pusilla, nor Taenia
PS.T <- merge_taxa(PS.T, taxa_names(PS.T)[which(PS.T@tax_table[,6]=="Cyclophyllidea")])
PS.T@tax_table[which(PS.T@tax_table[,6]=="Cyclophyllidea"),7] <- "Hymenolepis_sp"
PS.T@tax_table[which(PS.T@tax_table[,6]=="Cyclophyllidea"),6] <- "Hymenolepis"
sample_data(PS.T)$Hymenolepis_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Hymenolepis"))
#cor.test(log(1+sample_data(PS.T)$Hymenolepis_asv), log(1+sample_data(PS.T)$Hymenolepis_sp))
#plot(log(1+sample_data(PS.T)$Hymenolepis_asv), log(1+sample_data(PS.T)$Hymenolepis_sp))
#
# adding eimeria intensity to sample data
sample_data(PS.T)$Eimeria_asv <- sample_sums(subset_taxa(PS.T, Genus %in%"Eimeria"))
sample_data(PS.T)$Eimeria_ferrisi_asv <- sample_sums(subset_taxa(PS.T, Species %in%"ferrisi"))
sample_data(PS.T)$Eimeria_falciformis_asv <- sample_sums(subset_taxa(PS.T, Species %in%"falciformis"))
sample_data(PS.T)$Eimeria_vermiformis_asv <- sample_sums(subset_taxa(PS.T, Species %in%"vermiformis"))
   PS.T
}
#
# because we have different datasets, we need to do the same for fPS and PS.TSS
fPS <- Clean_parasites(fPS)
PS.TSS <- Clean_parasites(PS.TSS)

##########################################################
### let's clean up genus column in the tax table
all(sample_names(fPS)==sample_names(PS.TSS))
tax <- as.data.frame(tax_table(PS.TSS))
tax$Kingdom[is.na(tax$Kingdom)] <- "Unknown_domain"
tax$Kingdom[tax$Kingdom=="unidentified"] <- "Unknown_domain"
tax$Kingdom[!tax$Kingdom%in%c("Bacteria", "Archaea", "Unknown_domain")] <- "Eukarya"
tax[is.na(tax$Genus),]$Genus <- paste0("Unknown_genus_in_",tax[is.na(tax$Genus),]$Family)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Order)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Class)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Phylum)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Kingdom)

tax$Phylum[is.na(tax$Phylum)] <- "Unknown_phylum"
tax$Class[is.na(tax$Class)] <- "Unknown_class"
tax$Order[is.na(tax$Order)] <- "Unknown_order"
tax[which(tax$Genus=="uncultured"),"Genus"] <- paste(tax[which(tax$Genus=="uncultured"),"Order"], tax[which(tax$Genus=="uncultured"),"Genus"], sep="_")

fPS@tax_table <-tax_table(as.matrix(tax))
PS.TSS@tax_table <-tax_table(as.matrix(tax))

#### this takes too long to run, let's save intermediate files
saveRDS(PS.TSS, "tmp/PS.TSS.rds")
saveRDS(fPS, "tmp/fPS.rds")






