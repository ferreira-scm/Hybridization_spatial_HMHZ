library(ggplot2)
library(reshape)
library(microbiome)
library(vegan)
library(phyloseq)
library(DECIPHER)
library(phangorn)
library(ShortRead)
library(dada2)
library(dplyr)
library("Hmisc", lib.loc="/usr/local/lib/R/site-library/")
source("R/Correlation_net.R")

### this script does a quality filter (by amplicon) and then some transformations
# TSS --> Total sum scaling (per amplicon) and then merges all amplicon datasets into 1 (PS.TSS)
# fPS is the filered dataset with no normalization

# We also remove the silva handlers from the taxonomic table (all the "g__", "s__"...), so this
#makes the script run slower than expected.

# we do a bit of cleaning up in the sampla data, by adding inputed immune gene responses
# and by inputing the few missing values for BMI and HI

#### we clean up parasite ASVs by individual acessing if each parasite genus likely has several species

# we systematically merge ASV's that are likely from the same taxon. 

#output from this script:
# PS.TSS Relative abundances - total sum scaling
# fPS filtered dataset


#### preprocessing: filtering and transforming
PS.l <- readRDS(file = "data/PhyloSeqList_HMHZ_All.Rds")
# this is our filtering function
fil <- function(ps){
    x = phyloseq::taxa_sums(ps)
    # abundance filtering at 0.005%
    keepTaxa = (x / sum(x) > 0.00005)
    summary(keepTaxa)
    ps = phyloseq::prune_taxa(keepTaxa, ps)
# plus prevalnce filter at 1%
    KeepTaxap <- microbiome::prevalence(ps)>0.01
    ps <- phyloseq::prune_taxa(KeepTaxap, ps)
# subset samples based on total read count (100 reads)
#ps <- phyloseq::subset_samples(ps, phyloseq::sample_sums(ps) > 100)
    ps <- phyloseq::prune_samples(sample_sums(ps)>100, ps)
    ps
}
## filtering MA by amplicon
fPS.l <- list()
for (i in 1:length(PS.l)) {
    try(fPS.l[[i]] <- fil(PS.l[[i]]), silent=TRUE)
}
# and remove those handlers.
#test <- list()
for (i in 1:length(fPS.l)) {
    try(tax_table(fPS.l[[i]])[, colnames(tax_table(fPS.l[[i]]))] <- gsub(tax_table(fPS.l[[i]])[, colnames(tax_table(fPS.l[[i]]))], pattern="[a-z]__", replacement=""), silent=TRUE)
}
### now pool all amplicons
fPS <- fPS.l[[1]]
for (i in 2:length(fPS.l)){
    fPS <- try(merge_phyloseq(fPS,fPS.l[[i]]))
#    print(fPS)
}
#### let's transform by amplicon
PS.lTSS <- list()

for (i in 1:length(fPS.l)) {
    try(PS.lTSS[[i]] <- transform_sample_counts(fPS.l[[i]], function(x) x / sum(x)), silent=TRUE)
}
PS.TSS <- PS.lTSS[[1]]# and pool into 1 phyloseq object
for (i in 2:length(PS.lTSS)){
    PS.TSS <- try(merge_phyloseq(PS.TSS,PS.lTSS[[i]]))
}

# Now we want to relabel Eimeria ASV species as in https://github.com/ferreira-scm/Eimeria_AmpSeq/blob/master/R/Wild_8_AssignEim_tree_cor.R
Eim <- readRDS("data/EimeriaSpeciesAssign.RDS")
eim_rn <- which(rownames(PS.TSS@tax_table) %in% rownames(Eim@tax_table))
# sanity check
rownames(PS.TSS@tax_table[eim_rn])== rownames(Eim@tax_table)
fPS@tax_table[eim_rn, 7] <- Eim@tax_table[,7]
PS.TSS@tax_table[eim_rn, 7] <- Eim@tax_table[,7]
## ok now we agglomerate Eimeria species
PS.TSS1 <- subset_taxa(PS.TSS, !Genus%in%"Eimeria")
PS.TSSE <- subset_taxa(PS.TSS, Genus%in%"Eimeria")
PS.TSSE <- tax_glom(PS.TSSE, "Species")
PS.TSS <- merge_phyloseq(PS.TSSE, PS.TSS1)
fPS1 <- subset_taxa(fPS, !Genus%in%"Eimeria")
fPSE <- subset_taxa(fPS, Genus%in%"Eimeria")
fPSE <- tax_glom(fPSE, "Species")
fPS <- merge_phyloseq(fPSE, fPS1)

#### let's adjust the metadata too, to includ immune gene measures that have been inputed
Immune <- read.csv("https://raw.githubusercontent.com/fayweb/Eimeria_mouse_immunity/main/output_data/2.imputed_MICE_data_set.csv")
#subsetting to filed animals only
Immune <- Immune[Immune$origin=="Field",]
# keeping wanted variables only
keep <- c("Mouse_ID","IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")
Nkeep <- c("IFNy", "CXCR3", "IL.6", "IL.13", "IL1RN", "CASP1", "CXCL9", "IDO1", "IRGM1", "MPO", "MUC2", "MUC5AC", "MYD88", "NCR1", "PRF1", "RETNLB", "SOCS1", "TICAM1", "TNF")
Immune <- Immune[,keep]
head(Immune)
# fixing name structure
Immune$Mouse_ID <- gsub("AA", "AA_", Immune$Mouse_ID)
# removing samples that are not in our PS
Immune <- Immune[which(Immune$Mouse_ID%in%PS.TSS@sam_data$Mouse_ID),]
# sanity checj
Immune[match(PS.TSS@sam_data$Mouse_ID, Immune$Mouse_ID),"Mouse_ID"]==PS.TSS@sam_data$Mouse_ID
#and replace
fPS@sam_data[,Nkeep] <- Immune[match(PS.TSS@sam_data$Mouse_ID, Immune$Mouse_ID),Nkeep]
PS.TSS@sam_data[,Nkeep] <- Immune[match(PS.TSS@sam_data$Mouse_ID, Immune$Mouse_ID),Nkeep]

#hybridicity
PS.TSS@sam_data$hi <- 0.5-abs(PS.TSS@sam_data$HI-0.5)
fPS@sam_data$hi <- 0.5-abs(fPS@sam_data$HI-0.5)

##### we are going to input that BMI value missing, and the HI values missing too.
library(mice)
df <- data.frame(PS.TSS@sam_data[,c("Mouse_ID", "HI", "BMI")])
df.t <- mice(df, m=5, maxit=50, meth='pmm', seed=500)
df.t <- complete(df.t,1)
# sanity check
all(PS.TSS@sam_data$Mouse_ID==df.t$Mouse_ID)
PS.TSS@sam_data$BMI <- df.t$BMI
fPS@sam_data$BMI <- df.t$BMI
PS.TSS@sam_data$HI <- df.t$HI
fPS@sam_data$HI <- df.t$HI

#### now we clean up parasite ASVs by individual acessing if each parasite genus likely has several species
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

# now we want to merge ASV's that are likely from the same taxon. We expect this because of the multiple amplicons used that amplify same taxons and even within the same amplicon there could be several ASVs from the same taxon due to e.g. multiple gene copy numbers, sequencing errors... Briefly, we do correlation networs per genus, cluster based on positive, significant correlations and merge ASVs within clusters.
# This is different for parasites because we carefully manually evaluate all ASVs for each parasite genus and in here we do this systematically.
library(phyloseq)
library(Hmisc)
library(Matrix)
library(igraph)

genus <- get_taxa_unique(PS.TSS, "Genus")
# We don't do this for Eimeria
genus <- genus[!genus=="Eimeria"]

for (i in 1:length(genus)){
#mergingASV <- function(PS.T, gen){
    print(genus[i])
    Kaza <- prune_taxa(tax_table(PS.TSS)[,6]%in%genus[i], PS.TSS)
#    Kaza <- prune_samples(sample_sums(Kaza)>0, Kaza)
    kaza <- (Kaza@otu_table)
    tax <- data.frame(Kaza@tax_table)
############ correlation matrix################
    otu.cor <- rcorr(as.matrix(kaza), type="spearman")
# p value
    otu.pval <- forceSymmetric(otu.cor$P)
    cor.p <- p.adjust(otu.pval, method="BH") # adjusting for multiple testing
    otu.pval@x<- cor.p
    p.yes <- otu.pval<0.05 # only significant p values
    r.val = otu.cor$r # select all the correlation values
    p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion
# sanity check
#all(rownames(p.yes.r)==colnames(kaza))
############# network basded on the correlation adjancency matrix
    adjm <- as.matrix(p.yes.r)
#ignoring NAs
    adjm[is.na(adjm)] <- 0
    net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)
### remove negative edges
    net=delete.edges(net.grph, which(E(net.grph)$weight<0)) # here's my condition.
#plot(net,
#     vertex.label="")
    oc <- cluster_fast_greedy(net) # cluster
# and now we merge based on the clustered modules
    group <- list()
    for (i in 1:length(levels(as.factor(oc$membership)))){
        group[[i]] <- oc$names[which(oc$membership==i)]
        PS.TSS <- merge_taxa(PS.TSS, group[[i]])
    }
}

get_taxa_unique(PS.TSS, "Order")

PS.TSS@tax_table[97,5] <- "Unknown_family_in_Oscillospirales"
PS.TSS@tax_table[97,6] <- "Unknown_genus_in_Oscillospirales"

PS.TSS@tax_table[8,2] <- "Unknown_phylum_in_Eukarya"
PS.TSS@tax_table[8,3] <- "Unknown_class_in_Eukarya"
PS.TSS@tax_table[8,4] <- "Unknown_order_in_Eukarya"
PS.TSS@tax_table[8,5] <- "Unknown_family_in_Eukarya"
PS.TSS@tax_table[8,6] <- "Unknown_genus_in_Eukarya"


PS.TSS <- prune_taxa(taxa_sums(PS.TSS) > 0, PS.TSS) 

saveRDS(PS.TSS, "tmp/PS.TSS_filtered.rds")

