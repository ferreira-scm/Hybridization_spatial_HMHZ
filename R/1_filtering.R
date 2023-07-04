library(ggplot2)
library(reshape)
library(microbiome)
library(vegan)
library(phyloseq)

### this script does a quality filter (by amplicon) and then some transformations
# TSS --> Total sum scaling (per amplicon) and then merges all amplicon datasets into 1 (PS.TSS)
# fPS is the filered dataset with no normalization

# We also remove the silva handlers from the taxonomic table (all the "g__", "s__"...), so this
#makes the script run slower than expected.

# finally, we do a bit of cleaning up in the sampla data, by adding inputed immune gene responses
# and by inputing the few missing values for BMI and HI

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

