library(ggplot2)
library(microbiome)
library(vegan)
library(phyloseq)
library("Hmisc", lib.loc="/usr/local/lib/R/site-library/")
library(brms)
library(rstan)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(cowplot)
library(Matrix)
library(igraph)

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

labPS <- readRDS("/SAN/Susanas_den/gitProj/Eimeria_AmpSeq/tmp/Lab/PhyloSeqList_All_Tax_New.Rds")

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
labfPS <- list()
for (i in 1:length(labPS)) {
    try(labfPS[[i]] <- fil(labPS[[i]]), silent=TRUE)
}

# and remove those handlers.
#test <- list()
for (i in 1:length(labfPS)) {
    try(tax_table(labfPS[[i]])[, colnames(tax_table(labfPS[[i]]))] <- gsub(tax_table(labfPS[[i]])[, colnames(tax_table(labfPS[[i]]))], pattern="[a-z]__", replacement=""), silent=TRUE)
}

### now pool all amplicons
l.fPS <- labfPS[[1]]
for (i in 2:length(labfPS)){
    l.fPS <- try(merge_phyloseq(l.fPS,labfPS[[i]]))
    print(l.fPS)
}

#### let's transform by amplicon
l.PS.TSS.l <- list()

for (i in 1:length(labfPS)) {
    try(l.PS.TSS.l[[i]] <- transform_sample_counts(labfPS[[i]], function(x) x / sum(x)), silent=TRUE)
}

l.PS.TSS <- l.PS.TSS.l[[1]]# and pool into 1 phyloseq object
for (i in 2:length(l.PS.TSS.l)){
    l.PS.TSS <- try(merge_phyloseq(l.PS.TSS,l.PS.TSS.l[[i]]))
}

subset_samples(l.PS.TSS, dpi%in%c(0, 6))

### let's clean up genus column in the tax table
tax <- as.data.frame(tax_table(l.PS.TSS))

tax$Kingdom[is.na(tax$Kingdom)] <- "Unknown_domain"

#tax$Kingdom[!tax$Kingdom%in%c("Bacteria", "Archaea", "Unknown_domain")] <- "Eukarya"
tax[is.na(tax$Genus),]$Genus <- paste0("Unknown_genus_in_",tax[is.na(tax$Genus),]$Family)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Order)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Class)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Phylum)
tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Genus<-paste0("Unknown_genus_in_",tax[which(tax$Genus=="Unknown_genus_in_NA"),]$Kingdom)

tax$Phylum[is.na(tax$Phylum)] <- "Unknown_phylum"
tax$Class[is.na(tax$Class)] <- "Unknown_class"
tax$Order[is.na(tax$Order)] <- "Unknown_order"
tax[which(tax$Genus=="uncultured"),"Genus"] <- paste(tax[which(tax$Genus=="uncultured"),"Order"], tax[which(tax$Genus=="uncultured"),"Genus"], sep="_")

l.PS.TSS@tax_table <-tax_table(as.matrix(tax))


# now we want to merge ASV's that are likely from the same taxon. We expect this because of the multiple amplicons used that amplify same taxons and even within the same amplicon there could be several ASVs from the same taxon due to e.g. multiple gene copy numbers, sequencing errors... Briefly, we do correlation networs per genus, cluster based on positive, significant correlations and merge ASVs within clusters.
# This is different for parasites because we carefully manually evaluate all ASVs for each parasite genus and in here we do this systematically.

genus <- get_taxa_unique(l.PS.TSS, "Genus")

for (i in 1:length(genus)){
#mergingASV <- function(PS.T, gen){
    print(genus[i])
    Kaza <- prune_taxa(tax_table(l.PS.TSS)[,6]%in%genus[i], l.PS.TSS)
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
        l.PS.TSS <- merge_taxa(l.PS.TSS, group[[i]])
    }
}

##### sanity check

lab0 <- subset_samples(lab, dpi%in%c(0))
lab6 <- subset_samples(lab, dpi%in%c(6))

table(lab6@sam_data$Genome, lab6@sam_data$Strain)

############################################################## day 0 and day 6
lab <- subset_samples(l.PS.TSS, dpi%in%c(0, 6))

lab@sam_data$Genome[lab@sam_data$Genome%in%c("mus x mus", "mus")] <- "Mmm"
lab@sam_data$Genome[lab@sam_data$Genome%in%c("dom")] <- "Mmd"
lab@sam_data$Genome[lab@sam_data$Genome%in%c("mus x dom", "dom x mus")] <- "Hybrid"

lab@sam_data$Hyb[lab@sam_data$Genome%in%c("Mmm", "Mmd")] <- "Parental"
lab@sam_data$Hyb[lab@sam_data$Genome%in%c("Hybrid")] <- "Hybrid"


lab@sam_data$Infected[lab@sam_data$dpi==0] <- "NO"
lab@sam_data$Infected[lab@sam_data$dpi==6] <- "YES"

############ Decompose
Bac <- subset_taxa(lab, Kingdom %in%"Bacteria")
Parasite <- subset_taxa(lab, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas", "Oxyurida"))
Fungi <- subset_taxa(lab, Kingdom %in% c("Fungi"))
Diet <- subset_taxa(lab, Phylum %in% c("Anthophyta", "Phragmoplastophyta", "Charophyta", "Ochrophyta"))

############### Now we do dyadic models
# how many genera?
length(get_taxa_unique(Bac, "Genus"))-length(grep("Unknown", get_taxa_unique(Bac, "Genus")))

Euk <- subset_taxa(lab, Kingdom %in%c("Viridiplantae", "Fungi", "Chloroplastida", "Animalia", "Eukaryota", "Metazoa", "Alveolata"))

length(get_taxa_unique(Euk, "Genus"))-length(grep("Unknown", get_taxa_unique(Euk, "Genus")))

lab@sam_data$HI[lab@sam_data$Genome=="Mmd"] <- 0
lab@sam_data$HI[lab@sam_data$Genome=="Mmm"] <- 1
lab@sam_data$HI[lab@sam_data$Genome=="Hybrid"] <- 0.5

############# First create dyad data#######################
lab@sam_data$key <- paste(lab@sam_data$EH_ID, lab@sam_data$dpi, sep="_")

key <- data.frame(ID=sample_data(lab)$key)
metadt <- sample_data(lab)
####################
metadt$He <- 2*(metadt$HI)*(1-metadt$HI)


## 1) Jaccard distance for microbiome
JACM <- as.matrix(phyloseq::distance(lab, method="jaccard", type="samples", binary=T))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
jac <- c(as.dist(JACM))

## 2) Aitchison distance for microbiome
AIM <- as.matrix(vegan::vegdist(lab@otu_table, method="aitchison", pseudocount=1))
# transpose Aitchison disssimilary matrix to similarty matrix
AIM <- 1-AIM
ait <- c(as.dist(AIM))

## 3) Subspecies
HIM <- c(dist(metadt$HI))

## 4) Making HE distances
He <- c(dist(metadt$He))

## 5) Making Hx distances
#Create data frame with each sample name (character) and sampling time (numeric)
#Create an empty matrix to fill with distances
hxM<-array(0,c(nrow(metadt),nrow(metadt)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(metadt)){
    for (j in 1:nrow(metadt))
    {hxM[i,j]=abs(metadt$He[i] + metadt$He[j])
    }
}
Hx <- c(as.dist(hxM))

## 6) time of infection
tempm<- c(dist(metadt$dpi))

#Combine these vectors into a data frame
lab.dyad<-data.frame(Jac=jac, Ait=ait, HI=HIM, He=He, Hx=Hx, dpi=tempm)

#Now all we need to do is add the identities of both individuals in each dyad as separate columns into the data frame and exclude self-comparisons (as these are not meaningful).
# extracting Individual-combinations present in the matrices
list<-expand.grid(key$ID, key$ID)
# This created individual-to-same-individual pairs as well. Get rid of these:
list<-list[which(list$Var1!=list$Var2),]
# this still has both quantiles in--> add 'unique' key
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))

# add the names of both individuals participating in each dyad into the data frame
lab.dyad$IDA_s<-list$Var2
lab.dyad$IDB_s<-list$Var1

lab.dyad$IDA <- gsub("_.*", "", lab.dyad$IDA_s)
lab.dyad$IDB <- gsub("_.*", "", lab.dyad$IDB_s)

## 1) Distances for bacteria
jac_bac <- as.matrix(phyloseq::distance(Bac, method="jaccard", type="samples", binary=T))
jac_bac[is.na(jac_bac)] <- 0 # defining those as 0 distances
jac_bac <- 1-jac_bac
jac_bac <- c(as.dist(jac_bac))
ait_bac <- as.matrix(vegan::vegdist(Bac@otu_table, method="aitchison", pseudocount=1))
ait_bac <- 1-ait_bac
ait_bac <- c(as.dist(ait_bac))
################## parasite
jac_para <- as.matrix(phyloseq::distance(Parasite, method="jaccard", type="samples", binary=T))
jac_para[is.na(jac_para)] <- 0 # defining those as 0 distances
jac_para <- 1-jac_para
jac_para <- c(as.dist(jac_para))
ait_para <- as.matrix(vegan::vegdist(Parasite@otu_table, method="aitchison", pseudocount=1))
ait_para <- 1-ait_para
ait_para <- c(as.dist(ait_para))
################## fungi
jac_fun <- as.matrix(phyloseq::distance(Fungi, method="jaccard", type="samples", binary=T))
jac_fun[is.na(jac_fun)] <- 0 # defining those as 0 distances
jac_fun <- 1-jac_fun
jac_fun <- c(as.dist(jac_fun))
ait_fun <- as.matrix(vegan::vegdist(Fungi@otu_table, method="aitchison", pseudocount=1))
ait_fun <- 1-ait_fun
ait_fun <- c(as.dist(ait_fun))
################## plants
jac_pla <- as.matrix(phyloseq::distance(Diet, method="jaccard", type="samples", binary=T))
jac_pla[is.na(jac_pla)] <- 0 # defining those as 0 distances
jac_pla <- 1-jac_pla
jac_pla <- c(as.dist(jac_pla))
ait_pla <- as.matrix(vegan::vegdist(Diet@otu_table, method="aitchison", pseudocount=1))
ait_pla <- 1-ait_pla
ait_pla <- c(as.dist(ait_pla))

# adding to big dataset
lab.dyad$jac_bac <- jac_bac
lab.dyad$ait_bac <- ait_bac
lab.dyad$jac_para <- jac_para
lab.dyad$ait_para <- ait_para
lab.dyad$jac_fun <- jac_fun
lab.dyad$ait_fun <- ait_fun
lab.dyad$jac_pla <- jac_pla
lab.dyad$ait_pla <- ait_pla



# Make sure you have got rid of all self comparisons
lab.dyad<-lab.dyad[which(lab.dyad$IDA!=lab.dyad$IDB),]

nrow(lab.dyad)

#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
    range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("He", "dpi")
for(i in 1:ncol(lab.dyad[,which(colnames(lab.dyad)%in%scalecols)])){
    lab.dyad[,which(colnames(lab.dyad)%in%scalecols)][,i]<-range.use(lab.dyad[,which(colnames(lab.dyad)%in%scalecols)][,i],0,1)
    }


modelJ<-brm(Jac~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "zero_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelA<-brm(Ait~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelJ_bac<-brm(jac_bac~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "zero_one_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelA_bac<-brm(ait_bac~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelJ_fun<-brm(jac_fun~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "zero_one_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelA_fun<-brm(ait_fun~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelJ_para<-brm(jac_para~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "zero_one_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelA_para<-brm(ait_para~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelJ_pla<-brm(jac_pla~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "zero_one_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelA_pla<-brm(ait_pla~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = lab.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)


print(summary(modelJ), digits=3)

print(summary(modelA), digits=3)

print(summary(modelJ_fun), digits=3)

print(summary(modelA_fun), digits=3)

print(summary(modelJ_pla), digits=3)

print(summary(modelA_pla), digits=3)

print(summary(modelJ_para), digits=3)

print(summary(modelA_para), digits=3)

print(summary(modelJ_bac), digits=3)

print(summary(modelA_bac), digits=3)


####### Figure 2 #########################
resdf.fun<- function(model1_para, name, ASV){
    para <- summary(model1_para)$fixed
    data.frame(Domain=name,
               ASVs=ASV,
               HI_Estimate=para[rownames(para)=="HI", "Estimate"],
               HI_lCI=para[rownames(para)=="HI", "l-95% CI"],
               HI_uCI=para[rownames(para)=="HI", "u-95% CI"],
               He_Estimate=para[rownames(para)=="He", "Estimate"],
               He_lCI=para[rownames(para)=="He", "l-95% CI"],
               He_uCI=para[rownames(para)=="He", "u-95% CI"],
               dpi_Estimate=para$Estimate[rownames(para)=="dpi"],
               dpi_lCI=para[rownames(para)=="dpi", "l-95% CI"],
               dpi_uCI=para[rownames(para)=="dpi", "u-95% CI"]
               )
}

res.df <-resdf.fun(modelJ_para, "Parasite", 6)
res.df <- rbind(res.df, resdf.fun(modelJ_bac, "Bacteria", 207))
res.df <- rbind(res.df, resdf.fun(modelJ_pla, "Diet", 30))
res.df <- rbind(res.df, resdf.fun(modelJ_fun, "Fungi", 29))
res.df <- rbind(res.df, resdf.fun(modelJ, "Full model", 318))
res.df$Domain <- factor(res.df$Domain, level=c( "Diet", "Bacteria","Parasite", "Fungi", "Full model"))

res.dfA <-resdf.fun(modelA_para, "Parasite", 6)
res.dfA <- rbind(res.dfA, resdf.fun(modelA_bac, "Bacteria", 207))
res.dfA <- rbind(res.dfA, resdf.fun(modelA_pla, "Diet", 30))
res.dfA <- rbind(res.dfA, resdf.fun(modelA_fun, "Fungi", 29))
res.dfA <- rbind(res.dfA, resdf.fun(modelA, "Full model", 318))
res.dfA$Domain <- factor(res.dfA$Domain, level=c( "Diet", "Bacteria", "Parasite", "Fungi", "Full model"))

coul <- c("#136f63", "#032b43", "#3f88c5", "#ffba08", "#d00000")

genJ <- ggplot(res.df, aes(x=HI_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=HI_lCI, xmax=HI_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Genetic distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

genA <- ggplot(res.dfA, aes(x=HI_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=HI_lCI, xmax=HI_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Genetic distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HeJ <- ggplot(res.df, aes(x=He_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=He_lCI, xmax=He_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="hHe-dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HeA <- ggplot(res.dfA, aes(x=He_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=He_lCI, xmax=He_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="hHe-dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")


dpiJ <- ggplot(res.df, aes(x=dpi_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=dpi_lCI, xmax=dpi_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Experimental infection", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

dpiA <- ggplot(res.dfA, aes(x=dpi_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=dpi_lCI, xmax=dpi_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Experimental infection", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")


Fig5 <- plot_grid(genJ, genA, HeJ, HeA, dpiJ, dpiA,
                  labels="auto", ncol=2)


ggsave("fig/figure5.pdf", Fig5, width=170, height=100, units="mm", dpi=300)
