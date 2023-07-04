# the point of this script is to merge ASV's that are likely from the same taxon. We expect this because of the multiple amplicons used that amplify same taxons and even within the same amplicon there could be several ASVs from the same taxon due to e.g. multiple gene copy numbers, sequencing errors... Briefly, we do correlation networs per genus, cluster based on positive, significant correlations and merge ASVs within clusters.
# This is different for parasites because we carefully manually evaluate all ASVs for each parasite genus and in here we do this automatically.
library(phyloseq)
library(Hmisc)
library(Matrix)
library(igraph)

#fPS <- readRDS("tmp/fPS.rds")
PS.TSS <- readRDS("tmp/PS.TSS.rds")

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

