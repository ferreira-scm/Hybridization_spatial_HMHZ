#### This script has 2 functions:
# amp_func tells us the amplicon origin for each ASV
# Correlation_net makes correlation networks and colours nodes based on amplicon

amp_func <- function(PS.lT, parasite){
nmOxy <- list()
amp <- list()
for (i in 1:length(PS.lT)) {
    try(p <- subset_taxa(PS.lT[[i]],Genus==parasite), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Genus")
        print(paste(i, "- ", names(PS.l[i]), ": ", nrow(p@tax_table), sep=""))
        nmOxy[i] <- names(PS.l[i])
        amp[[i]] <- paste(rep(names(PS.l[i]), nrow(p@tax_table)), "ASV", seq(1, nrow(p@tax_table), 1), sep="_")
    }
    rm(p)
}

nmOxy <- unlist(nmOxy)
amp <- unlist(amp)
}


Correlation_net <- function(PS.lT, PS.l, PS.T, parasite){
# we want to now from which amplicon a particular ASV came from.
#We use this to colour nodes in our correlation networks
print(paste("amplicon NR - amplicon NAME: number of ASVs produced annotated as", parasite, sep=" "))
nmOxy <- list()
amp <- list()
for (i in 1:length(PS.lT)) {
    try(p <- subset_taxa(PS.lT[[i]],Genus==parasite), silent=TRUE)
    if (exists("p")) {
        a <- get_taxa_unique(p, "Genus")
        print(paste(i, "- ", names(PS.l[i]), ": ", nrow(p@tax_table), sep=""))
        nmOxy[i] <- names(PS.l[i])
        amp[[i]] <- paste(rep(names(PS.l[i]), nrow(p@tax_table)), "ASV", seq(1, nrow(p@tax_table), 1), sep="_")
    }
    rm(p)
}
nmOxy <- unlist(nmOxy)
amp <- unlist(amp)

############## co-occurrences network
##### Oxy, let's try and disentangle the species here
library(Hmisc)
library(Matrix)
library(igraph)
library(RColorBrewer)

Oxy <-subset_taxa(PS.T, Genus %in%parasite)
Oxy <- prune_samples(sample_sums(Oxy)>0, Oxy)
oxy <- (Oxy@otu_table)
tax <- data.frame(Oxy@tax_table)

otu.cor <- rcorr(as.matrix(oxy), type="spearman")
otu.pval <- forceSymmetric(otu.cor$P)
cor.p <- p.adjust(otu.pval, method="BH")
otu.pval@x <- cor.p
sel.tax <- tax[rownames(otu.pval),,drop=FALSE]
#sanity check
#all.equal(rownames(sel.tax), rownames(otu.pval))
p.yes <- otu.pval<0.05
r.val = otu.cor$r # select all the correlation values
p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion
adjm <- as.matrix(p.yes.r)
colnames(adjm) <- Oxy@tax_table[,7]
rownames(adjm) <- Oxy@tax_table[,7]
adjm[is.na(adjm)] <- 0

#we also want the node color to code for amplicon
amp_name <- as.factor(gsub("_ASV_[0-9]", "", amp))
#amp_name <- amp_name[-grep("D3A", amp_name)]

 nb.col <- length(levels(amp_name))
 coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)
 mc <- coul[as.numeric(amp_name)]

 net.grph=graph.adjacency(adjm,mode="undirected",weighted=TRUE,diag=FALSE)
 net.grph=delete.edges(net.grph, which(E(net.grph)$weight<0)) # removing negative edges

### negative correlations
    set.seed(1234)
    plot(net.grph,
     vertex.label=Oxy@tax_table[,7],
     vertex.color=adjustcolor(mc, 0.8),
     frame.col="grey")

# I want to cluster to inform our taxonomic annotation
    oc <- cluster_fast_greedy(net.grph) # cluster
    # and now we merge based on the clustered modules
group <- list()

for (i in 1:length(levels(as.factor(oc$membership)))){
    group[[i]] <- oc$names[which(oc$membership==i)]
#    PS.T <- merge_taxa(PS.T, group[[i]])
}

    cat("\nWe have 3 modules based on optimal cluster algorithm (igraph):\n")
print(group)

}
