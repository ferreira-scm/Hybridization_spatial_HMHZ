library(RColorBrewer) # needed for some extra colours in one of the graphs
library(vegan)
library(phyloseq)
library(igraph)
library(Hmisc)
library(Matrix)
library(SpiecEasi)

library(microbiome)

PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

Bac <- subset_taxa(PS.TSS, Kingdom %in%"Bacteria")
Parasite <- subset_taxa(PS.TSS, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))
Parasite@tax_table[,1] <- "Parasite"
Fungi <- subset_taxa(PS.TSS, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota"))
Fungi@tax_table[,1] <- "Fungi"


##### First a preliminary network
PS.f <- merge_phyloseq(Bac, Fungi)
PS.f <- merge_phyloseq(PS.f, Parasite)

PS.f

KeepTaxap <- microbiome::prevalence(PS.f)>0.10
PS.f <- phyloseq::prune_taxa(KeepTaxap, PS.f)

U.cor <- rcorr(as.matrix(PS.f@otu_table), type="spearman")
U.pval <- forceSymmetric(U.cor$P) # Self-correlation as NA
U.p <- p.adjust(U.pval, method="BH")# adjust for multiple testing
U.pval@x <- U.p
p.yes <- U.pval<0.01
r.val = U.cor$r # select all the correlation values
p.yes.r <- r.val*p.yes # only select correlation values based on p-value criterion
## select asv based on rho
yes.r <- abs(p.yes.r)>0.2 # output is logical vector
p.yes.rr <- p.yes.r*yes.r # use logical vector for subscripting.
adjm <- as.matrix(p.yes.rr)
all(taxa_names(PS.f)==colnames(adjm))
adjm[is.na(adjm)] <- 0
net.grph <- graph_from_adjacency_matrix(adjm, mode="undirected", weighted=T)
edgew<-E(net.grph)$weight
E(net.grph)$weight <- abs(E(net.grph)$weight)
V(net.grph)$Genus <- as.vector(PS.f@tax_table[,6])
V(net.grph)$Domain <- as.vector(PS.f@tax_table[,1])
V(net.grph)$Phylum <- as.vector(PS.f@tax_table[,2])
bad.vs<-V(net.grph)[degree(net.grph) == 0]
net.grph <- delete.vertices(net.grph, bad.vs)

# set edge color postive correlation pink color, negative blue.
E.color.Uni = edgew
E.color.Uni = ifelse(E.color.Uni>0, "pink",ifelse(E.color.Uni<0, "blue","grey"))
E(net.grph)$color = as.character(E.color.Uni)
#change edge width
E(net.grph)$width = abs(edgew)

summary(as.factor(V(net.grph)$Domain))

V(net.grph)$Stype <- "shape"
V(net.grph)$Stype[which(V(net.grph)$Domain=="Parasite")] <- "sphere"
V(net.grph)$Stype[which(V(net.grph)$Domain=="Bacteria")] <- "circle"
V(net.grph)$Stype[which(V(net.grph)$Domain=="Fungi")] <- "square"

V(net.grph)$color <- "black"
V(net.grph)$color[which(V(net.grph)$Domain=="Parasite")] <- "deeppink3"
V(net.grph)$color[which(V(net.grph)$Domain=="Bacteria")] <- "darkolivegreen"
V(net.grph)$color[which(V(net.grph)$Domain=="Fungi")] <- "blue"

col <- V(net.grph)$color

degS <- degree(net.grph)

set.seed(123)
plot(net.grph,
     vertex.size=3,
     vertex.shape=V(net.grph)$Stype,
     vertex.frame.color=col,
     vertex.color=col,
          vertex.label="")

############################################################################
#### let's try the spiec easi
Euk <- merge_phyloseq(Parasite, Fungi)
Bac

pargs <- list(rep.num=1000, seed=10010, ncores=90, thresh=0.05)
## mb
#t1 <- Sys.time()
#se.net <- spiec.easi(list(Bac, Euk), method="mb", pulsar.params=pargs)
#t2 <- Sys.time()
#t2-t1
#saveRDS(se.net, "tmp/se.net.rds")

se.net <- readRDS("tmp/se.net.rds")

KeepTaxap <- microbiome::prevalence(Bac)>0.05
Bac.f <- phyloseq::prune_taxa(KeepTaxap, Bac)

KeepTaxap <- microbiome::prevalence(Euk)>0.05
Euk.f <- phyloseq::prune_taxa(KeepTaxap, Euk)

## mb
se.fnet <- spiec.easi(list(Bac.f, Euk.f), method="mb", pulsar.params=pargs)
saveRDS(se.fnet, "tmp/se.fnet.rds")


se.fnet$select$stars$summary # lambda path

# coding bacteria/eukaryote nodes
dtype <- c(rep(1,ntaxa(Bac.f)), rep(2,ntaxa(Euk.f)))

bac.ids=taxa_names(Bac.f)
euk.ids= taxa_names(Euk.f)
net.ids <- c(bac.ids,euk.ids)

# plotting
bm=symBeta(getOptBeta(se.fnet), mode="maxabs")

diag(bm) <- 0

#weights <- Matrix::summary(t(bm))[,3] # includes negative weights
weights <- (1-Matrix::summary(t(bm))[,3])/2 # ort

net <- adj2igraph(Matrix::drop0(getRefit(se.fnet)),
                    edge.attr=list(weight=weights),
                    vertex.attr = list(name=net.ids))

betaMat=as.matrix(symBeta(getOptBeta(se.fnet)))


# we want positive edges to be green and negative to be red
edges <- E(net)
edge.colors=c()
for (e.index in 1:length(edges)){
    adj.nodes=ends(net, edges[e.index])
    xindex=which(net.ids==adj.nodes[1])
    yindex=which(net.ids==adj.nodes[2])
    beta=betaMat[xindex, yindex]
    if (beta>0){
        edge.colors=append(edge.colors, "#1B7837")
    }else if(beta<0){
        edge.colors=append(edge.colors, "#762A83")
    }
}
E(net)$color=edge.colors

### defining attributes
#V(net10)$family=c(PS.f@tax <- table[,4], ARG.f@tax <- table[,3])
#V(net10)$family2=c(PS.f@tax <- table[,4], ARG.f@tax <- table[,6])
V(net)$type=c(rep("Bacteria", length(taxa_names(Bac.f))), rep("Eukaryote", length(taxa_names(Euk.f))))
V(net)$genus=c(Bac.f@tax_table[,6], Euk.f@tax_table[,6])
V(net)$phylum=c(Bac.f@tax_table[,2], Euk.f@tax_table[,2])
V(net)$domain=c(Bac.f@tax_table[,1], Euk.f@tax_table[,1])

V(net)$stype <- c(rep("circle",ntaxa(Bac.f)), rep("square",ntaxa(Euk.f)))


hub.s <- hub_score(net)$vector

V(net)$lab.hub <- ""


V(net)$lab.hub[which(hub.s>0.3)] <- V(net)$genus[which(hub.s>0.3)]

V(net)$label.cex <- 0.5
V(net)$label.dist <- 0

V(net)$label.degree <- pi/2

# we also want the node color to code for phylum
nb.col <- length(levels(as.factor(V(net)$phylum)))
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)
mc <- coul[as.numeric(as.factor(V(net)$phylum))]

plot(net,
     layout=layout_with_fr(net),
     vertex.shape=V(net)$stype,
     vertex.label=V(net)$lab.hub,
     vertex.size=as.integer(cut(hub.s, breaks=10))+2,
     vertex.color=adjustcolor(mc,0.8),
     edge.width=as.integer(cut(E(net)$weight, breaks=6))/3,
     margin=c(0,1,0,0))
legend(x=-2, y=1, legend=levels(as.factor(V(net)$phylum)), col=coul, bty="n",x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)

################### modularity analysis
### modules with fast and greedy
modules =cluster_fast_greedy(net, weights=E(net)$weight)
modules.louvain <- cluster_louvain(net, weights=E(net)$weight)

modularity(modules)

sizes(modules)

V(net)$cluster=modules$membership

## should we plot this?
nodes <- V(net)$name

cluster_id <- V(net)$cluster

gen <- paste(V(net)$domain, V(net)$genus, sep="__")

nodes<-as.data.frame(cbind(nodes, cluster_id, gen))

#nodes$cluster_id <- as.numeric(nodes$cluster_id)
#
#nodes$cluster_id[order(nodes$cluster_id)]

nodes[cluster_id==3,]

nodes <- nodes[nodes$cluster_id<46,]

nodes[grep("Fungi", nodes$gen),]

fam.cl <- as.data.frame(table(nodes$fam, nodes$cluster <- id))

##############some viz of cluster 2
net10c <- delete <- edges(net10b, which(E(net10b)$color=="red"))
V(net10c)$c2 <- ""
V(net10c)$c2[which(V(net10c)$cluster==6)] <- V(net10c)$species[which(V(net10c)$cluster==6)]
V(net10c)$c2
V.cl2 <- V(net10b)$name[which(V(net10b)$cluster==6)]
net.cl2 <- induced <- subgraph(net10b, V.cl2)
nodes[cluster <- id==6,]
degc <- igraph::degree(net.cl2, mode="all")
