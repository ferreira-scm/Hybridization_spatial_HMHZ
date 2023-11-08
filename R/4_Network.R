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


Euk <- merge_phyloseq(Parasite, Fungi)
Bac


## prevalebce filtering of 5%
KeepTaxap <- microbiome::prevalence(Bac)>0.05
Bac <- phyloseq::prune_taxa(KeepTaxap, Bac)
KeepTaxap <- microbiome::prevalence(Euk)>0.05
Euk <- phyloseq::prune_taxa(KeepTaxap, Euk)
#### spiec easi
pargs <- list(rep.num=1000, seed=10010, ncores=90, thresh=0.05)
## mb
#t1 <- Sys.time()
#se.net <- spiec.easi(list(Bac, Euk), method="mb", pulsar.params=pargs)
#t2 <- Sys.time()
#t2-t1
#saveRDS(se.net, "tmp/se.fnet.rds")

se.net <- readRDS("tmp/se.fnet.rds")

se.net$select$stars$summary # lambda path

# coding bacteria/eukaryote nodes
dtype <- c(rep(1,ntaxa(Bac)), rep(2,ntaxa(Euk)))

bac.ids=taxa_names(Bac)
euk.ids= taxa_names(Euk)
net.ids <- c(bac.ids,euk.ids)

# plotting
bm=symBeta(getOptBeta(se.net), mode="maxabs")

diag(bm) <- 0

#weights <- Matrix::summary(t(bm))[,3] # includes negative weights
weights <- (1-Matrix::summary(t(bm))[,3])/2 # ort

weights

net <- SpiecEasi::adj2igraph(Matrix::drop0(getRefit(se.net)),
                    edge.attr=list(weight=weights),
                    vertex.attr = list(name=net.ids))

E(net)

betaMat=as.matrix(symBeta(getOptBeta(se.net)))

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
V(net)$type=c(rep("Bacteria", length(taxa_names(Bac))), rep("Eukaryote", length(taxa_names(Euk))))
V(net)$genus=c(Bac@tax_table[,6], Euk@tax_table[,6])
V(net)$family=c(Bac@tax_table[,5], Euk@tax_table[,5])
V(net)$phylum=c(Bac@tax_table[,2], Euk@tax_table[,2])
V(net)$domain=c(Bac@tax_table[,1], Euk@tax_table[,1])

V(net)$stype <- c(rep("circle",ntaxa(Bac)), rep("square",ntaxa(Euk)))

bad<-V(net)[degree(net) == 0]

net <-delete.vertices(net, bad)

hub.s <- hub_score(net)$vector

V(net)$lab.hub <- ""


V(net)$lab.hub <- V(net)$genus

V(net)$label.cex <- 0.5
V(net)$label.dist <- 0

V(net)$label.degree <- pi/2

# we also want the node color to code for phylum
nb.col <- length(levels(as.factor(V(net)$phylum)))
coul <- colorRampPalette(brewer.pal(8, "Accent"))(nb.col)
mc <- coul[as.numeric(as.factor(V(net)$phylum))]



pdf("fig/Network_prev05.pdf",
                width =10, height = 10)
set.seed(1002)
plot(net,
     layout=layout_with_fr(net),
     vertex.shape=V(net)$stype,
     vertex.label=V(net)$lab.hub,
     vertex.label.dist=0.4,
     vertex.label.degree=-pi/2,
     vertex.size=3,
     vertex.color=adjustcolor(mc,0.8),
     edge.width=as.integer(cut(E(net)$weight, breaks=6))/3,
     margin=c(0,1,0,0))
legend(x=-2, y=1, legend=levels(as.factor(V(net)$phylum)), col=coul, bty="n",x.intersp=0.25,text.width=0.045, pch=20, pt.cex=1.5)
dev.off()


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
