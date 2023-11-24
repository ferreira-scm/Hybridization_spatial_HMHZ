source("R/5_lab.R")

head(lab.dyad)



summary(dropmodel)

library("MCMCglmm")
library(doParallel)
library(parallel)
library(magrittr)
library(dplyr)
library(purrr)
library(forcats)



tax <- as.data.frame(Fungi@tax_table)

tax_G <- as.data.frame(table(tax[,6]))

genuses<- as.character(tax_G$Var1)
genuses <- c("none", genuses)

genuses

i=5

data.dyad_REAL <- lab.dyad


dropRes_list<-list()
cl <- parallel::makeCluster(60, type="FORK")
doParallel::registerDoParallel(cl)

dropRes_list<-foreach(i = 1:length(genuses)) %dopar% {

    fix.X11()
    library(phyloseq)
    #choose the genus to drop and prune the taxa to keep everything else
    gen.i<-genuses[i]
    taxa_tokeep<-rownames(tax[which(tax$Genus!=genuses[i]),])
    mic.i<-prune_taxa(taxa_tokeep, Fungi)
    #Calculate Jaccard and Aitchison microbiome dissimilarity for each mouse pair

    JACM.i <- as.matrix(phyloseq::distance(mic.i, method="jaccard", type="samples", binary=T))
    # transpose Jaccard disssimilary matrix to Jaccard similarty matrix
    JACM.i <- 1-JACM.i
    jac <- c(as.dist(JACM.i))

    ## 2) Aitchison distance for microbiome
    AIM.i <- as.matrix(vegan::vegdist(mic.i@otu_table, method="aitchison", pseudocount=1))
    # transpose Aitchison disssimilary matrix to similarty matrix
    AIM.i <- 1-AIM.i
    ait <- c(as.dist(AIM.i))

    
    #Make a new dyadic data frame from these vectors and order it to be in the same order as the original dyadic data frame
    data.dyad.i<-data.frame(Jaccard=jac,Aitchisons=ait)
    # extracting Sample_name-combinations of the matrix
    list<-expand.grid(key$ID,key$ID)
    # This created sample-to-same-sample pairs as well. Get rid of these:
    list<-list[which(list$Var1!=list$Var2),]
    # the resulting list still has both quantiles of the original matrix (i.e. all values are doubled) in--> add 'unique' key and subset to one quantile only
    list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
    list<-subset(list, !duplicated(list$key))
    # Sample_name combinations are now in the same order as the lower quantile value vector
    # So we can add dyad name and each participant ID to dyadic dataframe
    data.dyad.i$IDA_s<-list$Var2
    data.dyad.i$IDB_s<-list$Var1

    data.dyad.i$IDA <- gsub("_.*", "", data.dyad.i$IDA_s)
    data.dyad.i$IDB <- gsub("_.*", "", data.dyad.i$IDB_s)

    
    # Make sure you have no self comparisons in the data (This is the case by default here, since we are using just one sample per individual)
    data.dyad.i<-data.dyad.i[which(data.dyad.i$IDA!=data.dyad.i$IDB),]

### Combine new Jaccard variable with rest of dyadic data columns
    data.dyad<-data.dyad_REAL
    data.dyad$Jaccard_similarity<-data.dyad.i$Jaccard
    data.dyad$Aitchison_similarity<-data.dyad.i$Aitchison
    
    #factorize terms used for multimembership random structure and make sure levels are     same and in same order
    data.dyad$IDA<-as.factor(data.dyad$IDA)
    data.dyad$IDB<-as.factor(data.dyad$IDB)
    #all(levels(data.dyad$IDA)==levels(data.dyad$IDB))

    # Scale all predictors not between 0-1 already to be between 0-1
    scalecols<-c("He", "dpi")
    range.use <- function(x,min.use,max.use){(x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
    for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
        data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }

    #The MCMCglmm model

    
    dropmodel<-brm(Aitchison_similarity~1+ HI+He+dpi+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0,
                silent = 2,
                refresh=0)

    ASVs_dropped.i<-nrow(tax_table(lab))-nrow(tax_table(mic.i))
    para <- summary(dropmodel)$fixed

    resdf.i<-data.frame(Genus_dropped=gen.i,
                        ASVs_dropped=ASVs_dropped.i,
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
    return(resdf.i)
}

parallel::stopCluster(cl)

dropRes_list
