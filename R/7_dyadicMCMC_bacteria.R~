# This script follows the super cool approach from Aura Raulo
# https://github.com/nuorenarra/Analysing-dyadic-data-with-brms/blob/main/R_Making_dyadic_data/DYADIC_workshop_data_wrangling.Rmd

library("MCMCglmm")
library(ape)
library(brms)
library(rstan)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(ggmcmc)
library(ggthemes)
library(ggridges)
library(vegan)
library(phyloseq)
library(ggplot2)
library(bayesplot)
library(bayestestR)
library(brms)
library(MCMCglmm)
library(doParallel)
library(parallel)
library(magrittr)
library(dplyr)
library(purrr)
library(forcats)
library(tidyr)
library(modelr)
library(ggdist)
library(tidybayes)
library(ggplot2)
library(cowplot)
library(rstan)
library(brms)
library(ggrepel)
library(RColorBrewer)
library(gganimate)
library(posterior)
library(distributional)


### we don't include sex because it does not converge
PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

Bac <- subset_taxa(PS.TSS, Kingdom %in%"Bacteria")

Bac

############# First create dyad data#######################
key <- data.frame(ID=sample_data(Bac)$Mouse_ID)
metadt <- sample_data(Bac)

####################
## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(Bac, method="jaccard", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
dimnames(JACM)<- c(key, key)

## 2) Spatial distance matrix
distance.df <- metadt[,c("Mouse_ID", "Longitude", "Latitude")]
SPATM <- array(NA, c(length(distance.df$Mouse_ID),length(distance.df$Mouse_ID)))
# derive matrix with spatial distances between each location
for (i in 1:length(distance.df$Mouse_ID)){
    for (j in 1:length(distance.df$Mouse_ID))
    {SPATM[i,j]= sqrt((abs(distance.df$Longitude[i]-distance.df$Longitude[j]))^2+
                      (abs(distance.df$Latitude[i]-distance.df$Latitude[j]))^2)
    }
}
dimnames(SPATM)<- c(key, key)


## 1) Chisq distance
CHIM <- as.matrix(vegan::vegdist(Bac@otu_table, method="chisq"))

# transpose Chi square disssimilary matrix to similarty matrix
CHIM <- 1-CHIM
# sanity check
all(rownames(CHIM)==key)
dimnames(CHIM)<- c(key, key)

# 3) pairwise genetic distance based on genetic data
# I actually need to get this from sota again
sota <- read.csv("https://raw.githubusercontent.com/derele/Mouse_Eimeria_Field/master/data_products/SOTA_Data_Product.csv")
gen <- c("mtBamH", "YNPAR", "X332", "X347", "X65", "Tsx", "Btk", "Syap1",
         "Es1","Gpd1","Idh1","Mpi","Np", "Sod1", "Es1C", "Gpd1C", "Idh1C",
         "MpiC","NpC", "Sod1C")
sota <- sota[sota$Mouse_ID%in%distance.df$Mouse_ID,c("Mouse_ID", gen)]
sota <- sota[match(distance.df$Mouse_ID, sota$Mouse_ID),]
all(sota$Mouse_ID==distance.df$Mouse_ID)
rownames(sota) <- sota$Mouse_ID
sota$Mouse_ID <- NULL
gen.dis <- dist.gene(sota, method = "pairwise", pairwise.deletion = TRUE)
gen.dis <- as.matrix(gen.dis)
dimnames(gen.dis) <- c(key, key)

## 3) Sex pairs
Sex_frame<-metadt[,c("Mouse_ID","Sex")]
Sex_frame$Mouse_ID<-as.character(Sex_frame$Mouse_ID)
Sex_frame$Sex<-as.character(Sex_frame$Sex)
#Create an empty character matrix to fill with characters
SEXM<-array(as.character(NA),c(nrow(Sex_frame),nrow(Sex_frame)))

for(i in 1:nrow(Sex_frame)){
    for(j in 1:nrow(Sex_frame)){
        if(Sex_frame$Sex[i]=="F" & Sex_frame$Sex[i]==Sex_frame$Sex[j]){
            SEXM[i,j]= "FF"}
        if(Sex_frame$Sex[i]=="M" & Sex_frame$Sex[i]==Sex_frame$Sex[j]){
           SEXM[i,j]= "MM"}
        if( Sex_frame$Sex[i]!=Sex_frame$Sex[j]){
            SEXM[i,j]= "FM"}
    }
}
dimnames(SEXM)<-c(key, key)

# 4) Making BMI distances
#Create data frame with each sample name (character) and sampling time (numeric)
BMI_frame<-metadt[,c("Mouse_ID", "BMI")]
#Create an empty matrix to fill with distances
BMIM<-array(0,c(nrow(BMI_frame),nrow(BMI_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(BMI_frame)){
    for (j in 1:nrow(BMI_frame))
    {BMIM[i,j]=abs(BMI_frame$BMI[i] -BMI_frame$BMI[j])
    }
}
dimnames(BMIM) <- c(key, key)

# 6) Create farm/Locality matrix: 
#Create data frame with each Individual name (character) and their Age (Character)
Loc_frame<-metadt[,c("Mouse_ID","Locality")]
#Create an empty numeric matrix to fill with distances
LocM<-array(0,c(nrow(Loc_frame),nrow(Loc_frame)))
#Derive matrix with binary locality similarity between each sample
for(i in 1:nrow(Loc_frame)){
    for(j in 1:nrow(Loc_frame)){
        if(Loc_frame$Locality[i]==Loc_frame$Locality[j]){
            LocM[i,j]= "1"
        } else{
            LocM[i,j]= "0"
        }
    }
}
#Note that Locality similarity matrix has rownames and colnames in the same order as key
all(rownames(LocM)==key$ID)
dimnames(LocM) <- c(key, key)


# 6) this matrix will describe the distance in years between samples
#Transform dates into a numeric variable
metadt$Year <- as.numeric(metadt$Year)
#Create data frame with each sample name (character) and sampling time (numeric)
SampleTime_frame<-metadt[,c("Mouse_ID","Year")]
#Create an empty matrix to fill with distances
TEMPM<-array(0,c(nrow(SampleTime_frame),nrow(SampleTime_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(SampleTime_frame)){
 for (j in 1:nrow(SampleTime_frame))
{TEMPM[i,j]=abs(SampleTime_frame$Year[i] -SampleTime_frame$Year[j])
  }
}
dimnames(TEMPM)<-c(key,key)

# 7) Making HI distances
#Create data frame with each sample name (character) and sampling time (numeric)
HI_frame<-metadt[,c("Mouse_ID", "HI")]
#Create an empty matrix to fill with distances
HIM<-array(0,c(nrow(HI_frame),nrow(HI_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(HI_frame)){
    for (j in 1:nrow(HI_frame))
    {HIM[i,j]=abs(HI_frame$HI[i] -HI_frame$HI[j])
    }
}
dimnames(HIM) <- c(key, key)

# 7) Making hi distances
#Create data frame with each sample name (character) and sampling time (numeric)
hi_frame<-metadt[,c("Mouse_ID", "hi")]
#Create an empty matrix to fill with distances
hiM<-array(0,c(nrow(hi_frame),nrow(hi_frame)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(hi_frame)){
    for (j in 1:nrow(hi_frame))
    {hiM[i,j]=abs(hi_frame$hi[i] -hi_frame$hi[j])
    }
}
dimnames(hiM) <- c(key, key)


# here are our matrices
str(BMIM)
str(SPATM)
str(JACM)
str(gen.dis)
str(SEXM)
str(LocM)
str(HIM)
str(hiM)
str(TEMPM)
#First unravel the matrices into vectors matching the lower quantile of each matrix.

#From numeric matrices, this can be done by making a list (c()) of the distance object (dist()) derived from the matrix. as.dist() by default includes only the lower quantile of the matrix and excludes the diagonal.
#From categorical matrices, this can be done by making a list (c()) of the lower quantile of the matrix with lower.tri() -function.

jac<-c(as.dist(JACM))
bmi<-c(as.dist(BMIM))
spa<-c(as.dist(SPATM))
gen<-c(as.dist(gen.dis))
sex<-c(SEXM[lower.tri(SEXM)])
loc <- as.character(c(as.dist(LocM)))
HIm <- c(as.dist(HIM))
him <- c(as.dist(hiM))
tempm <- c(as.dist(TEMPM))
chi <- c(as.dist(CHIM))

#Combine these vectors into a data frame
data.dyad<-data.frame(BMI=bmi,Microbiome_similarity=jac,spatial=spa, genetic_dist=gen, locality=loc, HI=HIm, year=tempm, hi=him, sex=sex, MS=chi)

data.dyad$locality <- as.factor(data.dyad$locality)

#Now all we need to do is add the identities of both individuals in each dyad as separate columns into the data frame and exclude self-comparisons (as these are not meaningful).

# extracting Individual-combinations present in the matrices
list<-expand.grid(key$ID, key$ID)

str(list)

# This created individual-to-same-individual pairs as well. Get rid of these:
list<-list[which(list$Var1!=list$Var2),]

# this still has both quantiles in--> add 'unique' key
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))

# sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices
i=nrow(key)
JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]

# add the names of both individuals participating in each dyad into the data frame
data.dyad$IDA<-list$Var2
data.dyad$IDB<-list$Var1
# Make sure you have got rid of all self comparisons
data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]


######################### Now we model the data ####################

## sex combination into a factor
data.dyad$sex <- factor(data.dyad$sex, levels=c("MM", "FM", "FF"))

#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","genetic_dist", "BMI", "year", "hi")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }

#saveRDS(data.dyad, "tmp/data.dyad.RDS")
#data.dyad <- readRDS("tmp/data.dyad.RDS")

#hist(data.dyad$Microbiome_similarity)
# proportional values semi-normally distributed limited between 0 and 1 not including 1 and 0 --> best use betaregression, but gaussian would probably give similar estimates

#model1_bac<-brm(Microbiome_similarity~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(model1_bac, "tmp/BRMmodel1_bac.rds")

model1_bac_chi<-brm(MS~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(model1_bac_chi, "tmp/BRMmodel1_bac_chi.rds")

model1_bac

remove.packages('ggthemes')

############## Parasites

Parasite <- subset_taxa(PS.TSS, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))

all(key == data.frame(ID=sample_data(Parasite)$Mouse_ID))

####################
## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(Parasite, method="jaccard", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
dimnames(JACM)<- c(key, key)

jac<-c(as.dist(JACM))

## 1) Chisq distance
CHIM <- as.matrix(vegan::vegdist(Parasite@otu_table, method="chisq"))

# transpose Chi square disssimilary matrix to similarty matrix
CHIM <- 1-CHIM
# sanity check
all(rownames(CHIM)==key)
dimnames(CHIM)<- c(key, key)

chi<-c(as.dist(CHIM))
#Combine these vectors into a data frame
data.dyad<-data.frame(BMI=bmi,Microbiome_similarity=jac,spatial=spa, genetic_dist=gen, locality=loc, HI=HIm, year=tempm, hi=him, sex=sex, MS=chi)

data.dyad$locality <- as.factor(data.dyad$locality)

#Now all we need to do is add the identities of both individuals in each dyad as separate columns into the data frame and exclude self-comparisons (as these are not meaningful).

# extracting Individual-combinations present in the matrices
list<-expand.grid(key$ID, key$ID)

# This created individual-to-same-individual pairs as well. Get rid of these:
list<-list[which(list$Var1!=list$Var2),]

# this still has both quantiles in--> add 'unique' key
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))

# sanity check that the Individual name combinations are in the same exact order as the lower quantile value vector of the matrices
i=nrow(key)
JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]

# add the names of both individuals participating in each dyad into the data frame
data.dyad$IDA<-list$Var2
data.dyad$IDB<-list$Var1
# Make sure you have got rid of all self comparisons
data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]


######################### Now we model the data ####################

## sex combination into a factor
data.dyad$sex <- factor(data.dyad$sex, levels=c("MM", "FM", "FF"))

#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","genetic_dist", "BMI", "year", "hi")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }

#saveRDS(data.dyad, "tmp/data.dyad.RDS")
#data.dyad <- readRDS("tmp/data.dyad.RDS")

#hist(data.dyad$Microbiome_similarity)
# proportional values semi-normally distributed limited between 0 and 1 not including 1 and 0 --> best use betaregression, but gaussian would probably give similar estimates

## now let's take a look at correlation trends

#model1_para<-brm(Microbiome_similarity~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(model1_para, "tmp/BRMmodel1_para.rds")

model1_para_chi<-brm(MS~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(model1_para_chi, "tmp/BRMmodel1_para_chi.rds")


####### Diet
Diet <- subset_taxa(PS.TSS, Phylum %in% c("Anthophyta", "Phragmoplastophyta", "Charophyta", "Ochrophyta"))

####################
## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(Diet, method="jaccard", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
dimnames(JACM)<- c(key, key)
jac<-c(as.dist(JACM))
## 1) Chisq distance
CHIM <- as.matrix(vegan::vegdist(Diet@otu_table, method="chisq"))
# transpose Chi square disssimilary matrix to similarty matrix
CHIM <- 1-CHIM
# sanity check
all(rownames(CHIM)==key)
dimnames(CHIM)<- c(key, key)

chim<-c(as.dist(CHIM))

#Combine these vectors into a data frame
data.dyad<-data.frame(BMI=bmi,Microbiome_similarity=jac,spatial=spa, genetic_dist=gen, locality=loc, HI=HIm, year=tempm, hi=him, sex=sex, MS=chim)
data.dyad$locality <- as.factor(data.dyad$locality)
list<-expand.grid(key$ID, key$ID)
list<-list[which(list$Var1!=list$Var2),]
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))
i=nrow(key)
JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]
data.dyad$IDA<-list$Var2
data.dyad$IDB<-list$Var1
data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]
######################### Now we model the data ####################
data.dyad$sex <- factor(data.dyad$sex, levels=c("MM", "FM", "FF"))
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","genetic_dist", "BMI", "year", "hi")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }
#model1_diet<-brm(Microbiome_similarity~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(model1_diet, "tmp/BRMmodel1_diet.rds")


model1_diet_chi<-brm(MS~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(model1_diet_chi, "tmp/BRMmodel1_diet_chi.rds")


####### Fungi
Fungi <- subset_taxa(PS.TSS, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota"))
####################

## 1) Jaccard distance
JACM <- as.matrix(phyloseq::distance(Fungi, method="jaccard", type="samples"))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
dimnames(JACM)<- c(key, key)
jac<-c(as.dist(JACM))
## 1) Chisq distance
CHIM <- as.matrix(vegan::vegdist(Fungi@otu_table, method="chisq"))
# transpose Chi square disssimilary matrix to similarty matrix
CHIM <- 1-CHIM
# sanity check
all(rownames(CHIM)==key)
dimnames(CHIM)<- c(key, key)
chi<-c(as.dist(CHIM))

#Combine these vectors into a data frame
data.dyad<-data.frame(BMI=bmi,Microbiome_similarity=jac,spatial=spa, genetic_dist=gen, locality=loc, HI=HIm, year=tempm, hi=him, sex=sex, MS=chi)

data.dyad$locality <- as.factor(data.dyad$locality)
list<-expand.grid(key$ID, key$ID)
list<-list[which(list$Var1!=list$Var2),]
list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
list<-subset(list, !duplicated(list$key))
i=nrow(key)
JACM[which(rownames(JACM)==list$Var1[i]),which(colnames(JACM)==list$Var2[i])]==jac[i]
data.dyad$IDA<-list$Var2
data.dyad$IDB<-list$Var1
data.dyad<-data.dyad[which(data.dyad$IDA!=data.dyad$IDB),]
######################### Now we model the data ####################
data.dyad$sex <- factor(data.dyad$sex, levels=c("MM", "FM", "FF"))
range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","genetic_dist", "BMI", "year", "hi")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
}

saveRDS(data.dyad, "tmp/Fungi_data.dyad.RDS")

#model1_Fungi<-brm(Microbiome_similarity~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(model1_Fungi, "tmp/BRMmodel1_Fungi.rds")

#model1_Fungi_chi<-brm(MS~1+ spatial+locality+genetic_dist*hi+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(model1_Fungi_chi, "tmp/BRMmodel1_Fungi_chi.rds")

model1_Fungi_chi <- readRDS("tmp/BRMmodel1_Fungi_chi.rds")

model1_Fungi_chi

model1_Fungi

### uploading models
model1 <- readRDS("tmp/BRMmodel1.rds")
model1_para <- readRDS("tmp/BRMmodel1_para.rds")
model1_Fungi <- readRDS("tmp/BRMmodel1_Fungi.rds")
model1_diet <- readRDS("tmp/BRMmodel1_diet.rds")
model1_bac <- readRDS("tmp/BRMmodel1_bac.rds")

para <- summary(model1_para)$fixed

summary(model1)$fixed

rep("test", 10)

para$'l-95% CI'

res.fun <- function(model1_para, name, ASV){
    para <- summary(model1_para)$fixed
    data.frame(Domain=rep(name, 10),
           ASVs=rep(ASV,10),
           Effect=rownames(para),
           Estimate=para$Estimate,
           lCI=para$'l-95% CI',
           uCI=para$'u-95% CI')
}           
           
res <- res.fun(model1_para, "Parasite", 11)
res <- rbind(res, res.fun(model1_bac, "Bacteria", 383))
res <- rbind(res, res.fun(model1_diet, "Diet", 45))
res <- rbind(res, res.fun(model1_Fungi, "Fungi", 65))
res <- rbind(res, res.fun(model1, "Full model", 588))

res$Domain[res$Domain=="Diet"] <- "Plants"

res$Domain[res$Domain=="Parasite"] <- "Parasites"

res$Domain <- factor(res$Domain, levels=c( "Bacteria", "Parasites", "Fungi", "Plants", "Full model"))

library(scales)
#coul=c("#154360", "#b71c1c", "#512e5f", "#0e6251")
#coul=c("#edca82", "#097770", "#e0cdbe", "#a9c0a6")
coul=c("#F8B195","#F67280", "#6C5B7B", "#355C7D")

#Spatial, Locality, genetic, hi

res <- res[res$Effect%in%c("genetic_dist", "hi", "spatial", "locality1"),]

re.plot <- ggplot(res, aes(x = Estimate, y = Effect, fill = Effect)) +
    geom_errorbar(aes(xmin=lCI, xmax=uCI, colour=Effect), size=1, width=0.4)+
        geom_point(shape = 21, size=3) +
    scale_fill_manual(values = coul) +
    scale_colour_manual(values = coul) +
            xlab("Parameter estimate") +
    ylab("") +
    scale_y_discrete(labels = c("Genetic distances", "Hybridicity distances", "Shared locality", "Spatial distances")) +
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    facet_grid(Domain ~ ., scales = "free_y", space = "free_y")+
    theme_classic(base_size=12)+
    theme(
        strip.text = element_text(face = "bold"),
        panel.background = element_rect(fill = "white", color = NA),
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white", color = "black"),
        panel.border = element_rect(color = "black", fill = NA),
    #    axis.line = element_line(color = "black"),
    #    plot.title = element_text(size = 12, face = "bold"),
        legend.position="none")

re.plot
ggsave("fig/figure3.pdf", re.plot, width=150, height=180, units="mm", dpi=300)

############################################################################################
######## Decomposing the Fungi model
# Data wrangling step: Give all unclassified genera a name based on their family or order
tax <- as.data.frame(Fungi@tax_table)
#tax[is.na(tax$Genus),]$Genus
#tax[97,] <- "Unknown_genus_in_Oscillospirales" #  manual change here, need to go back and check
#tax[8,] <- "Unknown_genus_in_Eukarya"

tax_G<-as.data.frame(table(tax[,6]))
genuses<- as.character(tax_G$Var1)

genuses <- c("none", genuses)
#genuses <- append('none', genuses)

dropnumber <- length(genuses)#199
model_minutes<-60
cores <- 60
(dropnumber*model_minutes)/cores

key<-data.frame(ID=sample_data(Fungi)$Mouse_ID, Sample_name=sample_data(Fungi)$Mouse_ID)
data.dyad_REAL <- readRDS("tmp/Fungi_data.dyad.RDS")

names(data.dyad_REAL)

#summary(mcmcglmm_model)

#start cluster
dropRes_list<-list()
cl <- parallel::makeCluster(60, type="FORK")
doParallel::registerDoParallel(cl)

dropRes_list<-foreach(i = 1:length(genuses)) %dopar% {
    library(phyloseq)
    #choose the genus to drop and prune the taxa to keep everything else
    gen.i<-genuses[i]
    taxa_tokeep<-rownames(tax[which(tax$Genus!=genuses[i]),])
    mic.i<-prune_taxa(taxa_tokeep, Fungi)
    #Calculate Jaccard and Bray-Curtis microbiome dissimilarity for each mouse pair
    JACM.i<- as.matrix(phyloseq::distance(mic.i, method="jaccard"))
    BRAY.i<- as.matrix(phyloseq::distance(mic.i, method="bray"))
    #Unravel dissimilarity matrices into vectors
    bray<-c(as.dist(BRAY.i))
    jac<-c(as.dist(JACM.i))

    #Make a new dyadic data frame from these vectors and order it to be in the same order as       the original dyadic data frame
    data.dyad.i<-data.frame(Jaccard=jac,BrayCurtis=bray)
    # extracting Sample_name-combinations of the matrix
    list<-expand.grid(key$Sample_name,key$Sample_name)
    # This created sample-to-same-sample pairs as well. Get rid of these:
    list<-list[which(list$Var1!=list$Var2),]
    # the resulting list still has both quantiles of the original matrix (i.e. all values are doubled) in--> add 'unique' key and subset to one quantile only
    list$key <- apply(list, 1, function(x)paste(sort(x), collapse=''))
    list<-subset(list, !duplicated(list$key))
    # Sample_name combinations are now in the same order as the lower quantile value vector
    # So we can add dyad name and each participant ID to dyadic dataframe
    data.dyad.i$Sample_A<-list$Var2
    data.dyad.i$Sample_B<-list$Var1

    # extracting combinations of individual IDs for each pair
    keyA<-key[,c("ID","Sample_name")]
    colnames(keyA)<-c("IDA","Sample_A")
    keyB<-key[,c("ID","Sample_name")]
    colnames(keyB)<-c("IDB","Sample_B")

    keyA<-keyA[match(data.dyad.i$Sample_A,keyA$Sample_A),]
    keyB<-keyB[match(data.dyad.i$Sample_B,keyB$Sample_B),]

    data.dyad.i$IDA<-keyA$IDA
    data.dyad.i$IDB<-keyB$IDB

    # Make sure you have no self comparisons in the data (This is the case by default here,       since we are using just one sample per individual)
    data.dyad.i<-data.dyad.i[which(data.dyad.i$IDA!=data.dyad.i$IDB),] #

    ### Combine new Jaccard variable with rest of dyadic data columns
    data.dyad<-data.dyad_REAL
    data.dyad$Microbiome_similarity<-data.dyad.i$Jaccard

    #factorize terms used for multimembership random structure and make sure levels are     same and in same order
    data.dyad$IDA<-as.factor(data.dyad$IDA)
    data.dyad$IDB<-as.factor(data.dyad$IDB)
    all(levels(data.dyad$IDA)==levels(data.dyad$IDB))#T

# Scale all predictors not between 0-1 already to be between 0-1
    scalecols<-c("spatial","genetic_dist", "BMI", "year", "hi")
    range.use <- function(x,min.use,max.use){(x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
    for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
        data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }
#Transpose Jaccard dissimilarity to similarity
        data.dyad$Microbiome_similarity<-1-data.dyad$Microbiome_similarity
    
#The MCMCglmm model
dropmodel<-MCMCglmm(Microbiome_similarity~1+spatial+locality+genetic_dist*hi+year+BMI+sex,
                            data=data.dyad,
                            family= "gaussian",
                            random =~ mm(IDA+IDB),
                            verbose=FALSE)

   ASVs_dropped.i<-nrow(tax_table(Fungi))-nrow(tax_table(mic.i))
   resdf.i<-data.frame(Genus_dropped=gen.i,
                      ASVs_dropped=ASVs_dropped.i,
                genetic_dist_Estimate=summary(dropmodel)$solutions["genetic_dist",]["post.mean"],
                     genetic_dist_lCI=summary(dropmodel)$solutions["genetic_dist",]["l-95% CI"],
                     genetic_dist_uCI=summary(dropmodel)$solutions["genetic_dist",]["u-95% CI"],
                     spatial_Estimate=summary(dropmodel)$solutions["spatial",]["post.mean"],
                          spatial_lCI=summary(dropmodel)$solutions["spatial",]["l-95% CI"],
                          spatial_uCI=summary(dropmodel)$solutions["spatial",]["u-95% CI"],
                          hi_Estimate=summary(dropmodel)$solutions["hi",]["post.mean"],
                               hi_lCI=summary(dropmodel)$solutions["hi",]["l-95% CI"],
                               hi_uCI=summary(dropmodel)$solutions["hi",]["u-95% CI"],
                    locality_Estimate=summary(dropmodel)$solutions["locality1",]["post.mean"],
                         locality_lCI=summary(dropmodel)$solutions["locality1",]["l-95% CI"],
                         locality_uCI=summary(dropmodel)$solutions["locality1",]["u-95% CI"],
                        year_Estimate=summary(dropmodel)$solutions["year",]["post.mean"],
                             year_lCI=summary(dropmodel)$solutions["year",]["l-95% CI"],
                             year_uCI=summary(dropmodel)$solutions["year",]["u-95% CI"],
                         BMI_Estimate=summary(dropmodel)$solutions["BMI",]["post.mean"],
                              BMI_lCI=summary(dropmodel)$solutions["BMI",]["l-95% CI"],
                              BMI_uCI=summary(dropmodel)$solutions["BMI",]["u-95% CI"],
                gen_hi_Estimate=summary(dropmodel)$solutions["genetic_dist:hi",]["post.mean"],
                         gen_hi_lCI=summary(dropmodel)$solutions["genetic_dist:hi",]["l-95% CI"],
                gen_hi_uCI=summary(dropmodel)$solutions["genetic_dist:hi",]["u-95% CI"]
                      )
    return(resdf.i)
}
parallel::stopCluster(cl)
saveRDS(dropRes_list,"tmp/Fungi_dropRes_list.rds")

dropRes_list <- readRDS("tmp/Fungi_dropRes_list.rds")
    
#########################################################################
##rbind the resulting data frames to single master data frame
dropResults<-data.frame(Genus_dropped=NA,
                                            ASVs_dropped=NA,
                                            genetic_dist_Estimate=NA,
                                            genetic_dist_lCI=NA,
                                            genetic_dist_uCI=NA,
                                            spatial_Estimate=NA,
                                            spatial_lCI=NA,
                                            spatial_uCI=NA,
                                            hi_Estimate=NA,
                                            hi_lCI=NA,
                                            hi_uCI=NA,
                                            locality_Estimate=NA,
                                            locality_lCI=NA,
                                            locality_uCI=NA,
                                            year_Estimate=NA,
                                            year_lCI=NA,
                                            year_uCI=NA,
                                            BMI_Estimate=NA,
                                            BMI_lCI=NA,
                                            BMI_uCI=NA,
                                            gen_hi_Estimate=NA,
                                            gen_hi_lCI=NA,
                                            gen_hi_uCI=NA)


for(j in 1:length(genuses)){
    dropResults<-rbind(dropResults,dropRes_list[[j]])
  }

dropResults<-dropResults[2:nrow(dropResults),]
dropResults$Genus_dropped2<-genuses

saveRDS(dropResults,"tmp/dropResults_genus.rds")

dropResults <- readRDS("tmp/dropResults_genus.rds")

names(dropResults)

####### genetic dist
#For each microbial genera, calculate their importance on genetic and spatial effect on microbiome
dropResults$genetic_dist_CIbr<-abs(dropResults$genetic_dist_uCI-dropResults$genetic_dist_lCI)
genetic_dist_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$genetic_dist_CIbr
dropResults$genetic_dist_CIbr_increase<- dropResults$genetic_dist_CIbr-genetic_dist_CIbr_baseline


######################## spatial
#For each microbial genera, calculate their importance on social association effect on microbiome
dropResults$spatial_CIbr<-abs(dropResults$spatial_uCI-dropResults$spatial_lCI)
spatial_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$spatial_CIbr
dropResults$spatial_CIbr_increase<- dropResults$spatial_CIbr-spatial_CIbr_baseline

############### locality
######################## spatial
#For each microbial genera, calculate their importance on social association effect on microbiome
dropResults$locality_CIbr<-abs(dropResults$locality_uCI-dropResults$locality_lCI)
locality_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$locality_CIbr
dropResults$locality_CIbr_increase<- dropResults$locality_CIbr-locality_CIbr_baseline

############### hi
######################## spatial
#For each microbial genera, calculate their importance on social association effect on microbiome
dropResults$hi_CIbr<-abs(dropResults$hi_uCI-dropResults$hi_lCI)
hi_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$hi_CIbr
dropResults$hi_CIbr_increase<- dropResults$hi_CIbr-hi_CIbr_baseline

# remove species from tax to merge with dropResults table
tax$Species <- NULL
tax <- unique(tax)

plot_cor <- merge(dropResults, tax, by.x="Genus_dropped", by.y="Genus")

coul=c("#ff6f69", "#ffeead", "#5f8f79") 

plot_cor$genetic_dist_lCI

plot_cor$Genus_dropped

library(dplyr)

plot_cor$Genus_dropped <- factor(plot_cor$Genus_dropped, levels=plot_cor$Genus_dropped[order(plot_cor$Phylum)])

gen <- ggplot(data=plot_cor, aes(x=genetic_dist_Estimate, y=Genus_dropped, fill=Phylum))+
    geom_rect(aes(xmin=dropResults$genetic_dist_lCI[1], xmax=dropResults$genetic_dist_uCI[1], ymin=-Inf, ymax=Inf), fill="#7fb3d5", alpha=0.1)+
    geom_errorbar(aes(xmin = genetic_dist_lCI, xmax = genetic_dist_uCI), size=0.2, alpha=0.5, width=0.2)+
    geom_point(shape=21, size=2, alpha=0.8)+
    scale_fill_manual(values=coul)+
        scale_x_reverse()+
     labs(x="genetic distance", y="Dropped Genus")+ 
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none")
gen

hi <- ggplot(data=plot_cor, aes(x=hi_Estimate, y=Genus_dropped, fill=Phylum))+
    geom_rect(aes(xmin=dropResults$hi_lCI[1], xmax=dropResults$hi_uCI[1], ymin=-Inf, ymax=Inf), fill="#abebc6", alpha=0.1)+
    geom_errorbar(aes(xmin = hi_lCI, xmax = hi_uCI),size=0.2, alpha=0.5)+
    geom_point(shape=21, size=2, alpha=0.8)+
    scale_x_reverse()+
    scale_fill_manual(values=coul)+
    labs(x="hybridicity distance", y="Dropped genus")+  
    geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none",
          axis.title.y=element_blank(),
          axis.text.y=element_blank())
hi


lo <- ggplot(data=plot_cor, aes(x=locality_Estimate, y=Genus_dropped, fill=Phylum))+
    geom_rect(aes(xmin=dropResults$locality_lCI[1], xmax=dropResults$locality_uCI[1], ymin=-Inf, ymax=Inf), fill="#ba4a00", alpha=0.1)+
    geom_errorbar(aes(xmin = locality_lCI, xmax = locality_uCI),size=0.2, alpha=0.5)+
    geom_point(shape=21, size=2, alpha=0.8)+
    scale_fill_manual(values=coul)+
    labs(x="shared locality", y="Dropped genus")+  
    geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none",
          axis.title.y=element_blank(),
          axis.text.y=element_blank())
lo

spa <- ggplot(data=plot_cor, aes(x=spatial_Estimate, y=Genus_dropped, fill=Phylum))+
        geom_rect(aes(xmin=dropResults$spatial_lCI[1], xmax=dropResults$spatial_uCI[1], ymin=-Inf, ymax=Inf), fill="#f1c40f", alpha=0.01)+
    geom_errorbar(aes(xmin = spatial_lCI, xmax = spatial_uCI),size=0.2, alpha=0.5)+
    geom_point(shape=21, size=2, alpha=0.8)+
    scale_fill_manual(values=coul)+
#    scale_x_reverse()+
    labs(x="spatial distance", y="Dropped genus")+
    geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
                theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position="none",
          axis.title.y=element_blank(),
          axis.text.y=element_blank())

spa

fig4 <- plot_grid(gen, hi, spa, lo, nrow=1, rel_widths=c(3,1.5,1.5,1.5))

legend <- get_legend(hi+
                        theme(legend.position="bottom"))

Fig4 <- plot_grid(fig4, legend, rel_heights=c(0.8, 0.02), ncol=1)
Fig4


ggplot2::ggsave(file="fig/Fig4.pdf", Fig4, width = 230, height = 180, dpi = 300, units="mm")
