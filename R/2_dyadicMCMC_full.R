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
library(tidybayes)
library(cowplot)
library(ggpmisc)
library("ggmap")
library(sf, lib.loc="/usr/local/lib/R/site-library")
library("rnaturalearthdata")
library("rnaturalearth")
library(legendMap)
library("legendMap")
library(sf)
library(maps)
library(scico)
library("ggspatial")

PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

#### In this script we are going to model Microbiome beta diversity by geographic distance, shared location, genetic distance, hybridicity distance, sex combination, body mass distance, temporal distance. We do this with dyadic interactions using Bayesian multimodal models.

## I follow Aura Raulo approach
### Questions to asnwer:
#1) genetic distance is a better predictor than hybrid index distance? This has
#implications in specific allele combinations.
#2) spatial distance: is it just locality, or actual physical distance whithin localities?

############# First create dyad data#######################
key <- data.frame(ID=sample_data(PS.TSS)$Mouse_ID)
metadt <- sample_data(PS.TSS)
####################

metadt$He <- 2*(metadt$HI)*(1-metadt$HI)

doData <- TRUE

if(doData){

## 1) Jaccard distance for microbiome
JACM <- as.matrix(phyloseq::distance(PS.TSS, method="jaccard", type="samples", binary=T))
# transpose Jaccard disssimilary matrix to Jaccard similarty matrix
JACM <- 1-JACM
# sanity check
all(rownames(JACM)==key)
jac <- c(as.dist(JACM))

## 2) Spatial distance matrix
distance.df <- metadt[,c("Longitude", "Latitude")]
spa <- c(dist(distance.df, method="euclidean"))

## 3) Temporal (year) distance matrix
#Transform dates into a numeric variable
metadt$Year <- as.numeric(metadt$Year)
tempm <- c(dist(metadt$Year))

## 4) Aitchison distance for microbiome
AIM <- as.matrix(vegan::vegdist(PS.TSS@otu_table, method="aitchison", pseudocount=1))
# transpose Aitchison disssimilary matrix to similarty matrix
AIM <- 1-AIM
ait <- c(as.dist(AIM))

## 5) Making HI distances
HIM <- c(dist(metadt$HI))

## 6) Making HE distances
He <- c(dist(metadt$He))

## 7) Making Hx distances
#Create data frame with each sample name (character) and sampling time (numeric)
hx_frame<-metadt[,c("Mouse_ID", "HI")]
#Create an empty matrix to fill with distances
hxM<-array(0,c(nrow(metadt),nrow(metadt)))
#Derive matrix with time distances between each sample using abs()-function
for (i in 1:nrow(metadt)){
    for (j in 1:nrow(metadt))
    {hxM[i,j]=abs(metadt$He[i] + metadt$He[j])
    }
}
dimnames(hxM) <- c(key, key)
Hx <- c(as.dist(hxM))


#Combine these vectors into a data frame
data.dyad<-data.frame(Microbiome_similarity=jac,spatial=spa, HI=HIM, He=He, Hx=Hx, year=tempm, Microbiome_similarity_ai=ait)

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
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
    range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","He", "year")
for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
    data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }
saveRDS(data.dyad, "tmp/data.dyad.RDS")
}

data.dyad <- readRDS("tmp/data.dyad.RDS")

######################### Now we model the data ####################
## map
# plotting geographic location
capital <-world.cities[c(world.cities$capital==1),]
capital[grep("Republic", capital$country.etc),]
capital <- capital[capital$country.etc%in%"Germany",]
world <- ne_countries(scale="medium", returnclass="sf")
europe <- ne_countries(continent="Europe", returnclass="sf", scale="medium")
metadf <- PS.TSS@sam_data
coordf <- data.frame(Lon=metadf$Longitude, Lat=metadf$Latitude)

## shape of federal states
boundaries <- 
  st_read("tmp/VG250_Bundeslaender_esri.geojson")



sampling <-
    ggplot(data=europe)+
  geom_sf(data = boundaries, fill = NA, color = "black")+
  geom_point(data=metadf,aes(x=Longitude, y=Latitude, fill=HI), size=2, alpha=0.5, shape=21, colour="white")+
    coord_sf(xlim=c(6, 16), ylim=c(47, 56), expand=FALSE)+
 #   geom_point(shape=4, data=capital, aes(x=long, y=lat), size=2, col="black")+
    scale_fill_scico("Hybrid index", palette="roma", direction=-1)+
    scale_bar(lon = 6.5, lat = 55.5, arrow_length = 10, arrow_distance = 50,
                 distance_lon = 50, distance_lat = 7, distance_legend = 20,
                 dist_unit = "km", orientation = FALSE, legend_size = 2)+
    theme_bw(base_size=10)+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))+
  ggspatial::annotation_north_arrow(location = "tr")

hi_HI_dist <-ggplot(data = data.dyad, aes(x= HI, y= He))+
    geom_bin2d(bins=30)+
        scale_fill_scico(palette="bamako", direction=-1)+
#    stat_poly_line(method="lm",formula=y~x+I(x^2),  color="firebrick", size=1)+
#    stat_poly_eq(formula=y~x+I(x^2), aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#               parse = TRUE, size=2) +
    ylab("Expected heterozygousity distance")+
    xlab("Genetic distance")+
    theme_classic(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 10),
          axis.title.y = element_text(vjust = 2, size = 10),
          legend.position="bottom")


hx_HI_dist <-ggplot(data = data.dyad, aes(x= HI, y= Hx))+
    geom_bin2d(bins=30)+
        scale_fill_scico(palette="bamako", direction=-1)+
#    stat_poly_line(method="lm",formula=y~x+I(x^2),  color="firebrick", size=1)+
#    stat_poly_eq(formula=y~x+I(x^2), aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#               parse = TRUE, size=2) +
    ylab("Hybridicity")+
    xlab("Genetic distance")+
    theme_classic(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 10),
          axis.title.y = element_text(vjust = 2, size = 10),
          legend.position="bottom")

hx_He_dist <-ggplot(data = data.dyad, aes(x= He, y= Hx))+
    geom_bin2d(bins=30)+
        scale_fill_scico(palette="bamako", direction=-1)+
#    stat_poly_line(method="lm",formula=y~x+I(x^2),  color="firebrick", size=1)+
#    stat_poly_eq(formula=y~x+I(x^2), aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#               parse = TRUE, size=2) +
    ylab("Hybridicity")+
    xlab("Expected heteroyigousity distance")+
    theme_classic(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 10),
          axis.title.y = element_text(vjust = 2, size = 10),
          legend.position="bottom")

hi_HI <-ggplot(data = metadt, aes(x= HI, y= He))+
    geom_bin2d(bins=40)+
        scale_fill_scico(palette="bamako", direction=-1)+
    ylab("Expected heterozygousity")+
    xlab("Hybrid index")+
    theme_classic(base_size=10)+
    theme(axis.title.x = element_text(vjust = 0, size = 10),
          axis.title.y = element_text(vjust = 2, size = 10),
          legend.position="bottom")

Fig1 <- plot_grid(hi_HI, hi_HI_dist, hx_HI_dist, hx_He_dist, ncol=1)

ggsave("fig/Figure1.pdf", Fig1, width=75, height=350, units="mm", dpi=300)

ggsave("fig/Figure1_map.pdf", sampling, width=100, height=120, units="mm", dpi=300)

## Let's model
modelJ<-brm(Microbiome_similarity~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "zero_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
            inits=0)
modelJ <- add_criterion(modelJ, "loo")
saveRDS(modelJ, "tmp/BRMmodelJac.rds")
modelJ <- readRDS("tmp/BRMmodelJac.rds")
modelA<-brm(Microbiome_similarity_ai~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 6000,
                cores = 20, chains = 4,
                inits=0)
modelA <- add_criterion(modelA, "loo")
saveRDS(modelA, "tmp/BRMmodelA.rds")
modelA <- readRDS("tmp/BRMmodelA.rds")
print(summary(modelA), digits=3) 

############### sanity checks
#1) Is model consistent with mcmcglmm?
#mcmcglmm_modelA<-MCMCglmm(Microbiome_similarity_ai~1+spatial+HI*He+Hx+year,
#                            data=data.dyad,
#                            family= "gaussian",
#                            random =~ mm(IDA+IDB),
#                            verbose=FALSE)
#summary(mcmcglmm_modelA) # yes
#mcmcglmm_modelJ<-MCMCglmm(Microbiome_similarity~1+spatial+HI*He+Hx+year,
#                            data=data.dyad,
#                            family= "gaussian",
#                            random =~ mm(IDA+IDB),
#                            verbose=FALSE)
#summary(mcmcglmm_modelJ) # yes

#2) Denisty overlay = # Compare distribution of response variable to distributions of a set of predicted response variable values based on model -- are they a good fit?
#pp_check(modelA) # fine
#pp_check(modelJ) # fine

# 3) model convergence: catterpillar plots
model1_transformed <- ggs(modelJ)
cat <- filter(model1_transformed, Parameter %in% c("b_Intercept", "b_spatial", "b_HI", "b_Hx", "b_He", "b_year", "b_HI:He"))
par <- c("Intercept", "Spatial distance", "Genetic distance", "hHe-dist", "hHe-mean", "Year distance","genetic distance: hHe-dist")
names(par) <- (unique(model1_transformed$Parameter))[1:7]
ggplot(cat, aes(x=Iteration, y=value,col = as.factor(Chain)))+
    geom_line() +
    geom_vline(xintercept = 1000)+
        scale_color_brewer(palette="Dark2")+
        facet_grid(Parameter ~ . ,
                      scale  = 'free_y',
                   switch = 'y',
                   labeller=as_labeller(par))+
            labs(title = "Caterpillar Plots",
                 col   = "Chains")+
    theme_bw(base_size=10)+
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "right"
    )

model1_transformed_ai <- ggs(modelA)
cat_ai <- filter(model1_transformed_ai, Parameter %in% c("b_Intercept", "b_spatial", "b_HI", "b_Hx", "b_He", "b_year", "b_HI:He"))
ggplot(cat_ai, aes(x=Iteration, y=value, col=as.factor(Chain)))+
    geom_line() +
    geom_vline(xintercept = 1000)+
        scale_color_brewer(palette="Dark2")+
        facet_grid(Parameter ~ . ,
                      scale  = 'free_y',
                   switch = 'y',
                   labeller=as_labeller(par))+
            labs(title = "Caterpillar Plots",
                 col   = "Chains")+
    theme_bw(base_size=10)+
    theme(
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.background = element_blank(),
        axis.line = element_line(colour = "black"),
        legend.position = "none"
          )

########## end of sanity checks

