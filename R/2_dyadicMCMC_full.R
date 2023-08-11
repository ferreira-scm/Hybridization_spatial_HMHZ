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

doData <- FALSE

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
## 3) Chisq distance for microbiome
CHIM <- as.matrix(vegan::vegdist(PS.TSS@otu_table, method="chisq"))
# transpose Chi square disssimilary matrix to similarty matrix
CHIM <- 1-CHIM
# sanity check
all(rownames(CHIM)==key)
chi <- c(as.dist(CHIM))
## 4) Aitchison distance for microbiome
AIM <- as.matrix(vegan::vegdist(PS.TSS@otu_table, method="aitchison", pseudocount=1))
# transpose Aitchison disssimilary matrix to similarty matrix
AIM <- 1-AIM
ait <- c(as.dist(AIM))
## 5) Sex pairs
Sex_frame<-metadt[,c("Mouse_ID","Sex")]
Sex_frame$Mouse_ID<-as.character(Sex_frame$Mouse_ID)
Sex_frame$Sex<-as.character(Sex_frame$Sex)
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
sex<-c(SEXM[lower.tri(SEXM)])
# 6) Making BMI distances
bmi <- c(dist(metadt$BMI))
# 7) Create farm/Locality matrix: 
Loc_frame<-metadt[,c("Mouse_ID","Locality")]
LocM<-array(0,c(nrow(Loc_frame),nrow(Loc_frame)))
for(i in 1:nrow(Loc_frame)){
    for(j in 1:nrow(Loc_frame)){
        if(Loc_frame$Locality[i]==Loc_frame$Locality[j]){
            LocM[i,j]= "1"
        } else{
            LocM[i,j]= "0"
        }
    }
}
all(rownames(LocM)==key$ID)
loc <- as.character(c(as.dist(LocM)))
# 8) this matrix will describe the distance in years between samples
#Transform dates into a numeric variable
metadt$Year <- as.numeric(metadt$Year)
tempm <- c(dist(metadt$Year))
# 9) Making HI distances
HIM <- c(dist(metadt$HI))
#10) Making HE distances
HeM <- c(dist(metadt$He))

#Combine these vectors into a data frame
data.dyad<-data.frame(BMI=bmi,Microbiome_similarity=jac,spatial=spa, locality=loc, HI=HIM, He=HeM, year=tempm, sex=sex, Microbiome_similarity_chi=chi, Microbiome_similarity_ai=ait)
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
## sex combination into a factor
data.dyad$sex <- factor(data.dyad$sex, levels=c("MM", "FM", "FF"))
#scale all predictors to range between 0-1 if they are not already naturally on that scale
#define scaling function:
    range.use <- function(x,min.use,max.use){ (x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
scalecols<-c("spatial","He", "BMI", "year")
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

library("ggspatial")

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
sampling

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
hx_HI_dist

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
hx_He_dist


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
ggsave("fig/Figure1.pdf", Fig1, width=85, height=180, units="mm", dpi=300)

ggsave("fig/Figure1_map.pdf", sampling, width=100, height=120, units="mm", dpi=300)

## Let's model
#model1<-brm(Microbiome_similarity~1+ spatial+locality+HI*He+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#model1 <- add_criterion(model1, "loo")
#saveRDS(model1, "tmp/BRMmodel1.rds")
model1 <- readRDS("tmp/BRMmodel1.rds")
model1S<-brm(Microbiome_similarity~1+ spatial+locality+HI*He+Hx+year+BMI+sex+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
model1S <- add_criterion(model1S, "loo")
saveRDS(model1S, "tmp/BRMmodel1S.rds")
model1S <- readRDS("tmp/BRMmodel1S.rds")

#model1S2<-brm(Microbiome_similarity~1+ spatial+locality+HI+Hx+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#model1S2 <- add_criterion(model1S2, "loo")
#saveRDS(model1S2, "tmp/BRMmodel1S2.rds")
model1S2 <- readRDS("tmp/BRMmodel1S2.rds")

#model1_ai<-brm(Microbiome_similarity_ai~1+ spatial+locality+HI*He+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#model1_ai <- add_criterion(model1_ai, "loo")
#saveRDS(model1_ai, "tmp/BRMmodel1_ai.rds")
model1_ai <- readRDS("tmp/BRMmodel1_ai.rds")

model1_aiS<-brm(Microbiome_similarity_ai~1+ spatial+locality+HI*He+Hx+year+BMI+sex+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
                inits=0)
model1_aiS <- add_criterion(model1_aiS, "loo")
saveRDS(model1_aiS, "tmp/BRMmodel1_aiS.rds")
model1_aiS <- readRDS("tmp/BRMmodel1_aiS.rds")

#model1_aiS2<-brm(Microbiome_similarity_ai~1+ spatial+locality+HI+Hx+year+BMI+sex+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
#               inits=0)
#model1_aiS2 <- add_criterion(model1_aiS2, "loo")
#saveRDS(model1_aiS2, "tmp/BRMmodel1_aiS2.rds")
model1_aiS2 <- readRDS("tmp/BRMmodel1_aiS2.rds")


loo(model1)


#mcmcglmm_model<-MCMCglmm(Microbiome_similarity~1+spatial+locality+HI+He+year+BMI+sex,
#                            data=data.dyad,
#                            family= "gaussian",
#                            random =~ mm(IDA+IDB),
#                            verbose=FALSE)
mcmcglmm_model <- readRDS("tmp/mcmcglmm_model.rds")
summary(mcmcglmm_model)

#mcmcglmm_model_ai<-MCMCglmm(Microbiome_similarity_chi~1+spatial+locality+HI+He+year+BMI+sex,
#                            data=data.dyad,
#                            family= "gaussian",
#                            random =~ mm(IDA+IDB),
#                            verbose=FALSE)
#saveRDS(mcmcglmm_model, "tmp/mcmcglmm_model.rds")
mcmcglmm_model_ai <- readRDS("tmp/mcmcglmm_model.rds")
summary(mcmcglmm_model_chi)

#Denisty overlay = # Compare distribution of response variable to distributions of a set of predicted response variable values based on model -- are they a good fit?
#pp_check(model1) # fine

# model convergence for jaccard
model1_transformed <- ggs(model1)

cat <- filter(model1_transformed, Parameter %in% c("b_Intercept", "b_spatial", "b_locality1", "b_HI", "b_He", "b_year", "b_BMI", "b_sexFM", "b_sexFF", "b_HI:He"))
par <- c("Intercept", "Spatial distance", "Shared locality", "HI", "He", "Year distance",
         "BMI distance", "Female-Male", "Female-female", "HI*He")
names(par) <- (unique(model1_transformed$Parameter))[1:10]
caterpillar <- ggplot(filter(model1_transformed, Parameter %in% c("b_Intercept", "b_spatial", "b_locality1", "b_HI", "b_He", "b_year", "b_BMI", "b_sexFM", "b_sexFF", "b_HI:He")),
       aes(x   = Iteration,
           y   = value,
           col = as.factor(Chain)))+
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

# model convergence for chi
model1_transformed_ai <- ggs(model1_ai)
cat_ai <- filter(model1_transformed_ai, Parameter %in% c("b_Intercept", "b_spatial", "b_locality1", "b_HI", "b_He", "b_year", "b_BMI", "b_sexFM", "b_sexFF", "b_HI:He"))
par <- c("Intercept", "Spatial distance", "Shared locality", "HI", "He", "Year distance",
         "BMI distance", "Female-Male", "Female-female", "HI*He")
names(par) <- (unique(model1_transformed_ai$Parameter))[1:10]
caterpillar_ai <- ggplot(filter(model1_transformed_ai, Parameter %in% c("b_Intercept", "b_spatial", "b_locality1", "b_HI", "b_He", "b_year", "b_BMI", "b_sexFM", "b_sexFF", "b_HI:He")),
       aes(x   = Iteration,
           y   = value,
           col = as.factor(Chain)))+
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

caterp <- plot_grid(caterpillar, caterpillar_ai, ncol=2, labels="auto")

ggsave("fig/figureS1_caterpillar.pdf", caterpillar, width=170, height=250, units="mm", dpi=300)

## don't forget to subset here
cat2 <- cat[cat$Iteration>1000,]
cat2$Parameter <- droplevels(cat2$Parameter)
cat2$Parameter <- factor(cat2$Parameter, levels=c("b_Intercept", "b_spatial", "b_locality1", "b_HI", "b_He", "b_year", "b_BMI", "b_sexFM",  "b_sexFF", "b_HI:He"))

cat2_ai <- cat_ai[cat_ai$Iteration>1000,]
cat2_ai$Parameter <- droplevels(cat2_ai$Parameter)
cat2_ai$Parameter <- factor(cat2_ai$Parameter, levels=c("b_Intercept", "b_spatial", "b_locality1", "b_HI", "b_He", "b_year", "b_BMI", "b_sexFM",  "b_sexFF", "b_HI:He"))


newname <- c("Intercept", "Spatial distance", "Shared locality", "HI distance", "He distance", "Temporal distance", "Body mass index distance", "Female-Male", "Female-Female", "HI:He")# rename"Spatial distance", "Shared locality", "Genetic distance", "Hybridicity distance", "Temporal distance", "Body mass index distance", "Female-Male", "Female-Female", "Genetic*hybridicity distance")
name <- unique(cat2$Parameter)
for (i in 1:10){
    cat2$Parameter <- gsub(name[i], newname[i], cat2$Parameter)
}
for (i in 1:10){
    cat2_ai$Parameter <- gsub(name[i], newname[i], cat2_ai$Parameter)
}

cat2$Parameter[cat2$Parameter=="HI distance:He"] <- "HI distance: He distance"
cat2 <- cat2[!cat2$Parameter%in%"Intercept",]
cat2_ai$Parameter[cat2_ai$Parameter=="HI distance:He"] <- "HI distance: He distance"
cat2_ai <- cat2_ai[!cat2_ai$Parameter%in%"Intercept",]

cat2$Parameter <- factor(cat2$Parameter, levels=c("Temporal distance", "Body mass index distance", "Female-Male", "Female-Female", "HI distance: He distance", "Spatial distance", "Shared locality", "HI distance", "He distance"))

cat2_ai$Parameter <- factor(cat2_ai$Parameter, levels=c("Temporal distance", "Body mass index distance", "Female-Male", "Female-Female", "HI distance: He distance", "Spatial distance", "Shared locality", "HI distance", "He distance"))

colours <- c("#979a9a", "#979a9a","#d0d3d4","#d0d3d4","#d0d3d4","#355C7D", "#6C5B7B","#F8B195", "#F67280")

FigX <-    ggplot(cat2, aes(x = value, y=Parameter, fill=Parameter))+
    geom_density_ridges(rel_min_height = 0.005, scale=5, alpha=0.8)+
    geom_vline(xintercept = 0, col  = "black", size = 1, linetype="dashed")+
    geom_boxplot(outlier.shape = NA,
                 width=0.2)+
    scale_fill_manual(values=colours)+
    xlab("Posterior probability distribution")+
    ylab("")+
    theme_bw(base_size=12)+
    theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 12, face = "bold")
          )

FigX_ai <-    ggplot(cat2_ai, aes(x = value, y=Parameter, fill=Parameter))+
    geom_density_ridges(rel_min_height = 0.005, scale=5, alpha=0.8)+
    geom_vline(xintercept = 0, col  = "black", size = 1, linetype="dashed")+
    geom_boxplot(outlier.shape = NA,
                 width=0.2)+
    scale_fill_manual(values=colours)+
    xlab("Posterior probability distribution")+
    ylab("")+
    theme_bw(base_size=12)+
    theme(
        legend.position = "none",
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        axis.title.y = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(size = 9),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size = 12, face = "bold")
          )


FigureS2 <- plot_grid(FigX, FigX_ai, labels="auto", ncol=1)

ggsave("fig/figureS2.pdf", FigureS2, width=180, height=180, units="mm", dpi=300)


#Now we decompose the full models with jaccard and then chi-squared distances, to see how the results are dependent on particular genus
###################
# sanity check
all(sample_names(PS.TSS)==key$ID)

# Data wrangling step: Give all unclassified genera a name based on their family or order
tax <- as.data.frame(PS.TSS@tax_table)
#tax[is.na(tax$Genus),]$Genus
#tax[97,] <- "Unknown_genus_in_Oscillospirales" #  manual change here, need to go back and check
#tax[8,] <- "Unknown_genus_in_Eukarya"

tax_G<-as.data.frame(table(tax[,6]))
genuses<- as.character(tax_G$Var1)
genuses <- c("none", genuses)
#genuses <- append('none', genuses)

key<-data.frame(ID=sample_data(PS.TSS)$Mouse_ID, Sample_name=sample_data(PS.TSS)$Mouse_ID)
data.dyad_REAL <- readRDS("tmp/data.dyad.RDS")

doLeaveOneOut <- FALSE

if(doLeaveOneOut){
#start cluster
dropRes_list_J<-list()
cl <- parallel::makeCluster(40, type="FORK")
doParallel::registerDoParallel(cl)
dropRes_list_J<-foreach(i = 1:length(genuses)) %dopar% {
    library(phyloseq)
    #choose the genus to drop and prune the taxa to keep everything else
    gen.i<-genuses[i]
    taxa_tokeep<-rownames(tax[which(tax$Genus!=genuses[i]),])
    mic.i<-prune_taxa(taxa_tokeep, PS.TSS)
    #Calculate Jaccard and Aitchison dissimilarity for each mouse pair
    JACM.i<- as.matrix(phyloseq::distance(mic.i, method="jaccard", type="samples", binary=T))
    AIT.i <- as.matrix(vegan::vegdist(mic.i@otu_table, method="aitchison", pseudocount=1))
    #Unravel dissimilarity matrices into vectors
    jac<-c(as.dist(JACM.i))
    ait<-c(as.dist(AIT.i))
    #Make a new dyadic data frame from these vectors and order it to be in the same order as       the original dyadic data frame
    data.dyad.i<-data.frame(Jaccard=jac,Aitchison=ait)
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
    ### Combine new variables with rest of dyadic data columns
    data.dyad<-data.dyad_REAL
    data.dyad$MS_J<-data.dyad.i$Jaccard
    data.dyad$MS_A<-data.dyad.i$Aitchison
    #factorize terms used for multimembership random structure and make sure levels are     same and in same order
    data.dyad$IDA<-as.factor(data.dyad$IDA)
    data.dyad$IDB<-as.factor(data.dyad$IDB)
    all(levels(data.dyad$IDA)==levels(data.dyad$IDB))#T
# Scale all predictors not between 0-1 already to be between 0-1
    scalecols<-c("spatial", "BMI", "year", "He")
    range.use <- function(x,min.use,max.use){(x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
    for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
        data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }
#Transpose dissimilarity to similarity
    data.dyad$MS_J<-1-data.dyad$MS_J
    data.dyad$MS_A<-1-data.dyad$MS_A
#The MCMCglmm model
dropmodel_J<-MCMCglmm(MS_J~1+spatial+locality+HI*He+year+BMI+sex,
                            data=data.dyad,
                            family= "gaussian",
                            random =~ mm(IDA+IDB),
                            verbose=FALSE)
   ASVs_dropped.i<-nrow(tax_table(PS.TSS))-nrow(tax_table(mic.i))
   resdfJ.i<-data.frame(Genus_dropped=gen.i,
                      ASVs_dropped=ASVs_dropped.i,
                HI_Estimate=summary(dropmodel_J)$solutions["HI",]["post.mean"],
                     HI_lCI=summary(dropmodel_J)$solutions["HI",]["l-95% CI"],
                     HI_uCI=summary(dropmodel_J)$solutions["HI",]["u-95% CI"],
                     spatial_Estimate=summary(dropmodel_J)$solutions["spatial",]["post.mean"],
                          spatial_lCI=summary(dropmodel_J)$solutions["spatial",]["l-95% CI"],
                          spatial_uCI=summary(dropmodel_J)$solutions["spatial",]["u-95% CI"],
                          He_Estimate=summary(dropmodel_J)$solutions["He",]["post.mean"],
                               He_lCI=summary(dropmodel_J)$solutions["He",]["l-95% CI"],
                               He_uCI=summary(dropmodel_J)$solutions["He",]["u-95% CI"],
                    locality_Estimate=summary(dropmodel_J)$solutions["locality1",]["post.mean"],
                         locality_lCI=summary(dropmodel_J)$solutions["locality1",]["l-95% CI"],
                         locality_uCI=summary(dropmodel_J)$solutions["locality1",]["u-95% CI"],
                        year_Estimate=summary(dropmodel_J)$solutions["year",]["post.mean"],
                             year_lCI=summary(dropmodel_J)$solutions["year",]["l-95% CI"],
                             year_uCI=summary(dropmodel_J)$solutions["year",]["u-95% CI"],
                         BMI_Estimate=summary(dropmodel_J)$solutions["BMI",]["post.mean"],
                              BMI_lCI=summary(dropmodel_J)$solutions["BMI",]["l-95% CI"],
                              BMI_uCI=summary(dropmodel_J)$solutions["BMI",]["u-95% CI"],
                       HI_He_Estimate=summary(dropmodel_J)$solutions["HI:He",]["post.mean"],
                           HI_He_lCI=summary(dropmodel_J)$solutions["HI:He",]["l-95% CI"],
                            HI_He_uCI=summary(dropmodel_J)$solutions["HI:He",]["u-95% CI"]
                      )
    return(resdfJ.i)
}
parallel::stopCluster(cl)
    saveRDS(dropRes_list_J,"tmp/dropRes_list_J.rds")
}

dropRes_list <- readRDS("tmp/dropRes_list_J.rds")
#########################################################################
##rbind the resulting data frames to single master data frame
dropResults<-data.frame(Genus_dropped=NA,
                                            ASVs_dropped=NA,
                                            HI_Estimate=NA,
                                            HI_lCI=NA,
                                            HI_uCI=NA,
                                            spatial_Estimate=NA,
                                            spatial_lCI=NA,
                                            spatial_uCI=NA,
                                            He_Estimate=NA,
                                            He_lCI=NA,
                                            He_uCI=NA,
                                            locality_Estimate=NA,
                                            locality_lCI=NA,
                                            locality_uCI=NA,
                                            year_Estimate=NA,
                                            year_lCI=NA,
                                            year_uCI=NA,
                                            BMI_Estimate=NA,
                                            BMI_lCI=NA,
                                            BMI_uCI=NA,
                                            HI_He_Estimate=NA,
                                            HI_He_lCI=NA,
                                            HI_He_uCI=NA)


for(j in 1:length(genuses)){
    dropResults<-rbind(dropResults,dropRes_list[[j]])
  }
dropResults<-dropResults[2:nrow(dropResults),]
dropResults$Genus_dropped2<-genuses

####### genetic dist
#For each microbial genera, calculate their importance on genetic and spatial effect on microbiome
dropResults$HI_CIbr<-abs(dropResults$HI_uCI-dropResults$HI_lCI)
HI_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$HI_CIbr
dropResults$HI_CIbr_increase<- dropResults$HI_CIbr-HI_CIbr_baseline


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
dropResults$He_CIbr<-abs(dropResults$He_uCI-dropResults$He_lCI)
He_CIbr_baseline<-dropResults[which(dropResults$Genus_dropped=="none"),]$He_CIbr
dropResults$He_CIbr_increase<- dropResults$He_CIbr-He_CIbr_baseline

tax <- tax[, c("Genus", "Phylum", "Kingdom")]
tax <- unique(tax)
plot_cor <- merge(dropResults, tax, by.x="Genus_dropped", by.y="Genus")
plot_cor$Group <- "Other"
plot_cor$Group[plot_cor$Kingdom=="Bacteria"] <- "Bacteria"
plot_cor$Group[plot_cor$Phylum %in% c("Anthophyta", "Phragmoplastophyta", "Charophyta", "Ochrophyta")] <- "Plant"
plot_cor$Group[plot_cor$Genus_dropped %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas")] <- "Parasite"
plot_cor$Group[plot_cor$Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota")] <- "Fungi"
plot_cor$Group <- factor(plot_cor$Group, levels=c("Bacteria", "Parasite", "Fungi", "Plant", "Other"))
#quick fix
plot_cor$Phylum <- NULL
plot_cor <- (unique(plot_cor))

coul=c("#ffe599", "#a13030", "#d8a91c", "#29431d", "#783f04")

# cheching the interesting taxa
int_loc1 <- plot_cor[plot_cor$locality_uCI<dropResults$locality_lCI[dropResults$Genus_dropped=="none"],]
annotation <- data.frame(x=int_loc1$locality_Estimate, Group=int_loc1$Group, label=int_loc1$Genus_dropped)
annotation$x <- annotation$x-0.01
int_loc2 <- plot_cor[plot_cor$locality_lCI>dropResults$locality_uCI[dropResults$Genus_dropped=="none"],]
int_loc2$locality_Estimate <- int_loc2$locality_Estimate+0.01
annotation <- rbind(annotation, data.frame(x=int_loc2$locality_Estimate, Group=int_loc2$Group, label=int_loc2$Genus_dropped))

loc <- ggplot(data=plot_cor, aes(x=locality_Estimate, y=Group, fill=Group))+
    geom_rect(aes(xmin=dropResults$locality_lCI[dropResults$Genus_dropped=="none"], xmax=dropResults$locality_uCI[dropResults$Genus_dropped=="none"], ymin=-Inf, ymax=Inf),
              fill="#6C5B7B", alpha=0.05)+
    geom_errorbar(aes(xmin = locality_lCI, xmax = locality_uCI),size=0.7, alpha=0.3, width=0.1,  position=position_jitter(seed=40), colour="black")+
    geom_point(shape=21, size=2, alpha=0.8, position=position_jitter(seed=40))+
    scale_fill_manual(values=coul)+
    labs(x="Posterior estimate of shared locality", y="")+
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    geom_hline(yintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position= "none")+
    geom_text(data=annotation, aes(x=x, y=Group, label=label),
              colour="black", position=position_jitter(height=0.2, seed=80),
             size=2)

# cheching the interesting taxa
int <- plot_cor[plot_cor$spatial_uCI<dropResults$spatial_lCI[dropResults$Genus_dropped=="none"],]
annotation <- data.frame(x=int$spatial_Estimate, Group=int$Group, label=int$Genus_dropped)
#annotation$x <- annotation$x-0.01
int2 <- plot_cor[plot_cor$spatial_lCI>dropResults$spatial_uCI[dropResults$Genus_dropped=="none"],]
#int_loc2$spatial_Estimate <- int_loc2$spatial_Estimate+0.01
annotation <- rbind(annotation, data.frame(x=int2$spatial_Estimate, Group=int2$Group, label=int2$Genus_dropped))

spa <- ggplot(data=plot_cor, aes(x=spatial_Estimate, y=Group, fill=Group))+
    geom_rect(aes(xmin=dropResults$spatial_lCI[dropResults$Genus_dropped=="none"], xmax=dropResults$spatial_uCI[dropResults$Genus_dropped=="none"], ymin=-Inf, ymax=Inf),
              fill="#355C7D", alpha=0.05)+
    geom_errorbar(aes(xmin = spatial_lCI, xmax = spatial_uCI),size=0.7, alpha=0.3, width=0.1,  position=position_jitter(seed=40), colour="black")+
    geom_point(shape=21, size=2, alpha=0.8, position=position_jitter(seed=40))+
        scale_x_reverse()+
    scale_fill_manual(values=coul)+
    labs(x="Posterior estimate of spatial distance", y="")+
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    geom_hline(yintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position= "none")+
        geom_text(data=annotation, aes(x=x, y=Group, label=label),
              colour="black", position=position_jitter(height=0.2, seed=11),
             size=2)

int <- plot_cor[plot_cor$hi_uCI<dropResults$hi_lCI[dropResults$Genus_dropped=="none"],]

annotation <- data.frame(x=int$hi_Estimate, Group=int$Group, label=int$Genus_dropped)
#annotation$x <- annotation$x-0.01
int2 <- plot_cor[plot_cor$hi_lCI>dropResults$hi_uCI[dropResults$Genus_dropped=="none"],]
#int_loc2$spatial_Estimate <- int_loc2$spatial_Estimate+0.01
annotation <- rbind(annotation, data.frame(x=int2$hi_Estimate, Group=int2$Group, label=int2$Genus_dropped))

int <-plot_cor[plot_cor$He_lCI>0,]
annotation <- data.frame(x=int$He_Estimate, Group=int$Group, label=int$Genus_dropped)

He <- ggplot(data=plot_cor, aes(x=He_Estimate, y=Group, fill=Group))+
    geom_rect(aes(xmin=dropResults$He_lCI[dropResults$Genus_dropped=="none"], xmax=dropResults$He_uCI[dropResults$Genus_dropped=="none"], ymin=-Inf, ymax=Inf),
              fill="#f897a1", alpha=0.05)+
    geom_errorbar(aes(xmin = He_lCI, xmax = He_uCI),size=0.7, alpha=0.3, width=0.1,  position=position_jitter(seed=40), colour="black")+
    geom_point(shape=21, size=2, alpha=0.8, position=position_jitter(seed=40))+
    scale_fill_manual(values=coul)+
    labs(x="Posterior estimate of He distance", y="")+
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    geom_hline(yintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position= "none")+
            geom_text(data=annotation, aes(x=x+0.005, y=Group, label=label),
              colour="black", position=position_jitter(height=0.2, seed=10),
             size=2)


int <-plot_cor[plot_cor$HI_uCI<dropResults$HI_lCI[dropResults$Genus_dropped=="none"],]
annotation <- data.frame(x=int$HI_Estimate, Group=int$Group, label=int$Genus_dropped)
#annotation$x <- annotation$x-0.01
int2 <- plot_cor[plot_cor$HI_lCI>dropResults$HI_uCI[dropResults$Genus_dropped=="none"],]
#int_loc2$spatial_Estimate <- int_loc2$spatial_Estimate+0.01
annotation <- rbind(annotation, data.frame(x=int2$HI_Estimate, Group=int2$Group, label=int2$Genus_dropped))
HI <- ggplot(data=plot_cor, aes(x=HI_Estimate, y=Group, fill=Group))+
    geom_rect(aes(xmin=dropResults$HI_lCI[dropResults$Genus_dropped=="none"], xmax=dropResults$HI_uCI[dropResults$Genus_dropped=="none"], ymin=-Inf, ymax=Inf),
              fill="#F8B195", alpha=0.05)+
    geom_errorbar(aes(xmin = HI_lCI, xmax = HI_uCI),size=0.7, alpha=0.3, width=0.1,  position=position_jitter(seed=40), colour="black")+
    geom_point(shape=21, size=2, alpha=0.8, position=position_jitter(seed=40))+
    scale_fill_manual(values=coul)+
        scale_x_reverse()+
    labs(x="Posterior estimate of HI distances", y="")+
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    geom_hline(yintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position= "none")+
                geom_text(data=annotation, aes(x=x+0.005, y=Group, label=label),
              colour="black", position=position_jitter(height=0.2, seed=10),
             size=2)


Fig3 <- plot_grid(spa, loc, He, HI, labels="auto")

##############################################################################
############################## now for ait ####################################
###############################################################################

if(doLeaveOneOut){
#start cluster
dropRes_list_A<-list()
cl <- parallel::makeCluster(60, type="FORK")
doParallel::registerDoParallel(cl)
dropRes_list_A<-foreach(i = 1:length(genuses)) %dopar% {
    library(phyloseq)
    #choose the genus to drop and prune the taxa to keep everything else
    gen.i<-genuses[i]
    taxa_tokeep<-rownames(tax[which(tax$Genus!=genuses[i]),])
    mic.i<-prune_taxa(taxa_tokeep, PS.TSS)
    #Calculate Jaccard and Aitchison dissimilarity for each mouse pair
    JACM.i<- as.matrix(phyloseq::distance(mic.i, method="jaccard", type="samples", binary=T))
    AIT.i <- as.matrix(vegan::vegdist(mic.i@otu_table, method="aitchison", pseudocount=1))
    #Unravel dissimilarity matrices into vectors
    jac<-c(as.dist(JACM.i))
    ait<-c(as.dist(AIT.i))
    #Make a new dyadic data frame from these vectors and order it to be in the same order as       the original dyadic data frame
    data.dyad.i<-data.frame(Jaccard=jac,Aitchison=ait)
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
    ### Combine new variables with rest of dyadic data columns
    data.dyad<-data.dyad_REAL
    data.dyad$MS_J<-data.dyad.i$Jaccard
    data.dyad$MS_A<-data.dyad.i$Aitchison
    #factorize terms used for multimembership random structure and make sure levels are     same and in same order
    data.dyad$IDA<-as.factor(data.dyad$IDA)
    data.dyad$IDB<-as.factor(data.dyad$IDB)
    all(levels(data.dyad$IDA)==levels(data.dyad$IDB))#T
# Scale all predictors not between 0-1 already to be between 0-1
    scalecols<-c("spatial", "BMI", "year", "He")
    range.use <- function(x,min.use,max.use){(x - min(x,na.rm=T))/(max(x,na.rm=T)-min(x,na.rm=T)) * (max.use - min.use) + min.use }
    for(i in 1:ncol(data.dyad[,which(colnames(data.dyad)%in%scalecols)])){
        data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i]<-range.use(data.dyad[,which(colnames(data.dyad)%in%scalecols)][,i],0,1)
    }
#Transpose dissimilarity to similarity
    data.dyad$MS_J<-1-data.dyad$MS_J
    data.dyad$MS_A<-1-data.dyad$MS_A
#The MCMCglmm model
dropmodel_A<-MCMCglmm(MS_A~1+spatial+locality+HI*He+year+BMI+sex,
                            data=data.dyad,
                            family= "gaussian",
                            random =~ mm(IDA+IDB),
                            verbose=FALSE)
   ASVs_dropped.i<-nrow(tax_table(PS.TSS))-nrow(tax_table(mic.i))
   resdfA.i<-data.frame(Genus_dropped=gen.i,
                      ASVs_dropped=ASVs_dropped.i,
                HI_Estimate=summary(dropmodel_A)$solutions["HI",]["post.mean"],
                     HI_lCI=summary(dropmodel_A)$solutions["HI",]["l-95% CI"],
                     HI_uCI=summary(dropmodel_A)$solutions["HI",]["u-95% CI"],
                     spatial_Estimate=summary(dropmodel_A)$solutions["spatial",]["post.mean"],
                          spatial_lCI=summary(dropmodel_A)$solutions["spatial",]["l-95% CI"],
                          spatial_uCI=summary(dropmodel_A)$solutions["spatial",]["u-95% CI"],
                          He_Estimate=summary(dropmodel_A)$solutions["He",]["post.mean"],
                               He_lCI=summary(dropmodel_A)$solutions["He",]["l-95% CI"],
                               He_uCI=summary(dropmodel_A)$solutions["He",]["u-95% CI"],
                    locality_Estimate=summary(dropmodel_A)$solutions["locality1",]["post.mean"],
                         locality_lCI=summary(dropmodel_A)$solutions["locality1",]["l-95% CI"],
                         locality_uCI=summary(dropmodel_A)$solutions["locality1",]["u-95% CI"],
                        year_Estimate=summary(dropmodel_A)$solutions["year",]["post.mean"],
                             year_lCI=summary(dropmodel_A)$solutions["year",]["l-95% CI"],
                             year_uCI=summary(dropmodel_A)$solutions["year",]["u-95% CI"],
                         BMI_Estimate=summary(dropmodel_A)$solutions["BMI",]["post.mean"],
                              BMI_lCI=summary(dropmodel_A)$solutions["BMI",]["l-95% CI"],
                              BMI_uCI=summary(dropmodel_A)$solutions["BMI",]["u-95% CI"],
                       HI_He_Estimate=summary(dropmodel_A)$solutions["HI:He",]["post.mean"],
                            HI_He_lCI=summary(dropmodel_A)$solutions["HI:He",]["l-95% CI"],
                            HI_He_uCI=summary(dropmodel_A)$solutions["HI:He",]["u-95% CI"]
                      )
    return(resdfA.i)
}
parallel::stopCluster(cl)
    saveRDS(dropRes_list_A,"tmp/dropRes_list_A.rds")
}

dropRes_list_A <- readRDS("tmp/dropRes_list_A.rds")

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

tax <- tax[, c("Genus", "Phylum", "Kingdom")]
tax <- unique(tax)
plot_cor <- merge(dropResults, tax, by.x="Genus_dropped", by.y="Genus")
plot_cor$Group <- "Other"
plot_cor$Group[plot_cor$Kingdom=="Bacteria"] <- "Bacteria"
plot_cor$Group[plot_cor$Phylum %in% c("Anthophyta", "Phragmoplastophyta", "Charophyta", "Ochrophyta")] <- "Plant"
plot_cor$Group[plot_cor$Genus_dropped %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas")] <- "Parasite"
plot_cor$Group[plot_cor$Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota")] <- "Fungi"
plot_cor$Group <- factor(plot_cor$Group, levels=c("Bacteria", "Parasite", "Fungi", "Plant", "Other"))
#quick fix
plot_cor$Phylum <- NULL
plot_cor <- (unique(plot_cor))

coul=c("#ffe599", "#a13030", "#d8a91c", "#29431d", "#783f04")

# cheching the interesting taxa
plot_cor[plot_cor$locality_uCI<dropResults$locality_lCI[dropResults$Genus_dropped=="none"],]
loc <- ggplot(data=plot_cor, aes(x=locality_Estimate, y=Group, fill=Group))+
    geom_rect(aes(xmin=dropResults$locality_lCI[dropResults$Genus_dropped=="none"], xmax=dropResults$locality_uCI[dropResults$Genus_dropped=="none"], ymin=-Inf, ymax=Inf),
              fill="#6C5B7B", alpha=0.05)+
    geom_errorbar(aes(xmin = locality_lCI, xmax = locality_uCI),size=0.7, alpha=0.3, width=0.1,  position=position_jitter(seed=40), colour="black")+
    geom_point(shape=21, size=2, alpha=0.8, position=position_jitter(seed=40))+
    scale_fill_manual(values=coul)+
    labs(x="Posterior estimate of shared locality", y="")+
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    geom_hline(yintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position= "none")

# cheching the interesting taxa
plot_cor[plot_cor$spatial_uCI<dropResults$spatial_lCI[dropResults$Genus_dropped=="none"],]
spa <- ggplot(data=plot_cor, aes(x=spatial_Estimate, y=Group, fill=Group))+
    geom_rect(aes(xmin=dropResults$spatial_lCI[dropResults$Genus_dropped=="none"], xmax=dropResults$spatial_uCI[dropResults$Genus_dropped=="none"], ymin=-Inf, ymax=Inf),
              fill="#355C7D", alpha=0.05)+
    geom_errorbar(aes(xmin = spatial_lCI, xmax = spatial_uCI),size=0.7, alpha=0.3, width=0.1,  position=position_jitter(seed=40), colour="black")+
    geom_point(shape=21, size=2, alpha=0.8, position=position_jitter(seed=40))+
        scale_x_reverse()+
    scale_fill_manual(values=coul)+
    labs(x="Posterior estimate of spatial distance", y="")+
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    geom_hline(yintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position= "none")


plot_cor[plot_cor$hi_uCI<dropResults$hi_lCI[dropResults$Genus_dropped=="none"],]
plot_cor[plot_cor$hi_lCI>0,]

hi <- ggplot(data=plot_cor, aes(x=hi_Estimate, y=Group, fill=Group))+
    geom_rect(aes(xmin=dropResults$hi_lCI[dropResults$Genus_dropped=="none"], xmax=dropResults$hi_uCI[dropResults$Genus_dropped=="none"], ymin=-Inf, ymax=Inf),
              fill="#f897a1", alpha=0.05)+
    geom_errorbar(aes(xmin = hi_lCI, xmax = hi_uCI),size=0.7, alpha=0.3, width=0.1,  position=position_jitter(seed=40), colour="black")+
    geom_point(shape=21, size=2, alpha=0.8, position=position_jitter(seed=40))+
    scale_fill_manual(values=coul)+
    labs(x="Posterior estimate of hybridicity distance", y="")+
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    geom_hline(yintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position= "none")

plot_cor[plot_cor$genetic_dist_uCI<dropResults$genetic_dist_lCI[dropResults$Genus_dropped=="none"],]

gen <- ggplot(data=plot_cor, aes(x=genetic_dist_Estimate, y=Group, fill=Group))+
    geom_rect(aes(xmin=dropResults$genetic_dist_lCI[dropResults$Genus_dropped=="none"], xmax=dropResults$genetic_dist_uCI[dropResults$Genus_dropped=="none"], ymin=-Inf, ymax=Inf),
              fill="#F8B195", alpha=0.05)+
    geom_errorbar(aes(xmin = genetic_dist_lCI, xmax = genetic_dist_uCI),size=0.7, alpha=0.3, width=0.1,  position=position_jitter(seed=40), colour="black")+
    geom_point(shape=21, size=2, alpha=0.8, position=position_jitter(seed=40))+
    scale_fill_manual(values=coul)+
        scale_x_reverse()+
    labs(x="Posterior estimate of genetic distance", y="")+
        geom_vline(xintercept=0, colour="firebrick", linetype="dashed", size=1)+
    geom_hline(yintercept=0, colour="firebrick", linetype="dashed", size=1)+
    theme_bw(base_size=10)+
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.position= "none")

Fig3_chi <- plot_grid(spa, loc, hi, gen, labels="auto")

# save
ggsave("fig/figureS3_jac.pdf", Fig3, width=170, height=200, units="mm", dpi=300)
ggsave("fig/figureS4_chi.pdf", Fig3_chi, width=170, height=200, units="mm", dpi=300)

#### residulas

