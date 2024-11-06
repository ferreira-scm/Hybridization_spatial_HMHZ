# This script follows the super cool approach from Aura Raulo
# https://github.com/nuorenarra/Analysing-dyadic-data-with-brms/blob/main/R_Making_dyadic_data/DYADIC_workshop_data_wrangling.Rmd

library(brms)
library(rstan)
library(RColorBrewer) # needed for some extra colours in one of the graphs
library(vegan)
library(phyloseq)
library(ggplot2)
library(cowplot)

### we don't include sex because it does not converge
PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

Bac <- subset_taxa(PS.TSS, Kingdom %in%"Bacteria")
Parasite <- subset_taxa(PS.TSS, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))
Fungi <- subset_taxa(PS.TSS, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota"))
Diet <- subset_taxa(PS.TSS, Phylum %in% c("Anthophyta", "Phragmoplastophyta", "Charophyta", "Ochrophyta"))

get_taxa_unique(PS.TSS, "Kingdom")

# How many annotated genera?
Euk <- subset_taxa(PS.TSS, Kingdom %in%"Eukarya")
length(get_taxa_unique(Euk, "Genus"))-length(grep("Unknown", get_taxa_unique(Euk, "Genus")))
length(get_taxa_unique(Bac, "Genus"))-length(grep("Unknown", get_taxa_unique(Bac, "Genus")))

############# First create dyad data#######################
data.dyad <- readRDS("tmp/data.dyad.RDS")

names(data.dyad)

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
data.dyad$jac_bac <- jac_bac
data.dyad$ait_bac <- ait_bac
data.dyad$jac_para <- jac_para
data.dyad$ait_para <- ait_para
data.dyad$jac_fun <- jac_fun
data.dyad$ait_fun <- ait_fun
data.dyad$jac_pla <- jac_pla
data.dyad$ait_pla <- ait_pla

#now model
#modelJ_bac<-brm(jac_bac~1+ spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "zero_one_inflated_beta",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelJ_bac, "tmp/BRMmodelJ_bac.rds")

#modelA_bac<-brm(ait_bac~1+ spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelA_bac, "tmp/BRMmodelA_bac.rds")

#modelJ_para<-brm(jac_para~1+ spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "zero_one_inflated_beta",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelJ_para, "tmp/BRMmodelJ_para.rds")

#modelA_para<-brm(ait_para~1+ spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelA_para, "tmp/BRMmodelA_para.rds")

#modelJ_diet<-brm(jac_pla~1+ spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "zero_one_inflated_beta",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelJ_diet, "tmp/BRMmodelJ_diet.rds")

#modelA_diet<-brm(ait_pla~1+ spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelA_diet, "tmp/BRMmodelA_diet.rds")

#modelJ_fun<-brm(jac_fun~1+ spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "zero_one_inflated_beta",
#                warmup = 1000, iter = 6000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelJ_fun, "tmp/BRMmodelJ_fun.rds")

#modelA_fun<-brm(ait_fun~1+ spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelA_fun, "tmp/BRMmodelA_fun.rds")

############ effect of fungi on bacterial community
#modelFB_A<-brm(ait_bac~1+ ait_fun + spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelFB_A, "tmp/BRMmodelFB_A.rds")
#

#modelFB_J<-brm(jac_bac~1+ jac_fun + spatial+HI*He+Hx+year+
#                (1|mm(IDA,IDB)),
#                data = data.dyad,
#                family= "gaussian",
#                warmup = 1000, iter = 3000,
#                cores = 20, chains = 4,
#               inits=0)
#saveRDS(modelFB_J, "tmp/BRMmodelFB_J.rds")

#################################
### uploading models
modelA <- readRDS("tmp/BRMmodelA.rds")
modelJ <- readRDS("tmp/BRMmodelJac.rds")

modelJ_fun <- readRDS("tmp/BRMmodelJ_fun.rds")
modelA_fun <- readRDS("tmp/BRMmodelA_fun.rds")
modelJ_diet <- readRDS("tmp/BRMmodelJ_diet.rds")
modelA_diet <- readRDS("tmp/BRMmodelA_diet.rds")
modelA_para <- readRDS("tmp/BRMmodelA_para.rds")
modelJ_para <- readRDS("tmp/BRMmodelJ_para.rds")
modelJ_bac <- readRDS("tmp/BRMmodelJ_bac.rds")
modelA_bac <- readRDS("tmp/BRMmodelA_bac.rds")

modelFB_J <- readRDS("tmp/BRMmodelFB_J.rds")
modelFB_A <- readRDS("tmp/BRMmodelFB_A.rds")


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
               Hx_Estimate=para[rownames(para)=="Hx", "Estimate"],
               Hx_lCI=para[rownames(para)=="Hx", "l-95% CI"],
               Hx_uCI=para[rownames(para)=="Hx", "u-95% CI"],
               spatial_Estimate=para$Estimate[rownames(para)=="spatial"],
               spatial_lCI=para[rownames(para)=="spatial", "l-95% CI"],
               spatial_uCI=para[rownames(para)=="spatial", "u-95% CI"],
               year_Estimate=para$Estimate[rownames(para)=="year"],
               year_lCI=para[rownames(para)=="year", "l-95% CI"],
               year_uCI=para[rownames(para)=="year", "u-95% CI"],
               HI_He_Estimate=para$Estimate[rownames(para)=="HI:He"],
               HI_He_lCI=para[rownames(para)=="HI:He", "l-95% CI"],
               HI_He_uCI=para[rownames(para)=="HI:He", "u-95% CI"]
               )
}

res.df <-resdf.fun(modelJ_para, "Parasite", 11)
res.df <- rbind(res.df, resdf.fun(modelJ_bac, "Bacteria", 383))
res.df <- rbind(res.df, resdf.fun(modelJ_diet, "Diet", 45))
res.df <- rbind(res.df, resdf.fun(modelJ_fun, "Fungi", 65))
res.df <- rbind(res.df, resdf.fun(modelJ, "Full model", 588))
res.df$Domain <- factor(res.df$Domain, level=c( "Diet", "Bacteria","Parasite", "Fungi", "Full model"))

res.dfA <-resdf.fun(modelA_para, "Parasite", 11)
res.dfA <- rbind(res.dfA, resdf.fun(modelA_bac, "Bacteria", 383))
res.dfA <- rbind(res.dfA, resdf.fun(modelA_diet, "Diet", 45))
res.dfA <- rbind(res.dfA, resdf.fun(modelA_fun, "Fungi", 65))
res.dfA <- rbind(res.dfA, resdf.fun(modelA, "Full model", 588))
res.dfA$Domain <- factor(res.dfA$Domain, level=c( "Diet", "Bacteria", "Parasite", "Fungi", "Full model"))

coul <- c("#136f63", "#032b43", "#3f88c5", "#ffba08", "#d00000")

genJ <- ggplot(res.df, aes(x=HI_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=HI_lCI, xmax=HI_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Subspecies' genetic distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

genA <- ggplot(res.dfA, aes(x=HI_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=HI_lCI, xmax=HI_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Genetic distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HeJ <- ggplot(res.df, aes(x=He_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=He_lCI, xmax=He_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="hHe-distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HeA <- ggplot(res.dfA, aes(x=He_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=He_lCI, xmax=He_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="hHe-dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HIHeJ <- ggplot(res.df, aes(x=HI_He_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=HI_He_lCI, xmax=HI_He_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="genetic distance: hHe-dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HIHeA <- ggplot(res.dfA, aes(x=HI_He_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=HI_He_lCI, xmax=HI_He_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="genetic distance: hHe-dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HxJ <- ggplot(res.df, aes(x=Hx_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=Hx_lCI, xmax=Hx_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="hHe-mean", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HxA <- ggplot(res.dfA, aes(x=Hx_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=Hx_lCI, xmax=Hx_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="hHe-mean", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

spaJ <- ggplot(res.df, aes(x=spatial_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=spatial_lCI, xmax=spatial_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Spatial distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")
spaA <- ggplot(res.dfA, aes(x=spatial_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=spatial_lCI, xmax=spatial_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Spatial distance estimate", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")
yearJ <- ggplot(res.df, aes(x=year_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=year_lCI, xmax=year_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Temporal distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")
yearA <- ggplot(res.dfA, aes(x=year_Estimate, y=Domain, colour=Domain))+
        geom_vline(xintercept=0, linetype="dashed", linewidth=0.9)+
    geom_errorbar(aes(xmin=year_lCI, xmax=year_uCI, colour=Domain),
                  size=1, width=0.3)+
    geom_point(size=2)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Temporal distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")


Fig2 <- plot_grid(genJ, genA, HeJ, HeA, HxJ, HxA, spaJ, spaA, yearJ, yearA,
                  labels="auto", ncol=2)

#FigS1 <- plot_grid(spaJ, spaA, yearJ, yearA, labels="auto") 

ggsave("fig/figure2.pdf", Fig2, width=170, height=200, units="mm", dpi=300)

#ggsave("fig/figureS1.pdf", FigS1, width=170, height=150, units="mm", dpi=300)


library(modelr)
library(tidybayes)

######################## Figure 3: interaction
newdata0 <- data.frame(He=seq_range(0:1, n=51),
                       year=rep(0, n=51),
                       HI=rep(0.1, n=51),
                       Hx=rep(median(data.dyad$Hx), n=51),
                       IDA=rep("AA_0197", 51),
                       IDB=rep("AA_0089", 51),
                       spatial=rep(median(data.dyad$spatial)))

newdata0.5 <- data.frame(He=seq_range(0:1, n=51),
                       year=rep(0, n=51),
                       HI=rep(0.5, n=51),
                       Hx=rep(median(data.dyad$Hx), n=51),
                       IDA=rep("AA_0197", 51),
                       IDB=rep("AA_0089", 51),
                       spatial=rep(median(data.dyad$spatial)))

newdata0.63 <- data.frame(He=seq_range(0:1, n=51),
                       year=rep(0, n=51),
                       HI=rep(0.63, n=51),
                       Hx=rep(median(data.dyad$Hx), n=51),
                       IDA=rep("AA_0197", 51),
                       IDB=rep("AA_0089", 51),
                       spatial=rep(median(data.dyad$spatial)))

newdata1 <- data.frame(He=seq_range(0:1, n=51),
                       year=rep(0, n=51),
                       HI=rep(0.9, n=51),
                       Hx=rep(median(data.dyad$Hx), n=51),
                       IDA=rep("AA_0197", 51),
                       IDB=rep("AA_0089", 51),
                       spatial=rep(median(data.dyad$spatial)))

pred.df0 <- add_epred_draws(newdata0, modelA)
gen0 <- ggplot(pred.df0, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
    ylab("Gut community similarity")+
    xlab("He")+
    ylim(-2.03, -1.85)+
    xlim(min(data.dyad$He[data.dyad$HI<0.1]), max(data.dyad$He[data.dyad$HI<0.1]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.1")+
                theme_bw(base_size=10)

pred.df5 <- add_epred_draws(newdata0.5, modelA)
gen5 <-ggplot(pred.df5, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    ylim(-2.03, -1.85)+
    labs(fill="level:")+
    ggtitle("Genetic distance = 0.5")+
    theme_bw(base_size=10)

pred.df63 <- add_epred_draws(newdata0.63, modelA)
gen63 <-ggplot(pred.df63, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    ylim(-2.03, -1.85)+
    labs(fill="level:")+
    ggtitle("Genetic distance = 0.63")+
    theme_bw(base_size=10)


pred.df1 <- add_epred_draws(newdata1, modelA)
gen1 <-ggplot(pred.df1, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    ylim(-2.03, -1.85)+
    xlim(min(data.dyad$He[data.dyad$HI>0.9]), max(data.dyad$He[data.dyad$HI>0.9]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.9")+
    theme_bw(base_size=10)

predF0 <- add_epred_draws(newdata0, modelA_fun)
genF0 <- ggplot(predF0, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    ylim(0.3, 0.525)+
    xlim(min(data.dyad$He[data.dyad$HI<0.1]), max(data.dyad$He[data.dyad$HI<0.1]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.1")+
                theme_bw(base_size=10)

predF5 <- add_epred_draws(newdata0.5, modelA_fun)
genF5 <-ggplot(predF5, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
    ylim(0.3, 0.525)+
                ylab("Gut community similarity")+
                xlab("He")+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.5")+
                theme_bw(base_size=10)

predF63 <- add_epred_draws(newdata0.63, modelA_fun)
genF63 <-ggplot(predF63, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
    ylim(0.3, 0.525)+
                ylab("Gut community similarity")+
                xlab("He")+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.63")+
                theme_bw(base_size=10)


predF1 <- add_epred_draws(newdata1, modelA_fun)
genF1 <-ggplot(predF1, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    ylim(0.3, 0.525)+
    xlim(min(data.dyad$He[data.dyad$HI>0.9]), max(data.dyad$He[data.dyad$HI>0.9]))+
    labs(fill="level:")+
    ggtitle("Genetic distance = 0.9")+
    theme_bw(base_size=10)

All <- plot_grid(gen0, gen5, gen1, labels="auto", rel_widths=c(0.7,1,0.7), nrow=1)
Fun <- plot_grid(genF0, genF5, genF1, labels=c("d", "e", "f"), rel_widths=c(0.7,1,0.7), nrow=1)

Fig3 <- plot_grid(All, Fun, ncol=1)

ggplot2::ggsave(file="fig/Fig3.pdf", Fig3, width = 190, height = 140, dpi = 300, units="mm")

summary(data.dyad$HI>0.9)

########################Figure 4 ab
# plotting bacteria~fungi
nd <- data.frame(jac_fun=seq_range(data.dyad$jac_fun, n=51),
                 He=rep(median(data.dyad$He), n=51),
                 year=rep(0, n=51),
                 HI=rep(0.9, n=51),
                 Hx=rep(median(data.dyad$Hx), n=51),
                 IDA=rep("AA_0197", 51),
                 IDB=rep("AA_0089", 51),
                 spatial=rep(median(data.dyad$spatial)))
                 
add_epred_draws(nd, modelFB_J) %>%
    ggplot(aes(x = jac_fun, y = jac_bac)) +
    stat_lineribbon(aes(y = .epred), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    xlab("Fungi composition similarity")+
    ylab("Bacteria composition similarity")+
    theme_bw(base_size=10)-> Bac_funJ



nd2 <- data.frame(ait_fun=seq_range(data.dyad$ait_fun, n=51),
                 He=rep(median(data.dyad$He), n=51),
                 year=rep(0, n=51),
                 HI=rep(0.9, n=51),
                 Hx=rep(median(data.dyad$Hx), n=51),
                 IDA=rep("AA_0197", 51),
                 IDB=rep("AA_0089", 51),
                 spatial=rep(median(data.dyad$spatial)))
                 
add_epred_draws(nd2, modelFB_A) %>%
    ggplot(aes(x = ait_fun, y = ait_bac)) +
    stat_lineribbon(aes(y = .epred), .width = c(.95, .80, .50),  # regression line and CI
                    alpha = 0.5, colour = "black") +
    xlab("Fungi composition similarity")+
    ylab("Bacteria composition similarity")+
    theme_bw(base_size=10)-> Bac_funA

library(gridExtra)

Fig4ab <- grid.arrange(Bac_funJ, Bac_funA, ncol=2)

ggplot2::ggsave(file="fig/Fig4ab.pdf", Fig4ab, width = 185, height = 65, dpi = 300, units="mm")

