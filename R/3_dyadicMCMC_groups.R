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
Parasite <- subset_taxa(PS.TSS, Genus %in%c("Eimeria", "Cryptosporidium", "Syphacia", "Aspiculuris", "Ascaridida", "Mastophorus","Trichuris", "Hymenolepis", "Tritrichomonas"))
Fungi <- subset_taxa(PS.TSS, Phylum %in% c("Mucoromycota", "Ascomycota", "Basidiomycota"))
Diet <- subset_taxa(PS.TSS, Phylum %in% c("Anthophyta", "Phragmoplastophyta", "Charophyta", "Ochrophyta"))

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
modelJ_bac<-brm(jac_bac~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "zero_one_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelJ_bac, "tmp/BRMmodelJ_bac.rds")

modelA_bac<-brm(ait_bac~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelA_bac, "tmp/BRMmodelA_bac.rds")

modelJ_para<-brm(jac_para~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "zero_one_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelJ_para, "tmp/BRMmodelJ_para.rds")

modelA_para<-brm(ait_para~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelA_para, "tmp/BRMmodelA_para.rds")

modelJ_diet<-brm(jac_pla~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "zero_one_inflated_beta",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelJ_diet, "tmp/BRMmodelJ_diet.rds")

modelA_diet<-brm(ait_pla~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelA_diet, "tmp/BRMmodelA_diet.rds")

modelJ_fun<-brm(jac_fun~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "zero_one_inflated_beta",
                warmup = 1000, iter = 6000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelJ_fun, "tmp/BRMmodelJ_fun.rds")

modelA_fun<-brm(ait_fun~1+ spatial+HI*He+Hx+year+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelA_fun, "tmp/BRMmodelA_fun.rds")

############ effect of fungi on bacterial community

modelFB_A<-brm(ait_bac~1+ ait_fun+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 6000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelFB_A, "tmp/BRMmodelFB_A.rds")
#

modelFB_J<-brm(jac_bac~1+ jac_fun+
                (1|mm(IDA,IDB)),
                data = data.dyad,
                family= "gaussian",
                warmup = 1000, iter = 3000,
                cores = 20, chains = 4,
               inits=0)
saveRDS(modelFB_J, "tmp/BRMmodelFB_J.rds")



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

HIHeJ <- ggplot(res.df, aes(x=HI_He_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=HI_He_lCI, xmax=HI_He_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="genetic distance: hHe-dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HIHeA <- ggplot(res.dfA, aes(x=HI_He_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=HI_He_lCI, xmax=HI_He_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="genetic distance: hHe-dist", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HxJ <- ggplot(res.df, aes(x=Hx_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Hx_lCI, xmax=Hx_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="hHe-mean", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

HxA <- ggplot(res.dfA, aes(x=Hx_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=Hx_lCI, xmax=Hx_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="hHe-mean", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")

spaJ <- ggplot(res.df, aes(x=spatial_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=spatial_lCI, xmax=spatial_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Spatial distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")
spaA <- ggplot(res.dfA, aes(x=spatial_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=spatial_lCI, xmax=spatial_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Spatial distance estimate", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")
yearJ <- ggplot(res.df, aes(x=year_Estimate, y=Domain, colour=Domain))+
    geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=year_lCI, xmax=year_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Temporal distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")
yearA <- ggplot(res.dfA, aes(x=year_Estimate, y=Domain, colour=Domain))+
        geom_vline(xintercept=0, linetype="dashed", linewidth=1)+
    geom_errorbar(aes(xmin=year_lCI, xmax=year_uCI, colour=Domain),
                  size=1, width=0.4)+
    geom_point(size=3)+
#    scale_x_reverse()+
   scale_colour_manual(values=coul)+
#    scale_discrete_vi()+
    labs(x="Temporal distance", y="")+
    theme_classic(base_size=12)+
    theme(legend.position = "none")


Fig2 <- plot_grid(genJ, HeJ, HxJ, HIHeJ,
                  labels="auto", ncol=2)

FigS2 <- plot_grid(spaJ, spaA, yearJ, yearA, labels="auto") 

ggsave("fig/figure2.pdf", Fig2, width=170, height=120, units="mm", dpi=300)


############## ending analysis here for now ##################

## Figure 3: interaction of genetic

pp_check(modelA) # fine

# model convergence for jaccard
modelJ_transformed <- ggs(modelJ)

cat <- filter(modelJ_transformed, Parameter %in% c("b_Intercept", "b_spatial", "b_HI", "b_He", "b_Hx", "b_year", "b_HI:He"))
par <- c("Intercept", "Spatial distance","genetic distances", "hHe-dist", "hHe-mean", "Year distance",
         "Genetic distance*hHe-dist")
names(par) <- (unique(modelJ_transformed$Parameter))[1:7]

caterpillar <- ggplot(filter(modelJ_transformed, Parameter %in% c("b_Intercept", "b_spatial", "b_locality1", "b_HI", "b_He", "b_Hx", "b_year", "b_HI:He")),
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

#ggsave("fig/figureS1_caterpillar.pdf", caterpillar, width=170, height=250, units="mm", dpi=300)


modelA_fun

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
    xlim(min(data.dyad$He[data.dyad$HI<0.1]), max(data.dyad$He[data.dyad$HI<0.1]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.1")+
                theme_bw(base_size=12)

pred.df5 <- add_epred_draws(newdata0.5, modelA)
gen5 <-ggplot(pred.df5, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
                xlab("He")+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.5")+
                theme_bw(base_size=12)


pred.df1 <- add_epred_draws(newdata1, modelA)
gen1 <-ggplot(pred.df1, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    xlim(min(data.dyad$He[data.dyad$HI>0.9]), max(data.dyad$He[data.dyad$HI>0.9]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.9")+
                theme_bw(base_size=12)

predF0 <- add_epred_draws(newdata0, modelA_fun)
genF0 <- ggplot(predF0, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    xlim(min(data.dyad$He[data.dyad$HI<0.1]), max(data.dyad$He[data.dyad$HI<0.1]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.1")+
                theme_bw(base_size=12)

predF5 <- add_epred_draws(newdata0.5, modelA_fun)
genF5 <-ggplot(predF5, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
                xlab("He")+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.5")+
                theme_bw(base_size=12)


predF1 <- add_epred_draws(newdata1, modelA_fun)
genF1 <-ggplot(predF1, aes(x=He, y=.epred))+
    stat_lineribbon(size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    xlim(min(data.dyad$He[data.dyad$HI>0.9]), max(data.dyad$He[data.dyad$HI>0.9]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.9")+
                theme_bw(base_size=12)

All <- plot_grid(gen0, gen5, gen1, labels="auto", rel_widths=c(0.4,1,0.4), nrow=1)
Fun <- plot_grid(genF0, genF5, genF1, labels=c("d", "e", "f"), rel_widths=c(0.4,1,0.4), nrow=1)

Fig3 <- plot_grid(All, Fun, ncol=1)

ggplot2::ggsave(file="fig/Fig3.pdf", Fig3, width = 190, height = 170, dpi = 300, units="mm")

############### Fungi

newdata0 <- data.frame(He=seq_range(0:1, n=51),
                       year=rep(0, n=51),
                       HI=rep(0.1, n=51),
                       Hx=rep(median(data.dyad$Hx), n=51),
                       IDA=rep("AA_0197", 51),
                       IDB=rep("AA_0089", 51),
                       spatial=rep(median(data.dyad$spatial)))

pred.df0 <- add_epred_draws(newdata0, modelA_fun)

gen_pred0 <-ggplot(data.dyad, aes(x=He, y=ait_fun))+
    geom_jitter(width=0.01, shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df0, aes(y = .epred),
                    size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    xlim(min(data.dyad$He[data.dyad$HI<0.1]), max(data.dyad$He[data.dyad$HI<0.1]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.1")+
                theme_bw(base_size=12)

gen_pred0


newdata0.5 <- data.frame(He=seq_range(0:1, n=51),
                       year=rep(0, n=51),
                       HI=rep(0.5, n=51),
                       Hx=rep(median(data.dyad$Hx), n=51),
                       IDA=rep("AA_0197", 51),
                       IDB=rep("AA_0089", 51),
                       spatial=rep(median(data.dyad$spatial)))

pred.df5 <- add_epred_draws(newdata0.5, modelA_fun)
gen_pred0.5 <-ggplot(data.dyad, aes(x=He, y=Microbiome_similarit))+
    geom_jitter(width=0.01, shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df5, aes(y = .epred),
                    size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
                xlab("He")+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.5")+
                theme_bw(base_size=12)

newdata1 <- data.frame(He=seq_range(0:1, n=51),
                       year=rep(0, n=51),
                       HI=rep(1, n=51),
                       Hx=rep(median(data.dyad$Hx), n=51),
                       IDA=rep("AA_0197", 51),
                       IDB=rep("AA_0089", 51),
                       spatial=rep(median(data.dyad$spatial)))

pred.df1 <- add_epred_draws(newdata1, modelA_fun)

gen_pred1 <-ggplot(data.dyad, aes(x=He, y=ait_fun))+
    geom_jitter(width=0.01, shape=21, size=1, colour="gray", alpha=0.7)+
    stat_lineribbon(data=pred.df1, aes(y = .epred),
                    size=0.5, .width=c(.95, .8, .5), alpha=0.5) +
#    scale_fill_manual(values=microshades_palette("micro_purple"))+
                ylab("Gut community similarity")+
    xlab("He")+
    xlim(min(data.dyad$He[data.dyad$HI>0.9]), max(data.dyad$He[data.dyad$HI>0.9]))+
                labs(fill="level:")+
                ggtitle("Genetic distance = 0.5")+
                theme_bw(base_size=12)

gen_pred1

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


ggplot2::ggsave(file="fig/Fig4.pdf", Fig4, width = 230, height = 180, dpi = 300, units="mm")w

## Let's model
newdata0 <- data.frame(hi=seq_range(data.dyad$hi, n=51),
                      BMI=rep(median(data.dyad$BMI), n=51),
                      year=rep(0, n=51),
                      genetic_dist=rep(0, n=51),
#             hi=rep(median(data.dyad$hi), n=51),
                      locality=rep(0, n=51),
                      sex=rep("MM", 51),
                      IDA=rep("AA_0197", 51),
                      IDB=rep("AA_0089", 51),
                      spatial=rep(median(data.dyad$spatial)))
newdata1 <- data.frame(hi=seq_range(data.dyad$hi, n=51),
                      BMI=rep(median(data.dyad$BMI), n=51),
                      year=rep(0, n=51),
                      genetic_dist=rep(1, n=51),
#             hi=rep(median(data.dyad$hi), n=51),
                      locality=rep(0, n=51),
                      sex=rep("MM", 51),
                      IDA=rep("AA_0197", 51),
                      IDB=rep("AA_0089", 51),
                      spatial=rep(median(data.dyad$spatial)))
newdata05 <- data.frame(hi=seq_range(data.dyad$hi, n=51),
                      BMI=rep(median(data.dyad$BMI), n=51),
                      year=rep(0, n=51),
                      genetic_dist=rep(0.5, n=51),
#             hi=rep(median(data.dyad$hi), n=51),
                      locality=rep(0, n=51),
                      sex=rep("MM", 51),
                      IDA=rep("AA_0197", 51),
                      IDB=rep("AA_0089", 51),
                      spatial=rep(median(data.dyad$spatial)))

bac.df0 <- add_epred_draws(newdata0, model1_bac_chi)

bac0 <- ggplot(bac.df0, aes(x=hi, y=.epred))+
    stat_lineribbon(.width=c(0.5, 0.8, 0.95), colour="#136f63", fill="#136f63", alpha=0.5)+
    labs(title="Genetic distance = 0", x="Hybridicity distance", y="Predicted bacterial composition (Jaccard similarity)")+
    theme_classic(base_size=10)+
    theme(plot.title = element_text(hjust = 0.5, size = 12, face = "bold"),
               axis.title.x = element_text(size = 12),
               axis.title.y = element_text(size = 12),
               axis.text.x = element_text(size = 10),
               axis.text.y = element_text(size = 10),
               legend.position = "none")


bac0

bac.df1 <- add_epred_draws(newdata1, model1_bac_chi)
bac1 <- ggplot(bac.df1, aes(x=hi, y=.epred))+
    stat_lineribbon(.width=c(0.95), alpha=0.5)+
    labs(title="genetic distance = 1")+
    theme_classic()

bac.df05 <- add_epred_draws(newdata05, model1_bac_chi)
bac05 <- ggplot(bac.df05, aes(x=hi, y=.epred))+
    stat_lineribbon(.width=c(0.95), alpha=0.5)+
        labs(title="genetic distance = 0.5")+
    theme_classic()

plot_grid(bac0, bac05, bac1, nrow=1)


newdata0 <- data.frame(hi=seq_range(data.dyad$hi, n=51),
                      BMI=rep(median(data.dyad$BMI), n=51),
                      year=rep(0, n=51),
                      genetic_dist=rep(0, n=51),
#             hi=rep(median(data.dyad$hi), n=51),
                      locality=rep(0, n=51),
                      sex=rep("MM", 51),
                      IDA=rep("AA_0197", 51),
                      IDB=rep("AA_0089", 51),
                      spatial=rep(median(data.dyad$spatial)))
pred.df0 <- add_epred_draws(newdata0, model1_Fungi)

fun0 <- ggplot(pred.df0, aes(x=hi, y=.epred))+
    stat_lineribbon(.width=c(0.95), alpha=0.5)+
        labs(title="genetic distance = 0")+
    theme_classic()

newdata1 <- data.frame(hi=seq_range(data.dyad$hi, n=51),
                      BMI=rep(median(data.dyad$BMI), n=51),
                      year=rep(0, n=51),
                      genetic_dist=rep(1, n=51),
#             hi=rep(median(data.dyad$hi), n=51),
                      locality=rep(0, n=51),
                      sex=rep("MM", 51),
                      IDA=rep("AA_0197", 51),
                      IDB=rep("AA_0089", 51),
                      spatial=rep(median(data.dyad$spatial)))
pred.df1 <- add_epred_draws(newdata1, model1_Fungi)
fun1 <- ggplot(pred.df1, aes(x=hi, y=.epred))+
    stat_lineribbon(.width=c(0.95), alpha=0.5, colour="#ffba08")+
    labs(title="genetic distance = 1")+
    theme_classic()


newdata05 <- data.frame(hi=seq_range(data.dyad$hi, n=51),
                      BMI=rep(median(data.dyad$BMI), n=51),
                      year=rep(0, n=51),
                      genetic_dist=rep(0.5, n=51),
#             hi=rep(median(data.dyad$hi), n=51),
                      locality=rep(0, n=51),
                      sex=rep("MM", 51),
                      IDA=rep("AA_0197", 51),
                      IDB=rep("AA_0089", 51),
                      spatial=rep(median(data.dyad$spatial)))
pred.df05 <- add_epred_draws(newdata05, model1_Fungi)
fun05 <- ggplot(pred.df05, aes(x=hi, y=.epred))+
    stat_lineribbon(.width=c(0.95), alpha=0.5)+
        labs(title="genetic distance = 0.5")+
    theme_classic()


plot_grid(bac0, bac05, bac1, nrow=1)

plot_grid(fun0, fun05, fun1, nrow=1)

gen_pred0
