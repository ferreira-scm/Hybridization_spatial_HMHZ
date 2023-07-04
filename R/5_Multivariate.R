library("ggmap")
library(sf)
library("rnaturalearthdata")
library("rnaturalearth")
library(legendMap)
library(sf)
library(maps)
library("ggspatial")
library(phyloseq)
library(ggplot2)
library(scico)
library(vegan)
library(RColorBrewer)
library(microshades)

PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")

# plotting geographic location
capital <-world.cities[c(world.cities$capital==1),]
capital[grep("Republic", capital$country.etc),]
capital <- capital[capital$country.etc%in%c("Germany", "Czech Republic"),]
world <- ne_countries(scale="medium", returnclass="sf")
europe <- ne_countries(continent="Europe", returnclass="sf", scale="medium")
metadf <- PS.TSS@sam_data
coordf <- data.frame(Lon=metadf$Longitude, Lat=metadf$Latitude)
sampling <-
    ggplot(data=europe)+
#    ggmap(area)+
    geom_sf()+
    geom_point(data=metadf,aes(x=Longitude, y=Latitude, fill=HI), size=3, alpha=0.5, shape=21)+
    coord_sf(xlim=c(6, 16), ylim=c(47, 55), expand=FALSE)+
#    geom_text(label="Berlin", data=capital, aes(x=13, y=52.3), size=5, col="Black")+
    geom_point(shape=4, data=capital, aes(x=long, y=lat), size=5, col="black", fill="firebrick")+
    scale_fill_scico("Hybrid index", palette="roma", direction=-1)+
    #scale_fill_gradient2("Hybrid\nindex", high="red", low="navy", mid="white", midpoint=0.5)+
#    xlab("Longitude", ylab="Latitude")+
#    annotation_scale(location = "bl", width_hint = 0.5) +
#    annotation_north_arrow(location = "bl", which_north = "true",,
#                                         style = north_arrow_fancy_orienteering) +
    scale_bar(lon = 6.5, lat = 54.5, arrow_length = 10, arrow_distance = 50,
                 distance_lon = 50, distance_lat = 7, distance_legend = 20,
                 dist_unit = "km", orientation = FALSE, legend_size = 3)+
    theme_bw(base_size=12)+
    theme(panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
sampling
ggplot2::ggsave(file="fig/Figure1.pdf", sampling, width = 85, height = 85, dpi = 300, units="mm")

### Permanovas
sdata <- PS.TSS@sam_data
# Jaccard distances
jac <- vegdist(PS.TSS@otu_table, method="jaccard")
# Bray Curtis dissimilarity matrix
chi <- vegdist(PS.TSS@otu_table, method="chisq")
permaJac1 <- adonis2(jac~
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality,
                   by="margin")
permaChi1 <- adonis2(chi~
                   sdata$Sex+
                   sdata$hi+
                   sdata$BMI+
                   sdata$Year+
                   sdata$Locality,
                   by="margin")
permaJac1
permaChi1

