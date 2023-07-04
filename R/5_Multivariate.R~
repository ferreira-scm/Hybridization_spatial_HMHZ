library("ggmap")
library(sf)
library("rnaturalearthdata")
library("rnaturalearth")

library(legendMap)

library("legendMap")

library(sf)

library(maps)


head(world.cities)

capital <-world.cities[c(world.cities$capital==1),]

capital <- capital[capital$country.etc=="Germany",]


    world.cities %>%
      filter(country.etc %in% capital == 1)

PS.TSS <- readRDS("tmp/PS.TSS_filtered.rds")


world <- ne_countries(scale="medium", returnclass="sf")
class(world)

europe <- ne_countries(continent="Europe", returnclass="sf", scale="medium")

metadf <- PS.TSS@sam_data

coordf <- data.frame(Lon=metadf$Longitude, Lat=metadf$Latitude)

#sites <- st_as_sf(coordf, coords = c("Lon", "Lat"),
#                             crs = 4326, agr = "constant")

(capital)

world$name

library("ggspatial")

sampling <-
    ggplot(data=europe)+
#    ggmap(area)+
    geom_sf()+
    geom_point(data=metadf,aes(x=Longitude, y=Latitude, fill=HI), size=2, alpha=0.5, shape=21)+
    coord_sf(xlim=c(8, 16), ylim=c(46, 56), expand=FALSE)+
    geom_text(label="Berlin", data=capital, aes(x=13, y=52.3), size=5, col="Black")+
    geom_point(shape=4, data=capital, aes(x=long, y=lat), size=5, col="black", fill="firebrick")+
    scale_fill_scico("Hybrid index", palette="roma", direction=-1)+
    #scale_fill_gradient2("Hybrid\nindex", high="red", low="navy", mid="white", midpoint=0.5)+
#    xlab("Longitude", ylab="Latitude")+
#    annotation_scale(location = "bl", width_hint = 0.5) +
#    annotation_north_arrow(location = "bl", which_north = "true",,
#                                         style = north_arrow_fancy_orienteering) +
    scale_bar(lon = 9, lat = 55, arrow_length = 10, arrow_distance = 50,
                 distance_lon = 50, distance_lat = 7, distance_legend = 20,
                 dist_unit = "km", orientation = FALSE, legend_size = 3)+
    theme_bw(base_size=12)+
    theme(legend.position="top",
          panel.border = element_blank(), panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), axis.line = element_line(colour = "black"))
sampling

saveRDS(sampling, "tmp/map.RDS")

#saveRDS(sampling, "tmp/sampling_map.rds")
ggplot2::ggsave(file="fig/Figure1.pdf", sampling, width = 120, height = 120, dpi = 300, units="mm")
