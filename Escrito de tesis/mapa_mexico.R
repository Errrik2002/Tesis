install.packages("sf")
install.packages("dafishr")

library(sf)
library(ggplot2)
library(viridis)
library(dafishr)

##
nameshp <- system.file("shape/nc.shp", package = "sf")
View(nameshp)
d <- st_read(nameshp, quiet = T)
class(d)
d$vble <- d$SID74
d$vble2 <- d$SID79

mexicomapa <- dafishr::mx_eez
st_read(mexicomapa)


class(mexicomapa)

##
ggplot(d)+geom_sf(aes(fill=vble))+
  scale_fill_viridis()+ theme_bw()


st_crs(d)$epsg
d <- st_transform(d, 4326)

##################


