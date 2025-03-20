install.packages("sf")
install.packages("dafishr")

library(sf)
library(ggplot2)
library(viridis)
library(dafishr)
library(maps)
library(sf)
library(viridis)
library(leaflet)
d <- world.cities
# Select South Africa
d <- d[which(d$country.etc == "Mexico"), ]
# Transform data to sf object
d <- st_as_sf(d, coords = c("long", "lat"))
# Assign CRS
st_crs(d) <- 4326

d$vble <- d$pop
d$size <- sqrt(d$vble)/100
ggplot(d) + geom_sf(aes(col = vble, size = size)) +
  scale_color_viridis()

pal <- colorNumeric(palette = "viridis", domain = d$vble)
leaflet(d) %>% addTiles() %>%
  addCircles(lng = st_coordinates(d)[, 1],
             lat = st_coordinates(d)[, 2],
             radius = ~sqrt(vble)*10,
             color = ~pal(vble), popup = ~name) %>%
  addLegend(pal = pal, values = ~vble, position = "bottomright")



