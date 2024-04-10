library(tidyverse)
library(sf)

# This map was created by Koki Minoji based on the Global Map Japan,
# provided by Geospatial Information Authority of Japan (GSI).
# The data were downloaded from:
# https://www.gsi.go.jp/kankyochiri/gm_jpn.html (accessed 2024-03-24)

base_map <- read_sf("./gm-jpn-all_u_2_2/polbnda_jpn.shp") |> 
  st_transform(4612)
inwater_map <- read_sf("./gm-jpn-all_u_2_2/inwatera_jpn.shp") |> 
  st_transform(4612)
coastline_map <- read_sf("./gm-jpn-all_u_2_2/coastl_jpn.shp") |> 
  st_transform(4612)
river_map <- read_sf("./gm-jpn-all_u_2_2/riverl_jpn.shp") |>
  st_transform(4612) |>
  st_zm()

locations_tbl <- read.csv("locations.coords.csv", header = FALSE)
colnames(locations_tbl) <- c("Location", "latitude", "longitude")
locations_tbl$Location <- factor(
  locations_tbl$Location,
  levels = locations_tbl$Location
)

locations <- st_as_sf(
  locations_tbl,
  coords = c("longitude", "latitude"),
  agr = c("Location")
) |> 
  st_set_crs(4612)


ggplot() +
  geom_sf(data = base_map, fill = "white", color = "white") +
  geom_sf(data = coastline_map, linewidth = 0.2) +
  geom_sf(data = river_map, color = "lightblue", linewidth = 0.2) +
  geom_sf(
    data = inwater_map,
    fill = "aliceblue",
    color = "lightblue",
    linewidth = 0.2
  ) +
  geom_sf(
    data = locations,
    mapping = aes(color = Location, shape = Location),
    size = 3
  ) +
  coord_sf(
    xlim = c(134.5, 136.5),
    ylim = c(33.5, 35.7)
  ) +
  scale_shape_manual(values = 15:18) +
  scale_x_continuous(breaks = c(135, 136)) +
  scale_y_continuous(breaks = c(33.5, 34.5, 35.5)) +
  labs(x = "Longitude", y = "Latitude") +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_line(linewidth = 0.4),
    panel.border = element_rect(fill = NA, linewidth = 0.5),
    panel.background = element_rect(fill = "aliceblue"),
    legend.key = element_rect(fill = "white"),
    text = element_text(size = 14),
    legend.position = "none"
  )

ggsave(
  "locations_map.png",
   scale = 1, dpi = 400, width = 3.5, height = 3.5
)
