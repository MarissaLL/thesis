library(tidyverse)
library(cowplot)
library(osmplotr)
library(osmdata)


# Data from Williams 1956. Locations are very approximate, and largely just based on placenames
williams1956 <- read_delim("~/projects/kakapo-genomics/data/all_historic_kakapo_williams56.csv", delim = ",") %>% 
  filter(type != "Introduction")

# Data from Lentini 2018
lentini2018 <- read_delim("~/projects/kakapo-genomics/data/kakapo_sites_lentini_tidied.csv", delim = ",")

currentpop <- read_delim("~/projects/kakapo-genomics/data/current_kakapo_sites.csv", delim = ",", col_names = c("site", "latitude", "longitude", "type" ))

all_sites <- bind_rows(williams1956, lentini2018, currentpop)

# Load in map from maps package
# nz_map <- map_data("nz")

# Set inset map details ranges to download from OSM (bbox) and to plot (plotrange)

main_bbox <- c(top = -33.980, left = 164.993, right = 179.714, bottom = -47.769)

f_bbox <- c(top = -46.6683, left = 167.5031, right = 168.3215, bottom = -46.9540)
f_plotrange_x <- c(166.4, 167.0)
f_plotrange_y <- c(-45.6, -46.2)

r_bbox <- c(top = -46.6746, left = 167.0, right = 168.3054, bottom = -47.5)
r_plotrange_x <- c(167.5, 168.3)
r_plotrange_y <-  c(-46.67, -47.22)

hot_bbox <- c(top = -35.8156, left = 173.9905, right = 176.0738, bottom = -36.8829)
hot_plotrange_x <- c(174.1, 176)
hot_plotrange_y <- c(-35.9, -36.8)



# Download map data from OSM and reformat for plotting
# Have to use osmdata here for the main map because osmplotr query times out
# osmdata sometimes times out too, but usually works after a couple of retries
osm_main <- opq(main_bbox) %>%
  add_osm_feature(key = "natural", value = "coastline") %>%
  osmdata_sp(quiet = FALSE)

main_ff <- fortify(osm_main$osm_lines) 


# Using osmplotr for inset maps
osm_fiordland <-  extract_osm_objects(bbox = f_bbox, key = "boundary", sf = FALSE)
fiordland_ff <- fortify(osm_fiordland)

osm_rakiura <- extract_osm_objects(bbox =r_bbox , key = "natural", value = "coastline", return_type = "line", sf = FALSE)
rakiura_ff <- fortify(osm_rakiura) 

osm_hauturu_o_toi <- extract_osm_objects(bbox =hot_bbox , key = "natural", value = "coastline", return_type = "line", sf = FALSE)
hauturu_o_toi_ff <- fortify(osm_hauturu_o_toi) 




# Kakapo colours from manu: dark green "#7D9D33" light green "#CED38C" yellow "#DCC949" grey-brown "#BCA888" red-brown"#CD8862" dark-brown "#775B24"
# Plot the map
colours_order <- c("#775B24", "#7D9D33", "#CED38C") # Old = green, current = brown
# colours_order <- c("#7D9D33", "#775B24", "#BCA888") # Old = brown, current = green
# colours_order <- c("#CD8862", "#7D9D33", "#CED38C") # Old = green, current = red
# colours_order <- c("#CD8862", "#775B24", "#BCA888") # Old = brown, current = red


main_map <- ggplot() + 
  geom_path(data = main_ff, 
            mapping = aes(x=long, 
                          y = lat, 
                          group = group), 
            colour = 'black' ) +
  coord_fixed(1.4) + 
  geom_point(data = all_sites,
             mapping = aes(x = longitude,
                           y = latitude, 
                           colour = type),
             size = 3.5, alpha = 0.7) +
  scale_color_manual(name = c("Current population", "Historical sighting", "Subfossil remains"), aesthetics = "colour", values = colours_order) +
  geom_point(data = subset(all_sites, type == "Current population"),
             mapping = aes(x = longitude,
                           y = latitude),
             colour = colours_order[1],
             size = 3.5, alpha = 0.7)+
  geom_rect(aes(xmax = f_plotrange_x[1], xmin = f_plotrange_x[2], ymin = f_plotrange_y[1], ymax = f_plotrange_y[2]), fill = NA, colour = "red") +
  geom_rect(aes(xmax = r_plotrange_x[1], xmin = r_plotrange_x[2], ymin = r_plotrange_y[1], ymax = r_plotrange_y[2]), fill = NA, colour = "red") +
  geom_rect(aes(xmax = hot_plotrange_x[1], xmin = hot_plotrange_x[2], ymin = hot_plotrange_y[1], ymax = hot_plotrange_y[2]), fill = NA, colour = "red") +
  theme_classic() +
  theme(legend.position = c(0.8, 0.3),
        axis.line.x.top = element_line(),
        axis.line.y.right = element_line(),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ),
        axis.title = element_text(size = 18, margin = margin(r = 12)),
        axis.text = element_text(size = 12),
        legend.title = element_blank()) +
  labs(x = 'Longitude',
       y = 'Latitude')


fiordland_inset <- ggplot() +
  geom_polygon(data = fiordland_ff , aes(x = long, y = lat, group = group), fill = NA, color = 'black')+
  geom_point(data = all_sites,
             mapping = aes(x = longitude,
                           y = latitude,
                           colour = type),
             size = 3.5, alpha = 0.7) +
  scale_color_manual(name = c("Current population", "Historical sighting", "Subfossil remains"), aesthetics = "colour", values = colours_order) +
  coord_fixed(1.4, xlim = f_plotrange_x, ylim = f_plotrange_y)+
  theme_classic()+ 
  theme(legend.position = "none",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ))



rakiura_inset <- ggplot() +
  geom_path(data = rakiura_ff , aes(x = long, y = lat, group = group), color = 'black')+
  geom_point(data = all_sites,
             mapping = aes(x = longitude,
                           y = latitude,
                           colour = type),
             size = 3.5, alpha = 0.7) +
  scale_color_manual(name = c("Current population", "Historical sighting", "Subfossil remains"), aesthetics = "colour", values = colours_order) +
  coord_fixed(1.4, xlim = r_plotrange_x, ylim = r_plotrange_y)+
  theme_classic()+ 
  theme(legend.position = "blank",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ))



hot_inset <- ggplot() +
  geom_path(data = hauturu_o_toi_ff , aes(x = long, y = lat, group = group), color = 'black')+
  geom_point(data = all_sites,
             mapping = aes(x = longitude,
                           y = latitude,
                           colour = type),
             size = 3.5, alpha = 0.7) +
  scale_color_manual(name = c("Current population", "Historical sighting", "Subfossil remains"), aesthetics = "colour", values = colours_order) +
  coord_fixed(1.4, xlim = hot_plotrange_x, ylim = hot_plotrange_y)+
  theme_classic()+ 
  theme(legend.position = "blank",
        axis.line = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(),
        axis.title = element_blank(),
        panel.background = element_rect(fill = NA,
                                        colour = "black",
                                        size = 0.7 ))


full_map <-  ggdraw() +
  draw_plot(main_map) +
  draw_plot(fiordland_inset, x = 0.15, y = 0.17, width = 0.32)+
  draw_plot(rakiura_inset, x = 0.53, y = -0.315, width = 0.32)+
  draw_plot(hot_inset, x = 0.15, y = 0.4, width = 0.32)

full_map

ggsave("~/Pictures/kakapo_map_4.png", full_map, width = 15, height = 20, units = "cm")

