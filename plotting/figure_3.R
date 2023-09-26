library(tidyverse)
library(sf)
library(giscoR)
library(fs)
library(terra)
library(ggtext)
library(rmapshaper)
library(scales)
library(patchwork)
# Fig 3

fig3 <- dir_ls(regexp = "Fig3") %>% map(
  function(.data){
    load(.data)
    
    val <- as.vector(toplot)
    lonlat <- expand.grid(lon, lat)
    
    lonlat <- mutate(lonlat, val = val)
    
    r <- rast(lonlat)
    rm(lon, lat, toplot)
    return(r)
    
  }
  
)


names(fig3) <- letters[1:6]

australia <- rast(fig3[1:2])
chile <- rast(fig3[3:4])
europa <- rast(fig3[5:6])

crs(australia) <- "EPSG:4326"
crs(chile) <- "EPSG:4326"
crs(europa) <- "EPSG:4326"

europa <- project(europa, "EPSG:3035")
australia <- project(australia, "EPSG:4462")
chile <- project(chile, "EPSG:32718")

# EUROPA
europa_df <- as.data.frame(europa, xy = T)
europa_df <- pivot_longer(europa_df, 3:4, names_to = "subfig", values_to = "val")

limits <- gisco_get_countries(resolution = "10") %>% 
  st_crop(xmin = -15, xmax = 34,
          ymin = 33, ymax = 75) %>%
  st_transform("EPSG:3035")

limits_inner <- ms_innerlines(limits)

eu <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = europa_df, aes(x, y, fill = cut(val, c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 7000)))) +
  geom_sf(data = limits_inner, linewidth = .3, colour = "white") +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  facet_wrap(~subfig, ncol = 3) +
  labs(fill = "km<sup>2</sup>") +
  scale_fill_viridis_d(option = "magma", direction = -1, drop = F,
                       label = function(lab){ ifelse(lab == 2.5, number(lab, accuracy = .1),
                                                     number(lab, accuracy = 1))}) +
  guides(fill = guide_colorsteps(barwidth = 40, 
                                 barheight = 1,
                                 title.position = "right",
                                 title.vjust = -.1)) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_markdown(),
        legend.text = element_text(size = 22),  # Cambia 14 por el tamaño de texto
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 22))

# Australia

australia_df <- as.data.frame(australia, xy = T)
australia_df <- pivot_longer(australia_df, 3:4, names_to = "subfig", values_to = "val")

limits <- gisco_get_countries(resolution = "10") %>%
         st_crop(xmin = 112, xmax = 160,
                 ymin = -43, ymax = -9) %>% 
            st_transform("EPSG:4462")

austr <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = australia_df, aes(x, y, fill = cut(val, c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 7000)))) +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  labs(fill = NULL, shape = NULL) +
  facet_wrap(~subfig, ncol = 3) +
  labs(fill = "km<sup>2</sup>") +
  scale_fill_viridis_d(option = "magma", direction = -1,
                       label = function(lab){ ifelse(lab == 2.5, number(lab, accuracy = .1),
                                                     number(lab, accuracy = 1))}) +
  guides(fill = guide_colorsteps(barwidth = 40, 
                                 barheight = 1,
                                 title.position = "right",
                                 title.vjust = -.1)) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_markdown(),
        legend.text = element_text(size = 22),  # Cambia 14 por el tamaño de texto
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 22))


# CHILE

chile_df <- as.data.frame(chile, xy = T)
chile_df <- pivot_longer(chile_df, 3:4, names_to = "subfig", values_to = "val")

limits <- gisco_get_countries(resolution = "10") %>%
  st_crop(xmin = -80, xmax = -65,
          ymin = -56, ymax = -15) %>% 
  st_transform("EPSG:32718")

limits_inner <- ms_innerlines(limits)

chi <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = chile_df, aes(x, y, fill = cut(val, c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 7000)))) +
  geom_sf(data = limits_inner, linewidth = .3, colour = "white") +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
   facet_wrap(~subfig, ncol = 3) +
  labs(fill = "km<sup>2</sup>") +
  scale_fill_viridis_d(option = "magma", direction = -1, drop = F,
                       label = function(lab){ ifelse(lab == 2.5, number(lab, accuracy = .1),
                                                     number(lab, accuracy = 1))}) +
  guides(fill = guide_colorsteps(barwidth = 40, 
                                 barheight = 1,
                                 title.position = "right",
                                 title.vjust = -.1)) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_markdown(),
        legend.text = element_text(size = 22),  # Cambia 14 por el tamaño de texto
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 22))




(austr/eu | chi) + 
  plot_layout(guides = "collect",
              widths = c(10, 5.5)) &
  theme(legend.position = "bottom")
  

ggsave("Figure_3.png",
       width = 17,
       height = 12,
       units = "in",
       bg = "white",
       device = png,
       type = "cairo")

ggsave("Figure_3.pdf",
       width = 17,
       height = 12,
       units = "in",
       bg = "white",
       device = "pdf")




