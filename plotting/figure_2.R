library(tidyverse)
library(sf)
library(giscoR)
library(fs)
library(terra)
library(ggtext)
library(rmapshaper)
library(scales)
# Fig 2

fig2 <- dir_ls(regexp = "Fig2") %>% map(
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

fig2 <- rast(fig2)
names(fig2) <- letters[1:5]
crs(fig2) <- "EPSG:4326"



fig2 <- project(fig2, "ESRI:102008")


#
fig2_df <- as.data.frame(fig2, xy = T)
fig2_df <- pivot_longer(fig2_df, 3:7, names_to = "subfig", values_to = "val")

limits <- gisco_get_countries(resolution = "10") %>% 
            st_crop(xmin = -172, xmax = -50,
                    ymin = 11.5, ymax = 75) %>%
             st_transform("ESRI:102008")

limits_inner <- ms_innerlines(limits)

ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = fig2_df, aes(x, y, fill = cut(val, c(0, 1, 2.5, 5, 10, 25, 50, 100, 250, 500, 7000)))) +
  geom_sf(data = limits_inner, linewidth = .3, colour = "white") +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  facet_wrap(~subfig, ncol = 3) +
  labs(fill = "km<sup>2</sup>") +
  scale_fill_viridis_d(option = "magma", direction = -1,
                       label = function(lab){ ifelse(lab == 2.5, number(lab, accuracy = .1),
                                                     number(lab, accuracy = 1))}) +
  guides(fill = guide_colorsteps(barwidth = 18, 
                                 barheight = .5,
                                 title.position = "right",
                                 title.vjust = -.1)) +
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.title = element_markdown(),
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 10))
      


ggsave("Figure_2.png",
       width = 9,
       height = 5,
       units = "in",
       bg = "white",
       device = png,
       type = "cairo")

