library(tidyverse)
library(sf)
library(giscoR)
library(fs)
library(terra)
library(ggtext)
library(rmapshaper)
library(scales)
library(patchwork)
# Fig 5

fig5 <- dir_ls(regexp = "Fig5") %>% map(
  function(.data){
    load(.data)
    
    corre <- as.vector(corre)
    pval <- as.vector(pvalue)
    lonlat <- expand.grid(lon, lat)
    
    lonlat <- mutate(lonlat, corr = corre,
                     pval = pval)
    
    r <- rast(lonlat)
    crs(r) <- "EPSG:4326"
    rm(lon, lat, toplot)
    return(r)
    
  }
  
)


names(fig5) <- letters[1:3]

europa <- project(fig5[[3]], "EPSG:3035")
australia <- project(fig5[[1]], "EPSG:4462")
chile <- project(fig5[[2]], "EPSG:32718")


# EUROPA
europa_df <- as.data.frame(europa, xy = T)

europa_df <- mutate(europa_df,
                   pval = ifelse(pval > .05, "nonsig", "sig"),
                   corr_cat = cut(ifelse(corr < 0, 0, corr), c(0, 0.20, 0.40, .60, .80, 1), include.lowest = T))



limits <- gisco_get_countries(resolution = "10") %>% 
  st_crop(xmin = -10, xmax = 32,
          ymin = 33, ymax = 70) %>%
  st_transform("EPSG:3035")

limits_inner <- ms_innerlines(limits)

#col_def <- RColorBrewer::brewer.pal(5, "RdYlGn")
col_def <- c("#e0f0a5", "#a8d8b2", "#73c6aa", "#499cb0", "#2670a1")

eu <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = europa_df, aes(x, y, fill = corr_cat), alpha = .8) +
  geom_sf(data = limits_inner, linewidth = .3, colour = "white") +
  geom_point(data = filter(europa_df, pval == "sig"), aes(x, y, shape = "Statistical significant"), size = 1) +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  scale_fill_manual(values = col_def) +
  scale_shape_manual(values = 21) +
  guides(fill = guide_colorsteps(barwidth = 15, 
                                 barheight = 0.5,
                                 title.position = "right",
                                 title.vjust = -.1,
                                 label.vjust = -5,
                                 show.limits = T),
         shape = guide_legend(override.aes = list(shape = 21, size = 2),
                              label.position = "left",
                              order = 1)) +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(-5, "pt"),
        legend.box.margin = margin(),
        legend.title = element_markdown(),
        legend.text = element_text(size = 12),
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 50),
        plot.margin = margin(10, 10, 15, 5))


# Australia

australia_df <- as.data.frame(australia, xy = T)

australia_df <- mutate(australia_df,
                    pval = ifelse(pval > .05, "nonsig", "sig"),
                    corr_cat = cut(ifelse(corr < 0, 0, corr), c(0, 0.20, 0.40, .60, .80, 1), include.lowest = T))




limits <- gisco_get_countries(resolution = "10") %>%
  st_crop(xmin = 112, xmax = 160,
          ymin = -43, ymax = -9) %>% 
  st_transform("EPSG:4462")

austr <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = australia_df, aes(x, y, fill = corr_cat), alpha = .8) +
  geom_point(data = filter(australia_df, pval == "sig"), aes(x, y, shape = "Statistical significant"), size = 1) +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  scale_fill_manual(values = col_def) +
  scale_shape_manual(values = 21) +
  guides(fill = guide_colorsteps(barwidth = 15, 
                                 barheight = 0.5,
                                 title.position = "right",
                                 title.vjust = -.1,
                                 label.vjust = -5,
                                 show.limits = T),
         shape = guide_legend(override.aes = list(shape = 21, size = 2),
                              label.position = "left",
                              order = 1)) +
  
  
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(-5, "pt"),
        legend.box.margin = margin(),
        legend.title = element_markdown(),
        legend.text = element_text(size = 12),
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 50),
        plot.margin = margin(10, 10, 15, 1))



# CHILE

chile_df <- as.data.frame(chile, xy = T)

chile_df <- mutate(chile_df,
                       pval = ifelse(pval > .05, "nonsig", "sig"),
                       corr_cat = cut(ifelse(corr < 0, 0, corr), c(0, 0.20, 0.40, .60, .80, 1), include.lowest = T))


limits <- gisco_get_countries(resolution = "10") %>%
  st_crop(xmin = -80, xmax = -65,
          ymin = -56, ymax = -15) %>% 
  st_transform("EPSG:32718")

limits_inner <- ms_innerlines(limits)

chi <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = chile_df, aes(x, y, fill = corr_cat), alpha = .8) +
  geom_point(data = filter(chile_df, pval == "sig"), aes(x, y, shape = "Statistical significant"), size = 1) +
  geom_sf(data = limits_inner, linewidth = .3, colour = "white") +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  scale_fill_manual(values = col_def) +
  scale_shape_manual(values = 21) +
  guides(fill = guide_colorsteps(barwidth = 15, 
                                 barheight = 0.5,
                                 title.position = "right",
                                 title.vjust = -.1,
                                 label.vjust = -5,
                                 show.limits = T),
         shape = guide_legend(override.aes = list(shape = 21, size = 2),
                              label.position = "left",
                              order = 1)) +
  
  
  coord_sf() +
  theme_minimal() +
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(-5, "pt"),
        legend.box.margin = margin(),
        legend.title = element_markdown(),
        legend.text = element_text(size = 12),
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 50),
        plot.margin = margin(10, 20, 15, 5))


(austr | eu | chi) + 
  plot_annotation(tag_levels = 'a') +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(-5, "pt"),
        legend.box.margin = margin(),
        legend.title = element_markdown(),
        legend.justification = .6,
        legend.margin = margin(t = 10))


ggsave("Figure_5_.png",
       width = 16.75,
       height = 7.4,
       units = "in",
       bg = "white",
       device = png,
       type = "cairo")



