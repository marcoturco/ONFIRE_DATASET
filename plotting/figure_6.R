library(tidyverse)
library(sf)
library(giscoR)
library(fs)
library(terra)
library(ggtext)
library(rmapshaper)
library(scales)
# Fig 6

fig6 <- dir_ls(regexp = "Fig6") %>% map(
  function(.data){
    
    load("ForFig4a.RData") 
    lonlat <- expand.grid(lon, lat)
    
    load(.data) # error en dim latlon 
    corre <- as.vector(corre)
    pval <- as.vector(pvalue)
    
    lonlat <- mutate(lonlat, corr = corre,
                     pval = pval)
    
    r <- rast(lonlat)
    crs(r) <- "EPSG:4326"
    rm(lon, lat, toplot)
    return(r)
    
  }
  
)

fig6_prj <- map(fig6, project, y = "ESRI:102008")
names(fig6_prj) <- letters[1:3]
fig6_prj <- map(fig6_prj, as.data.frame, xy = T) %>% 
  list_rbind(names_to = "subfig")

fig6_prj <- mutate(fig6_prj,
                   pval = ifelse(pval > .05, "nonsig", "sig"),
                   corr_cat = cut(ifelse(corr < 0, 0, corr), c(0, 0.20, 0.40, .60, .80, 1), include.lowest = T))


limits <- gisco_get_countries(resolution = "10") %>% 
  st_crop(xmin = -172, xmax = -50,
          ymin = 11.5, ymax = 75) %>%
  st_transform("ESRI:102008")

limits_inner <- ms_innerlines(limits)

#col_def <- RColorBrewer::brewer.pal(5, "RdYlGn")
col_def <- c("#e0f0a5", "#a8d8b2", "#73c6aa", "#499cb0", "#2670a1")
ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = fig6_prj, aes(x, y, fill = corr_cat), alpha = .8) +
  geom_point(data = filter(fig6_prj, pval == "sig"), aes(x, y, shape = "Statistical significant"), size = 2) +
  geom_sf(data = limits_inner, linewidth = .3, colour = "white") +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  facet_wrap(~subfig, ncol = 3) +
  labs(fill = NULL, shape = NULL) +
  scale_fill_manual(values = col_def) +
  scale_shape_manual(values = ".") +
  guides(fill = guide_colorsteps(barwidth = 18, 
                                 barheight = .5,
                                 title.position = "right",
                                 title.vjust = -.1,
                                 label.vjust = -5,
                                 show.limits = T),
         shape = guide_legend(override.aes = list(shape = 16),
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
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 10),
        plot.margin = margin(10, 10, 15, 10))


ggsave("Figure_6.png",
       width = 9,
       height = 5,
       units = "in",
       bg = "white",
       device = png,
       type = "cairo")

ggsave("Figure_6.pdf",
       width = 9,
       height = 5,
       units = "in",
       bg = "white",
       device = "pdf")


