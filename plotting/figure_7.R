library(tidyverse)
library(sf)
library(giscoR)
library(fs)
library(terra)
library(ggtext)
library(rmapshaper)
library(scales)
library(patchwork)
# Fig 7

fig7 <- dir_ls(regexp = "Fig7") %>% map(
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


names(fig7) <- letters[1:6]
europa <- fig7[c(3,6)]
europa <- map(europa, project, y = "EPSG:3035")

australia <- fig7[c(1,4)]
australia <- map(australia, project, y = "EPSG:4462")

chile <- fig7[c(2,5)]
chile <- map(chile, project, y = "EPSG:32718")

# EUROPA
europa_df <- as.data.frame(rast(europa), xy = T)
names(europa_df)[3:6] <- c("c_corr", "c_pval",
                           "f_corr", "f_pval")

europa_df <- pivot_longer(europa_df,
             3:6,
             names_to = c("subfig", "var"),
             values_to = "val",
             names_sep = "_") %>%
             pivot_wider(names_from = var, values_from = val)


europa_df <- mutate(europa_df,
                    pval = ifelse(pval > .05, "nonsig", "sig"),
                    corr_cat = cut(ifelse(corr < 0, 0, corr), c(0, 0.20, 0.40, .60, .80, 1), include.lowest = T)) %>%
              filter(!is.na(corr))



limits <- gisco_get_countries(resolution = "10") %>% 
  st_crop(xmin = -15, xmax = 34,
          ymin = 33, ymax = 75) %>%
  st_transform("EPSG:3035")

limits_inner <- ms_innerlines(limits)

#col_def <- RColorBrewer::brewer.pal(5, "RdYlGn")
col_def <- c("#e0f0a5", "#a8d8b2", "#73c6aa", "#499cb0", "#2670a1")

eu <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = europa_df, aes(x, y, fill = corr_cat), alpha = .8)  +
  geom_sf(data = limits_inner, linewidth = .3, colour = "white") +
  geom_point(data = filter(europa_df, pval == "sig"), aes(x, y, shape = "Statistical significant"), size = .8) +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  facet_wrap(~subfig, ncol = 1) +
  labs(fill = NULL, shape = NULL) +
  scale_fill_manual(values = col_def) +
  scale_shape_manual(values = 21) +
  guides(fill = guide_colorsteps(barwidth = 15, 
                                 barheight = .5,
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
        legend.text = element_text(size = 14),  # Cambia 14 por el tamaño de texto
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 20),
        plot.margin = margin(5, 5, 15, 10))


# Australia

australia_df <- as.data.frame(rast(australia), xy = T)
names(australia_df)[3:6] <- c("a_corr", "a_pval",
                           "d_corr", "d_pval")

australia_df <- pivot_longer(australia_df,
                          3:6,
                          names_to = c("subfig", "var"),
                          values_to = "val",
                          names_sep = "_") %>%
  pivot_wider(names_from = var, values_from = val)

australia_df <- mutate(australia_df,
                       pval = ifelse(pval > .05, "nonsig", "sig"),
                       corr_cat = cut(ifelse(corr < 0, 0, corr), c(0, 0.20, 0.40, .60, .80, 1), include.lowest = T)) %>%
               filter(!is.na(corr))




limits <- gisco_get_countries(resolution = "10") %>%
  st_crop(xmin = 112, xmax = 160,
          ymin = -43, ymax = -9) %>% 
  st_transform("EPSG:4462")

austr <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = australia_df, aes(x, y, fill = corr_cat), alpha = .8)  +
  geom_point(data = filter(australia_df, pval == "sig"), aes(x, y, shape = "Statistical significant"), size = .8) +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  labs(fill = NULL, shape = NULL) +
  facet_wrap(~subfig, ncol = 1) +
  scale_fill_manual(values = col_def) +
  scale_shape_manual(values = 21) +
  guides(fill = guide_colorsteps(barwidth = 15, 
                                 barheight = .5,
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
        legend.text = element_text(size = 14),  # Cambia 14 por el tamaño de texto
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 20),
        plot.margin = margin(1, 10, 15, 10))



# CHILE

chile_df <- as.data.frame(rast(chile), xy = T)
names(chile_df)[3:6] <- c("b_corr", "b_pval",
                              "e_corr", "e_pval")

chile_df <- pivot_longer(chile_df,
                             3:6,
                             names_to = c("subfig", "var"),
                             values_to = "val",
                             names_sep = "_") %>%
  pivot_wider(names_from = var, values_from = val)


chile_df <- mutate(chile_df,
                   pval = ifelse(pval > .05, "nonsig", "sig"),
                   corr_cat = cut(ifelse(corr < 0, 0, corr), c(0, 0.20, 0.40, .60, .80, 1), include.lowest = T)) %>%
                 filter(!is.na(corr))


limits <- gisco_get_countries(resolution = "10") %>%
  st_crop(xmin = -80, xmax = -65,
          ymin = -56, ymax = -15) %>% 
  st_transform("EPSG:32718")

limits_inner <- ms_innerlines(limits)

chi <- ggplot() +
  geom_sf(data = limits, fill = "grey50", colour = NA) +
  geom_tile(data = chile_df, aes(x, y, fill = corr_cat), alpha = .8)  +
  geom_sf(data = limits_inner, linewidth = .3, colour = "white") +
  geom_point(data = filter(chile_df, pval == "sig"), aes(x, y, shape = "Statistical significant"), size = .8) +
  labs(fill = NULL, shape = NULL, x = NULL, y = NULL) +
  labs(fill = NULL, shape = NULL) +
  facet_wrap(~subfig, ncol = 1) +
  scale_fill_manual(values = col_def) +
  scale_shape_manual(values = 21) +
  guides(fill = guide_colorsteps(barwidth = 15, 
                                 barheight = .5,
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
        legend.text = element_text(size = 14),  # Cambia 14 por el tamaño de texto
        legend.justification = .6,
        legend.margin = margin(t = 10),
        strip.text.x = element_text(hjust = 0.1,
                                    face = "bold",
                                    size = 20),
        plot.margin = margin(5, 10, 15, 10))

(austr | eu | chi) + 
  plot_layout(guides = "collect") &
 theme(legend.position = "bottom",
        legend.box = "vertical",
        legend.box.just = "left",
        legend.spacing.y = unit(-5, "pt"),
        legend.box.margin = margin(),
        legend.title = element_markdown(),
        legend.justification = .6,
        legend.margin = margin(t = 10, b = 10))


ggsave("Figure_7_.png",
       width = 16,
       height = 10,
       units = "in",
       bg = "white",
       device = png,
       type = "cairo")



