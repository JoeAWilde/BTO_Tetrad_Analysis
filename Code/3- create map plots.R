library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(ggpubr)
library(progress)
library(ggspatial)

counties <- st_read("Data/county_outlines/counties_shape/CTYUA_MAY_2023_UK_BGC.shp") 

#Read in the data ####
df <- read.table("Outputs/script_1/ModelData.csv", sep = ",", header = T, fill = T) %>%
  na.omit() %>%
  pivot_longer(cols = c(pheasant_breed, pheasant_wint), names_to = "time_period", values_to = "abund") %>%
  filter(abund != "M") %>%
  mutate(rel_abund = 0)

#Convert the categorical abundance measures to continuous ####
for(i in 1:nrow(df)) {
  if(i == 1) {
    pb <- progress_bar$new(total = nrow(df), format = "[:bar] :percent eta::eta", clear = F)
    pb$tick(0)
  }
  if(df$County[i] == "Cornwall") {
    df$rel_abund[i] <- case_when(df$abund[i] == "0" ~ 0,
                                 df$abund[i] == "0.11-1" ~ 0.55,
                                 df$abund[i] == "0.2-1" ~ 0.6,
                                 df$abund[i] == "1" ~ 1,
                                 df$abund[i] == "1.33-1.5" ~ 1.415,
                                 df$abund[i] == "1.25-2" ~ 1.625,
                                 df$abund[i] == "1.6-3" ~ 2.3,
                                 df$abund[i] == "2" ~ 2,
                                 df$abund[i] == "2.25-3" ~ 2.625,
                                 df$abund[i] == "3" ~ 3,
                                 df$abund[i] == "3.25-5" ~ 4.125,
                                 df$abund[i] == "3.5-6" ~ 4.75,
                                 df$abund[i] == "4" ~ 4,
                                 df$abund[i] == "5-32" ~ 18.5,
                                 df$abund[i] == "5.33-35" ~ 20.165, 
                                 df$abund[i] == "7-85" ~ 46)
  } else if(df$County[i] == "Devon") {
    df$rel_abund[i] <- case_when(df$abund[i] == "0" ~ 0,
                                 df$abund[i] == "1" ~ 0.5,
                                 df$abund[i] == "2" ~ 1.5,
                                 df$abund[i] == "3" ~ 3.5,
                                 df$abund[i] == "4" ~ 7.5,
                                 df$abund[i] == "5" ~ 15,
                                 df$abund[i] == "6" ~ 35,
                                 df$abund[i] == "7" ~ 50)
  } else if(df$County[i] == "Berkshire") {
    df$rel_abund[i] <- case_when(df$abund[i] == "0" ~ (0*7.3),
                                 df$abund[i] == "1" ~ (0.1*7.3),
                                 df$abund[i] == "2" ~ (0.25*7.3),
                                 df$abund[i] == "3" ~ (0.75*7.3),
                                 df$abund[i] == "4" ~ (1.25*7.3),
                                 df$abund[i] == "5" ~ (1.75*7.3),
                                 df$abund[i] == "6" ~ (2.25*7.3),
                                 df$abund[i] == "7" ~ (2.5*7.3))
  } else {
    df$rel_abund[i] <- as.numeric(df$abund[i])
  }
  
  pb$tick()
}

#Tidy the dataframe to be used in the analysis ####
post_df <- df %>%
  mutate(
    tet = as.numeric(factor(X, levels = unique(X))),
    time_period = factor(time_period, levels = c("pheasant_wint", "pheasant_breed"))) %>%
  group_by(County, time_period) %>% 
  mutate(
    N_Birds = sum(rel_abund), 
    final_relabund = rel_abund / N_Birds
    ) %>%
  ungroup()

PAs <- st_read("Outputs/script_1/Suitable_Tetrads.shp")

#Devon plot - winter
dev_uas <- c("Devon", "Plymouth", "Torbay")

dev_shp <- counties %>%
    filter(CTYUA23NM %in% dev_uas) %>%
    summarise(geometry = st_union(geometry))

devwint_df <- post_df %>%
    filter(County == "Devon" & time_period == "pheasant_wint")

dev_PA <- PAs %>%
    st_filter(dev_shp, .pred = st_intersects)

devwint_plot <- ggplot() + 
    geom_point(data = devwint_df, aes(x = X, y = Y, colour = final_relabund), size = 1.5) + 
    geom_sf(data = dev_shp, fill = "transparent", linewidth = 0.2, colour = "black") + 
    geom_sf(data = dev_PA, fill = "grey", alpha = 0.6) +
    scale_colour_gradientn(name = "Relative abundance", 
        colours = c("low" = "purple4", "high" = "orange2")) + 
    annotation_north_arrow(
which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.1, "cm"), 
        text_cex = 1) + 
    theme_bw(base_size = 18) + 
    theme(legend.key.width = unit(0.3, "cm"), 
        axis.title = element_blank())

# ggsave(dev_plot, path = "Outputs/plots/MS_final_plots/", filename = "dev_wint_plot.png", 
#     units = "px", height = 1280, width = 1280)

#Cornwall plot - winter
cw_shp <- counties %>%
    filter(CTYUA23NM == "Cornwall")

cwwint_df <- post_df %>%
    filter(County == "Cornwall" & time_period == "pheasant_wint")

cw_PA <- PAs %>%
    st_filter(cw_shp, .pred = st_intersects)

cwwint_plot <- ggplot() + 
    geom_point(data = cwwint_df, aes(x = X, y = Y, colour = final_relabund), size = 1.5) + 
    geom_sf(data = cw_shp, fill = "transparent", linewidth = 0.2, colour = "black") + 
    geom_sf(data = cw_PA, fill = "grey", alpha = 0.6) +
    scale_colour_gradientn(name = "Relative abundance", 
        colours = c("low" = "purple4", "high" = "orange2")) + 
    annotation_north_arrow(
which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.1, "cm"), 
        text_cex = 1) + 
    theme_bw(base_size = 18) + 
    theme(legend.key.width = unit(0.3, "cm"), 
        axis.title = element_blank())

# ggsave(cw_plot, path = "Outputs/plots/MS_final_plots/", filename = "cw_wint_plot.png", 
#     units = "px", height = 1280, width = 1280)


#Berkshire plot - winter
bk_uas <- c("Bracknell Forest", "Reading", "Slough", "West Berkshire", "Windsor and Maidenhead", "Wokingham")

bk_shp <- counties %>%
    filter(CTYUA23NM %in% bk_uas) %>%
    summarise(geometry = st_union(geometry))

bkwint_df <- post_df %>%
    filter(County == "Berkshire" & time_period == "pheasant_wint")

bk_PA <- PAs %>%
    st_filter(bk_shp, .pred = st_intersects)

bkwint_plot <- ggplot() + 
    geom_point(data = bkwint_df, aes(x = X, y = Y, colour = final_relabund), size = 1.9) + 
    geom_sf(data = bk_shp, fill = "transparent", linewidth = 0.2, colour = "black") + 
    geom_sf(data = bk_PA, fill = "grey", alpha = 0.6) +
    scale_colour_gradientn(name = "Relative abundance", 
        colours = c("low" = "purple4", "high" = "orange2")) + 
    annotation_north_arrow(
which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.05, "cm"), 
        text_cex = 1) + 
    theme_bw(base_size = 18) + 
    theme(legend.key.width = unit(0.3, "cm"), 
        axis.title = element_blank())

# ggsave(bk_plot, path = "Outputs/plots/MS_final_plots/", filename = "bk_wint_plot.png", 
#     units = "px", height = 680, width = 1280)


#Hertfordshire plot - winter
hf_shp <- counties %>%
    filter(CTYUA23NM == "Hertfordshire")

hfwint_df <- post_df %>%
    filter(County == "Hertfordshire" & time_period == "pheasant_wint")

hf_PA <- PAs %>%
    st_filter(hf_shp, .pred = st_intersects)

hfwint_plot <- ggplot() + 
    geom_point(data = hfwint_df, aes(x = X, y = Y, colour = final_relabund), size = 2.3) + 
    geom_sf(data = hf_shp, fill = "transparent", linewidth = 0.2, colour = "black") + 
    geom_sf(data = hf_PA, fill = "grey", alpha = 0.6) + 
    scale_colour_gradientn(name = "Relative abundance", 
        colours = c("low" = "purple4", "high" = "orange2")) + 
    annotation_north_arrow(
which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.05, "cm"), 
        text_cex = 1) + 
    theme_bw(base_size = 18) + 
    theme(legend.key.width = unit(0.3, "cm"), 
        axis.title = element_blank())

# ggsave(hf_plot, path = "Outputs/plots/MS_final_plots/", filename = "hf_wint_plot.png", 
#     units = "px", height = 1280, width = 1280)

# all_wint <- ggarrange(bk_plot, cw_plot, dev_plot, hf_plot, 
#     ncol = 2, 
#     nrow = 2, 
#     labels = c("a) Berkshire", "b) Cornwall", "c) Devon", "d) Hertfordshire"), 
#     common.legend = T, 
#     legend = "right"
# )

# ggsave(all_wint, path = "Outputs/plots/MS_final_plots", filename = "all_wint.png", 
#     units = "px", height = 1280 * 4, width = 1280 * 4)






#Devon plot - breeding
dev_uas <- c("Devon", "Plymouth", "Torbay")

dev_shp <- counties %>%
    filter(CTYUA23NM %in% dev_uas) %>%
    summarise(geometry = st_union(geometry))

devbreed_df <- post_df %>%
    filter(County == "Devon" & time_period == "pheasant_breed")

dev_PA <- PAs %>%
    st_filter(dev_shp, .pred = st_intersects)

devbreed_plot <- ggplot() + 
    geom_point(data = devbreed_df, aes(x = X, y = Y, colour = final_relabund), size = 1.5) + 
    geom_sf(data = dev_shp, fill = "transparent", linewidth = 0.2, colour = "black") + 
    geom_sf(data = dev_PA, fill = "grey", alpha = 0.6) +
    scale_colour_gradientn(name = "Relative abundance", 
        colours = c("low" = "purple4", "high" = "orange2")) + 
    annotation_north_arrow(
which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.1, "cm"), 
        text_cex = 1) + 
    theme_bw(base_size = 18) + 
    theme(
        legend.key.width = unit(0.3, "cm"), 
        axis.title = element_blank()
    )

# ggsave(dev_plot, path = "Outputs/plots/MS_final_plots/", filename = "dev_breed_plot.png", 
#     units = "px", height = 1280, width = 1280)

#Cornwall plot - breeding
cw_shp <- counties %>%
    filter(CTYUA23NM == "Cornwall")

cwbreed_df <- post_df %>%
    filter(County == "Cornwall" & time_period == "pheasant_breed")

cw_PA <- PAs %>%
    st_filter(cw_shp, .pred = st_intersects)

cwbreed_plot <- ggplot() + 
    geom_point(data = cwbreed_df, aes(x = X, y = Y, colour = final_relabund), size = 1.5) + 
    geom_sf(data = cw_shp, fill = "transparent", linewidth = 0.2, colour = "black") + 
    geom_sf(data = cw_PA, fill = "grey", alpha = 0.6) +
    scale_colour_gradientn(name = "Relative abundance", 
        colours = c("low" = "purple4", "high" = "orange2")) + 
    annotation_north_arrow(
which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.1, "cm"), 
        text_cex = 1) + 
    theme_bw(base_size = 18) + 
    theme(legend.key.width = unit(0.3, "cm"), 
        axis.title = element_blank())

# ggsave(cw_plot, path = "Outputs/plots/MS_final_plots/", filename = "cw_breed_plot.png", 
#     units = "px", height = 1280, width = 1280)


#Berkshire plot - breeding
bk_uas <- c("Bracknell Forest", "Reading", "Slough", "West Berkshire", "Windsor and Maidenhead", "Wokingham")

bk_shp <- counties %>%
    filter(CTYUA23NM %in% bk_uas) %>%
    summarise(geometry = st_union(geometry))

bkbreed_df <- post_df %>%
    filter(County == "Berkshire" & time_period == "pheasant_breed")

bk_PA <- PAs %>%
    st_filter(bk_shp, .pred = st_intersects)

bkbreed_plot <- ggplot() + 
    geom_point(data = bkbreed_df, aes(x = X, y = Y, colour = final_relabund), size = 1.9) + 
    geom_sf(data = bk_shp, fill = "transparent", linewidth = 0.2, colour = "black") + 
    geom_sf(data = bk_PA, fill = "grey", alpha = 0.6) +
    scale_colour_gradientn(name = "Relative abundance", 
        colours = c("low" = "purple4", "high" = "orange2")) + 
    annotation_north_arrow(
which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.05, "cm"), 
        text_cex = 1) + 
    theme_bw(base_size = 18) + 
    theme(legend.key.width = unit(0.3, "cm"), 
        axis.title = element_blank())

# ggsave(bk_plot, path = "Outputs/plots/MS_final_plots/", filename = "bk_breed_plot.png", 
#     units = "px", height = 680, width = 1280)


#Hertfordshire plot - breeding
hf_shp <- counties %>%
    filter(CTYUA23NM == "Hertfordshire")

hfbreed_df <- post_df %>%
    filter(County == "Hertfordshire" & time_period == "pheasant_breed")

hf_PA <- PAs %>%
    st_filter(hf_shp, .pred = st_intersects)

hfbreed_plot <- ggplot() + 
    geom_point(data = hfbreed_df, aes(x = X, y = Y, colour = final_relabund), size = 2.3) + 
    geom_sf(data = hf_shp, fill = "transparent", linewidth = 0.2, colour = "black") + 
    geom_sf(data = hf_PA, fill = "grey", alpha = 0.6) + 
    scale_colour_gradientn(name = "Relative abundance", 
        colours = c("low" = "purple4", "high" = "orange2")) + 
    annotation_north_arrow(
        which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.05, "cm"), 
        text_cex = 1) + 
    theme_bw(base_size = 18) + 
    theme(legend.key.width = unit(0.3, "cm"), 
        axis.title = element_blank())

# ggsave(hf_plot, path = "Outputs/plots/MS_final_plots/", filename = "hf_breed_plot.png", 
#     units = "px", height = 1280, width = 1280)

all <- ggarrange(bkwint_plot, bkbreed_plot, cwwint_plot, cwbreed_plot, devwint_plot, devbreed_plot, hfwint_plot, hfbreed_plot,
    ncol = 2, 
    nrow = 4, 
    labels = c("a) Berkshire - winter", "b) Berkshire - breeding",
                "c) Cornwall - winter", "d) Cornwall - breeding", 
                "e) Devon - winter", "f) Devon - breeding", 
                "g) Hertfordshire - winter", "h) Hertfordshire - breeding"), 
    common.legend = T, 
    legend = "right"
)

ggsave(all, path = "Outputs/plots/MS_final_plots", filename = "all_plots.png", 
    units = "px", height = 1280 * 8, width = 1280 * 4)

UK_outline <- st_read("Data/CoastOutline/UK_Coastline.shp")

uk_plot <- ggplot() + 
    geom_sf(data = UK_outline, fill = "transparent") + 
    geom_sf(data = dev_shp, fill = "darkgrey") + 
    geom_text(aes(x = -3.8, y = 50.75, label = "Devon"), size = 3) + 
    geom_sf(data = cw_shp, fill = "darkgrey") + 
    geom_text(aes(x = -4.9, y = 50.38, label = "Cornwall"), size = 2.6) + 
    geom_sf(data = bk_shp, fill = "darkgrey") + 
    geom_text(aes(x = -1.072329, y = 51.435, label = "Berkshire"), size = 3) + 
    geom_sf(data = hf_shp, fill = "darkgrey") + 
    geom_text(aes(x = -0.23, y = 51.8088, label = "Hertfordshire"), size = 2.6) + 
    annotation_north_arrow(
        which_north = "true",
        location = "tr", 
        height = unit(0.3, "cm"), 
        width = unit(0.3, "cm"), 
        style = north_arrow_orienteering(text_size = 3)
    ) + 
    annotation_scale(
        height = unit(0.3, "cm"), 
        pad_x = unit(0.2, "cm"),
        pad_y = unit(0.05, "cm"), 
        text_cex = 1) +
    theme_bw() + 
    theme(axis.title = element_blank()) + 
    coord_sf(ylim = c(49.99, 52.05), xlim = c(-6, 0.05))

ggsave(uk_plot, path = "Outputs/plots/MS_final_plots/", filename = "county_map.png", 
    units = "px", height = 1080, width = 1920)
