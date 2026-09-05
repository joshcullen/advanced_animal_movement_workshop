
### Compare behavioral state estimates among models ###

library(tidyverse)
library(lubridate)
library(bayesmove)
library(sf)
library(rnaturalearth)
library(plotly)

source("R/utils.R")


#################################################
#### Load the model results from each method ####
#################################################

### SSM
ssm_1h <- read_csv("processed_data/Session_3/rw_behav_1h.csv")
ssm_8h <- read_csv("processed_data/Session_3/mp_behav_8h.csv")

### HMM
hmm_simple <- read_csv("processed_data/Session_3/HMM_3state_simple.csv")
hmm_mi <- read_csv("processed_data/Session_3/HMM_3state_MultImp.csv")
hmm_mi_covar <- read_csv("processed_data/Session_3/HMM_3state_MultImp_covar.csv")



##########################################################
#### Wrangle model results to compile into data.frame ####
##########################################################

ssm_res <- rbind(ssm_1h |> select(id, id_orig, date, lon, lat, g) |> mutate(method = "RW_1h"),
                 ssm_8h |> select(id, id_orig, date, lon, lat, g) |> mutate(method = "MP_8h"))

hmm_res <- rbind(hmm_simple |> select(ID = id_orig, date, x, y, disp, state_vit, state_fb) |> 
                   mutate(method = "Simple") |> 
                   add_trans_coords(coords = c('x','y'),
                                    proj = "+proj=utm +zone=36 +ellps=WGS84 +units=km +no_defs +south",
                                    new_proj = 4326),
                 hmm_mi |> select(ID, date, x, y, disp, state_vit, state_fb) |> 
                   mutate(method = "MI") |> 
                   add_trans_coords(coords = c('x','y'), proj = 32736, new_proj = 4326),
                 hmm_mi_covar |> select(ID, date, x, y, disp, state_vit, state_fb) |> 
                   mutate(method = "MI_covar") |> 
                   add_trans_coords(coords = c('x','y'), proj = 32736, new_proj = 4326))




############################################
#### Viz time series of state estimates ####
############################################

### SSM

# Full time series
ggplot() +
  geom_point(data = ssm_res, aes(date, g, color = method), size = 1) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  facet_wrap(~ id, scales = "free_x")

# First 30 days (i.e., zoomed in)
ggplot() +
  geom_point(data = ssm_res |> 
               group_by(id) |> 
               filter(date < (first(date) + days(30))), aes(date, g, color = method), size = 1) +
  scale_color_brewer(palette = "Set2") +
  theme_bw() +
  facet_wrap(~ id, scales = "free_x")

# Zoom-in w/ plotly
ggplotly(
  ggplot() +
    geom_point(data = ssm_res, aes(date, g, color = method), size = 1) +
    scale_color_brewer(palette = "Set2") +
    theme_bw() +
    facet_wrap(~ id, scales = "free_x")
)


### HMM

# Full time series
ggplot() +
  geom_point(data = hmm_res, aes(date, state_fb, color = method), size = 1) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(method ~ ID, scales = "free_x")

# First 30 days (i.e., zoomed in)
ggplot() +
  geom_point(data = hmm_res |> 
               group_by(ID) |> 
               filter(date < (first(date) + days(30))), aes(date, state_fb, color = method), size = 1) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(method ~ ID, scales = "free_x")

# Zoom-in w/ plotly
ggplotly(
  ggplot() +
    geom_point(data = hmm_res, aes(date, state_fb, color = method), size = 1) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    facet_grid(method ~ ID, scales = "free_x")
)







##############################################################
#### Viz map of behavioral state estimates across methods ####
##############################################################

africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")


### SSM
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = ssm_res, aes(lon, lat, group = id), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = ssm_res, aes(lon, lat, color = g), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  # scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
  scale_color_viridis_c("Move/nPersistence", option = "inferno") +
  theme_bw() +
  coord_sf(xlim = range(ssm_res$lon),
           ylim = range(ssm_res$lat)) +
  facet_grid(id_orig ~ method)


# Focus on ID 6469
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = ssm_res |> filter(id_orig == 6469), aes(lon, lat, group = id), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = ssm_res |> filter(id_orig == 6469), aes(lon, lat, color = g), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_viridis_c("Move/nPersistence", option = "inferno") +
  theme_bw() +
  coord_sf(xlim = ssm_res |> filter(id_orig == 6469) |> pull(lon) |> range(),
           ylim = ssm_res |> filter(id_orig == 6469) |> pull(lat) |> range()) +
  facet_grid(id_orig ~ method)
#both time scales seem to do a decent job at classifying behavioral patterns, just depends on time scale of interest


# Interactive focal map for RW_1h results
ggplotly(
  ggplot() +
    geom_sf(data = africa) +
    geom_path(data = ssm_res |> filter(id_orig == 6469, method == 'RW_1h'),
              aes(lon, lat, group = id), alpha = 0.5, linewidth = 0.25) +
    geom_point(data = ssm_res |> filter(id_orig == 6469, method == 'RW_1h'),
               aes(lon, lat, color = g), alpha = 0.5, size = 1) +
    geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
    scale_color_viridis_c("Move/nPersistence", option = "inferno") +
    theme_bw() +
    coord_sf(xlim = ssm_res |> filter(id_orig == 6469) |> pull(lon) |> range(),
             ylim = ssm_res |> filter(id_orig == 6469) |> pull(lat) |> range())
)




### HMM
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = hmm_res, aes(lon, lat, group = ID), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = hmm_res, aes(lon, lat, color = state_fb), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
  theme_bw() +
  coord_sf(xlim = range(hmm_res$lon),
           ylim = range(hmm_res$lat)) +
  facet_grid(ID ~ method)

# Focus on ID 6469
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = hmm_res |> filter(ID == 6469), aes(lon, lat, group = ID), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = hmm_res |> filter(ID == 6469), aes(lon, lat, color = state_fb), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
  theme_bw() +
  coord_sf(xlim = hmm_res |> filter(ID == 6469) |> pull(lon) |> range(),
           ylim = hmm_res |> filter(ID == 6469) |> pull(lat) |> range()) +
  facet_grid(method ~ ID)

# Interactive focal map
ggplotly(
  ggplot() +
    geom_sf(data = africa) +
    geom_path(data = hmm_res |> filter(ID == 6469), aes(lon, lat, group = ID), alpha = 0.5, linewidth = 0.25) +
    geom_point(data = hmm_res |> filter(ID == 6469), aes(lon, lat, color = state_fb), alpha = 0.5, size = 1) +
    geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
    scale_color_manual("State", values = c(RColorBrewer::brewer.pal(n = 3, "Dark2"), "grey")) +
    theme_bw() +
    coord_sf(xlim = hmm_res |> filter(ID == 6469) |> pull(lon) |> range(),
             ylim = hmm_res |> filter(ID == 6469) |> pull(lat) |> range()) +
    facet_grid(method ~ ID)
)



#####################################
### Explore results interactively ###
#####################################

#-- Can use Shiny app from {bayesmove} that's already in R, or can load results (as .csv, .shp, etc) into GIS software such as GeoLibre, QGIS, or ArcGIS --#

### SSM
ssm_res |> 
  rename(x = lon, y = lat) |> 
  shiny_tracks(epsg = 4326)


### HMM
hmm_res |> 
  select(-c(x,y)) |> 
  rename(x = lon, y = lat, id = ID) |> 
  shiny_tracks(epsg = 4326)



#### Main takeaways ####

#-- SSMs and HMMs both have ways of accounting for location error and time series irregularity when estimating behavioral states. SSMs estimate continuous state variable, which may allow greater flexibility for certain spp, time intervals, and movement patterns. By comparison, HMMs are convenient when practitioners have a reasonable idea what states may be present a priori, but would like a quantitative method for classifying behaviors for observations. While other packages and custom SSMs allow for inclusion of covariates on the modeling process, HMMs fit via {momentuHMM} directly allow for this both on the t.p.m., as well as the parameters of the movement metrics (e.g., step lengths, turning angles) themselves. Another nice feature of HMMs is that they provide an easier means of labeling observation given a finite number of states; however, these state names require a detailed understanding of the species and what these estimate state-dependent density distributions may be representing (which will also vary with changing time interval) --#