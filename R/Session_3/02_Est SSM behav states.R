
### Estimate behavioral states from SSM (continuous states) ###

library(aniMotum)
library(tidyverse)
library(bayesmove)
library(rnaturalearth)
library(sf)
library(tictoc)



###################
#### Load data ####
###################

### Load tracks

# Load tracks w/ bursts
dat <- read_csv("processed_data/Session_3/track_bursts.csv")

# Load regularized tracks (mean est.)
ssm_tracks <- read_csv("processed_data/Session_3/regularized_tracks.csv") |> 
  arrange(id, date) |>  #make sure data is properly sorted
  mutate(id_orig = as.character(id_orig))  #IDs are better handled as 'character'

# Load fitted SSM objects
load(file = "processed_data/Session_3/ssm_fits.RData")



### Load spatial layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")






###############################################################
### Estimate behavioral states (as move persistence; gamma) ###
###############################################################

### From average tracks

## Individual move persistence model ('mpm') estimates behavioral states separately across IDs
#Only really needed if using the "rw" or "crw" models
tic()
rw_fit_g <- fit_mpm(rw_fit, what = "predicted", model = "mpm",
                    control = mpm_control(verbose = 1))
toc()  #took 1 min to fit
#if issues w/ model not converging, try changing time step of SSM (or optimizer for fit_mpm) and re-run

print(rw_fit_g)  #all models converged
plot(rw_fit_g)
# not really clear what's happening at this fine temporal scale


# Zoom in on time series to explore diel patterns
tmp <- grab(rw_fit_g, normalise = TRUE) |> 
  mutate(id_orig = as.vector(str_match(id, "[0-9]+")),
         .after = id)

plotly::ggplotly(
  ggplot(tmp |> 
           filter(id_orig == 6469)) +
    geom_point(aes(date, g, color = g)) +
    scale_color_viridis_c(expression(gamma), option = "plasma") +
    theme_bw(base_size = 14)
)



## Try fitting SSM at coarser time scale and re-evaluating

# Remove bursts w/ few points (e.g., n < 30)
dat2 <- dat |>
  group_by(id, burst_id) |>
  filter(n() >= 30) |>
  group_by(id) |>
  mutate(burst_id = dense_rank(burst_id)) |>  #ensure all bursts are consecutive and begin at 1
  ungroup() |>
  mutate(burst_id = paste(id, burst_id, sep = "_"))  #create unique burst IDs for model fitting

dat_sf <- dat2 |> 
  rename(id_orig = id, id = burst_id) |>  #treat 'burst_id' as primary ID for model fitting
  relocate(lon, .before = lat) |>  #fix column order for aniMotum::fit_ssm()
  mutate(lc = 'G', .after = date) |>  #need to specify "location class" (G = GPS)
  st_as_sf(coords = c('lon','lat'), crs = 4326, remove = FALSE) |>  #convert to spatial object
  st_transform(32736)  #convert to a projected ref system (e.g., UTM)

tic()
rw_fit2 <- fit_ssm(dat_sf, 
                   model = "rw", 
                   time.step = 8,  #8 hr
                   spdf = FALSE,
                   control = ssm_control(verbose = 1, tdist = "norm"))
toc()  #took 1.5 min to fit


print(rw_fit2)  #all indiv. models converged
summary(rw_fit2)

# Viz the fitted vs observed pts; blue = observed; gold = fitted
plot(rw_fit2, what = "predicted", type = 1, ask = TRUE)  #plot time series of coords
plot(rw_fit2, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)  #plot maps of tracks

# Estimate move persistence at 8 hr scale
tic()
rw_fit_g8h <- fit_mpm(rw_fit2, what = "predicted", model = "mpm",
                      control = mpm_control(verbose = 1))
toc()  #took 9 sec to fit

print(rw_fit_g8h)  #all models converged
plot(rw_fit_g8h)





### Using previously fitted 'move persistence' model
plot(mp_fit, type = 3, normalise = TRUE)
plot(mp_fit, type = 4, normalise = TRUE)
# more variability than other 2-step approach

# Try re-fitting at coarser time scale
tic()
mp_fit2 <- fit_ssm(dat_sf, 
                   model = "mp", 
                   time.step = 8,  #8 hr
                   spdf = FALSE,
                   control = ssm_control(verbose = 1, tdist = "norm"))
toc()  #took 2 min to fit


print(mp_fit2)  #all indiv. models converged
summary(mp_fit2)

# Viz the fitted vs observed pts; blue = observed; gold = fitted
plot(mp_fit2, what = "predicted", type = 1, ask = TRUE)  #plot time series of coords
plot(mp_fit2, what = "predicted", type = 2, alpha = 0.1, ask = TRUE)  #plot maps of tracks

plot(mp_fit2, type = 3, normalise = TRUE)
plot(mp_fit2, type = 4, normalise = TRUE)
#behavioral patterns match up much better now at this coarser scale





###################################################
### Extract annotated tracks from model objects ###
###################################################

# Via 2-stage method
rw_behav <- join(ssm = rw_fit,
                 mpm = rw_fit_g,
                 what.ssm = "predicted",
                 normalise = TRUE) |> 
  mutate(id_orig = as.vector(str_match(id, "[0-9]+")),
         .after = id)

# Via single stage 'mp' model
mp_behav <- grab(mp_fit2, what = "predicted", normalise = TRUE) |> 
  mutate(id_orig = as.vector(str_match(id, "[0-9]+")),
         .after = id)



### Viz mapped tracks

# RW model
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = rw_behav, aes(lon, lat, group = id), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = rw_behav, aes(lon, lat, color = g), alpha = 0.5, size = 0.5) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_viridis_c("Move\nPersistence", option = "viridis", limits = c(0,1)) +
  theme_bw() +
  coord_sf(xlim = range(dat2$lon),
           ylim = range(dat2$lat))

# MP model
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = mp_behav, aes(lon, lat, group = id), alpha = 0.5, linewidth = 0.25) +
  geom_point(data = mp_behav, aes(lon, lat, color = g), alpha = 0.5, size = 1) +
  geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
  scale_color_viridis_c("Move\nPersistence", option = "viridis", limits = c(0,1)) +
  theme_bw() +
  coord_sf(xlim = range(dat2$lon),
           ylim = range(dat2$lat))


# Interactive mapping
bayesmove::shiny_tracks(data = rw_behav,
                        epsg = "+proj=utm +zone=36 +ellps=WGS84 +units=km +no_defs +south")


#-- Not explored in detail here, but time scale has large impact on ecological inferences. Care should be taken to address research questions when selecting a time scale and behavioral state method. For example, the behavioral state estimates at the 1 hr time step did not appear to provide informative behavioral states from this model, but DID perform better for coarser behavioral patterns. However, custom Bayesian versions of this SSM may produce more informative results, such as through the inclusion of covariates. --#




###############################
### Export annotated tracks ###
###############################

write_csv(rw_behav, "processed_data/Session_3/rw_behav_1h.csv")
write_csv(mp_behav, "processed_data/Session_3/mp_behav_8h.csv")
