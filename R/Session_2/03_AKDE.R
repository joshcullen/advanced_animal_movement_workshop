
### Calculate autocorrelated kernel density estimates (AKDE) ###

library(tidyverse)
library(ctmm)
library(sf)
library(rnaturalearth)
library(tictoc)
library(terra)
library(progressr)
library(furrr)
library(plotly)

source("R/utils.R")  #load in custom functions


###################
#### Load data ####
###################

dat <- read_csv('processed_data/Session_1/cleaned_tracks.csv')

glimpse(dat)
summary(dat)


# Load spatial layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")





###############################
#### Wrangle and prep data ####
###############################

# Interactively explore whether any hard boundaries blocking movement
plotly::ggplotly(
  ggplot() +
    geom_sf(data = africa) +
    geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), alpha = 0.5, linewidth = 0.25) +
    geom_sf(data = gl_pa, color = "black", fill = NA, linewidth = 0.25) +
    scale_color_brewer("ID", palette = "Dark2") +
    theme_bw() +
    coord_sf(xlim = range(dat$lon),
             ylim = range(dat$lat))
)
# Some hard boundaries limiting, but not fully blocking, movements of IDs 5605 and 6469


# Explore movement pattern (i.e., resident vs migratory) w/ net displacement
dat <- dat |> 
  add_trans_coords(coords = c('lon','lat'), proj = 4326, new_proj = 32736) |> 
  group_by(id) |> 
  arrange(date, .by_group = TRUE) |> 
  mutate(disp = sqrt((x - x[1])^2 + (y - y[1])^2),  #calc distance from first location
         disp = disp / 1000) |>  #convert units from m to km
  ungroup()

ggplot(dat) +
  geom_path(aes(date, disp, color = factor(id)), linewidth = 0.5) +
  scale_color_brewer("ID", palette = "Dark2") +
  theme_bw() +
  facet_wrap(~ id, scales = "free", ncol = 2)
#you can make the case that each ID makes one or more migrations away from primary home range
#however, these ranging bouts are relatively short-lived (~1-3 months); good to also check variograms to confirm

# Let's also do quick check w/ interactive map
bayesmove::shiny_tracks(dat, epsg = 32736)



# Switch to Movebank naming conventions
dat2 <- dat |> 
  rename(tag.local.identifier = id,
         timestamp = date,
         location.long = lon,
         location.lat = lat) |> 
  mutate(GPS.Horizontal.Error = 30)  #add column to store presumed location error (30 m)

# Convert to 'telemetry' class
dat.telem <- as.telemetry(dat2, projection = "EPSG:32736")

dat.telem
summary(dat.telem)
plot(dat.telem, col = rainbow(length(dat.telem)))
#in this case, there seems to be some dispersing/migratory behavior



### Explore time intervals

# For one ID
dt.plot(dat.telem$`6471`)
#mostly at 1 hr 

# For all tracks
dt.plot(dat.telem)
#mostly at 1 hr, but also small number at 2 and 4 hrs


# 1, 2, 4 hour sampling intervals
dt <- c(1,2,4) %#% "hour"





#########################
#### Plot variograms ####
#########################

# Generate variogram per track (using primary observed sampling intervals)
svf <- map(dat.telem, variogram, dt = dt)
level <- c(0.5, 0.95) # 50% and 95% CIs
xlim <- c(0,12) %#% "hour" # 0-12 hour window

### Explore individually
plot(svf$`5605`, xlim = xlim, level = level)
title("zoomed in")

plot(svf$`5605`, fraction = 0.65, level = level)
title("zoomed out")



### Explore together
par(mfrow = c(1,3))

# Zoomed-in
map(svf, \(x) plot(x, xlim = xlim, level = level))

# Zoomed-out
map(svf, \(x) plot(x, fraction = 0.65, level = level))

par(mfrow = c(1,1))
#all variograms look like there's some shift happening at 4-10 months
#if dispersal/migratory movements are present, tracks need to be segmented to only analyze the range-resident components






##################################################
### Segment tracks into range-resident periods ###
##################################################

#-- There's a number of options for doing this quantitatively, such as through segmentation models in the {segclust2d}, {bayesmove}, and {mcp} packages, but this doesn't always return great results. Alternatively, you could fit a SSM and then use coarse time steps to roughly delineate when these migratory periods occur. Otherwise, doing this manually may be necessary --#

#Viz displacement time series again
ggplot(dat) +
  geom_path(aes(date, disp, color = factor(id)), linewidth = 0.5) +
  scale_color_brewer("ID", palette = "Dark2") +
  theme_bw() +
  facet_wrap(~ id, scales = "free", ncol = 2)


# Let's manually identify start/stop times of resident and migratory periods
bayesmove::shiny_tracks(dat, epsg = 32736)

burst_windows <- tibble(id = rep(c(5605, 6469, 6471), c(9,12,12)),
                        burst = c(1:9, 1:12, 1:12),
                        start = c("2022-03-17 02:53:00", "2022-03-30 18:00:00", "2022-04-01 19:00:00", 
                                  "2022-04-17 21:00:00", "2022-04-22 12:00:00", "2022-12-08 04:00:00",
                                  "2022-12-09 19:00:00", "2023-01-09 02:00:00", "2023-01-13 14:00:00",
                                  
                                  "2022-11-23 16:21:00", "2022-12-25 14:54:00", "2023-02-22 15:02:00",
                                  "2023-04-02 18:52:00", "2023-11-29 23:53:00", "2023-12-14 23:45:00",
                                  "2024-02-19 03:57:00", "2024-03-03 19:09:00", "2024-11-29 16:52:00",
                                  "2024-12-27 05:15:00", "2025-03-24 22:30:00", "2025-04-09 18:50:00",
                                  
                                  "2022-11-23 16:12:00", "2022-11-26 16:27:00", "2023-01-11 12:18:00",
                                  "2023-01-12 21:28:00", "2023-01-14 11:39:00", "2023-01-24 01:49:00",
                                  "2023-04-14 10:36:00", "2023-04-16 02:50:00", "2023-06-18 10:30:00",
                                  "2023-06-20 02:43:00", "2025-03-07 06:43:00", "2025-03-21 21:28:00"),
                        
                        end = c("2022-03-30 18:00:00", "2022-04-01 19:00:00", "2022-04-17 21:00:00",
                                "2022-04-22 12:00:00", "2022-12-08 04:00:00", "2022-12-09 19:00:00",
                                "2023-01-09 02:00:00", "2023-01-13 14:00:00", "2023-03-28 08:00:00",
                                
                                "2022-12-25 14:54:00", "2023-02-22 15:02:00", "2023-04-02 18:52:00",
                                "2023-11-29 23:53:00", "2023-12-14 23:45:00", "2024-02-19 03:57:00",
                                "2024-03-03 19:09:00", "2024-11-29 16:52:00", "2024-12-27 05:15:00",
                                "2025-03-24 22:30:00", "2025-04-09 18:50:00", "2025-07-04 05:09:00",
                                
                                "2022-11-26 16:27:00", "2023-01-11 12:18:00", "2023-01-12 21:28:00",
                                "2023-01-14 11:39:00", "2023-01-24 01:49:00", "2023-04-14 10:36:00",
                                "2023-04-16 02:50:00", "2023-06-18 10:30:00", "2023-06-20 02:43:00",
                                "2025-03-07 06:43:00", "2025-03-21 21:28:00", "2025-07-04 05:01:00")
                        )


# Format df
burst_windows2 <- burst_windows |> 
  mutate(across(start:end, as_datetime))

# Calc variance in 'disp' and use to assign as "resident" or "migratory"
burst_windows3 <- dat |> 
  inner_join(burst_windows2,
             by = join_by(id, between(date, start, end))) |> 
  summarize(.by = c(id, burst, start, end),
            disp_var = var(disp)) |> 
  # Use threshold for IDs 5605 and 6469
  mutate(phase = case_when(id %in% c(5605, 6469) & disp_var < 100 ~ "resident",
                           id %in% c(5605, 6469) & disp_var >= 100 ~ "migratory",
                           TRUE ~ NA)) |> 
  # Manually assign give known order for ID 6471
  mutate(phase = case_when(id == 6471 & (burst %% 2) == 0 ~ "resident",
                           id == 6471 & (burst %% 2) != 0 ~ "migratory",
                           TRUE ~ phase))


# Add phase classes to data and remove "migratory" segments
dat3 <- dat |> 
  inner_join(burst_windows3,
            by = join_by(id, between(date, start, end))) |> 
  filter(phase == 'resident') |> 
  mutate(burst_id = paste(id, burst, sep = "_"),
         .after = id)



### Prep segmented tracks for analysis ###

# Switch to Movebank naming conventions
dat4 <- dat3 |> 
  select(-id) |> 
  rename(tag.local.identifier = burst_id,  #use the burst_id
         timestamp = date,
         location.long = lon,
         location.lat = lat) |> 
  mutate(GPS.Horizontal.Error = 30)  #add column to store presumed location error (30 m)

# Convert to 'telemetry' class
dat.telem2 <- as.telemetry(dat4, projection = "EPSG:32736")

dat.telem2
summary(dat.telem2)
plot(dat.telem2, col = rainbow(length(dat.telem2)))
#some forays present, but much better


# For all tracks
dt.plot(dat.telem2)
#mostly at 1 hr, but also small number at 2 and 4 hrs


# 1, 2, 4 hour sampling intervals
dt <- c(1,2,4) %#% "hour"



# Generate variogram per track (using primary observed sampling intervals)
svf2 <- map(dat.telem2, variogram, dt = dt)
level <- c(0.5, 0.95) # 50% and 95% CIs
xlim <- c(0,12) %#% "hour" # 0-12 hour window

### Explore individually
plot(svf2$`5605_1`, xlim = xlim, level = level)
title("zoomed in")

plot(svf2$`5605_1`, fraction = 0.65, level = level)
title("zoomed out")



### Explore together
par(mfrow = c(2,3), ask = TRUE)

# Zoomed-out
map(svf2, \(x) plot(x, fraction = 0.65, level = level))

par(mfrow = c(1,1), ask = FALSE)
#in general, variograms looking much better
#some segments may be a little too short, but these can be discarded later if necessary







##################################################
### Fit continuous-time movement models (CTMM) ###
##################################################

# Guess at the variogram and model to be used per track
ctmm_guess <- map2(.x = dat.telem2,
                   .y = svf2,
                   .f = ~ctmm.guess(data = .x, variogram = .y, interactive = FALSE)
                   )
# map(ctmm_guess, plot)


### Fit multiple CTMMs and select the best-fitting model

# Example for 1 track/burst
tic()
ctmm_5605_3 <- ctmm.select(data = dat.telem2$`5605_3`, CTMM = ctmm_guess$`5605_3`, IC = 'AICc',
                         verbose = TRUE, trace = 1, cores = 10, method = 'pHREML')
toc()  #took 13 sec

summary(ctmm_5605_3)
#shows model comparison by AICc (when `verbose = TRUE`), otherwise returns best model ONLY




# Example for multiple tracks/bursts
plan(multisession, workers = 10)  #run segments in parallel
handlers(handler_progress(format = ":spin :current/:total [:bar] :percent in :elapsed",
                                                width = 100,
                                                clear = FALSE))


with_progress({
  
  p <- progressor(along = dat.telem2)  #create progress bar
  
  ctmm_fit <- future_map2(.x = dat.telem2,
                          .y = ctmm_guess,
                          .f = ~{
                            tmp <- ctmm.select(data = .x, CTMM = .y, verbose = TRUE,
                                               method = 'pHREML', cores = 4)
                            p()  #print progress
                            tmp  #print results
                          },
                          .options = furrr_options(seed = 2026)
  )
  
})

plan(sequential)  #return to single core
# takes 6 min to run w/ 10 cores

map(ctmm_fit, summary)
#OUF anisotropic model typically best for all track bursts (of the 4 models compared)
#DOF (effective degrees of freedom) generally reflects number of range crossings; this should ideally be >3 for use in AKDE
##In this case, I'll just drop any segments where DOF[area] < 3 for simplicity
#To compare only 'dispersive' models, specify `CTMM = ctmm(range = FALSE)` within ctmm::ctmm.guess()
#Can use ctmm.fit() to simply fit an IID model (essentially the same as KDE w/ href bandwidth estimator)




# Which segments have large enough DOF?
idx <- map(ctmm_fit, pluck, 1) |>  #"plucks" first list element per track
  map(summary) |> 
  map(~{.x$DOF["area"]}) |>  #pull out DOF value
  unlist() |> 
  (\(z) which(z > 3))()
  

# Pull the best fitting models per track (listed first)
ctmm_fit_best <- map(ctmm_fit[idx], pluck, 1)  #"plucks" first list element per track



### Explore movement metrics from fitted models

# Home range area ('area')
# Range-crossing time (tau_position)
# Directional persistence timescale (tau_velocity)
# Average speed
# Expected square displacement over given time period (diffusion)

map(ctmm_fit_best,
    ~{.x |> 
        summary() |> 
        pluck("CI")})






#######################################
### Generate and viz AKDE per track ###
#######################################

### Fit standard AKDE model for 1 burst
tic()
akde_5605_1 <- akde(data = dat.telem2$`5605_1`, CTMM = ctmm_fit_best$`5605_1`,
             grid = list(dr = 1000, align.to.origin = TRUE, dr.fn = max))
toc()  #took 0.5 sec
#option 'dr' for arg `grid` specifies output spatial resolution (in meters)

summary(akde_5605_1)
plot(dat.telem2$`5605_1`, UD = akde_5605_1)



### Fit weighted AKDE model for 1 track (to account for sampling bias from irregular time series)
dt1 <- 0.5 %#% "hour"  #setting to half of median time step

tic()
akde_5605_1w <- akde(data = dat.telem2$`5605_1`, CTMM = ctmm_fit_best$`5605_1`, weights = TRUE, dt = dt1,
                   grid = list(dr = 1000, align.to.origin = TRUE, dr.fn = max))
toc()  #took 1 sec

summary(akde_5605_1w)
plot(dat.telem2$`5605_1`, UD = akde_5605_1w)
#very similar results



### Fit weighted AKDE across all tracks
tic()
akde <- akde(data = dat.telem2[idx], CTMM = ctmm_fit_best, weights = TRUE, dt = dt1,
             grid = list(dr = 1000, align.to.origin = TRUE, dr.fn = max)
             )
toc()  #took 18 min

map(akde, summary)
plot(dat.telem2[idx], UD = akde, col = rainbow(length(dat.telem2[idx])))




### Extract 50 and 95% contours (core area, home range)
akde_sf <- map(akde,
               ~as.sf(.x, level.UD = c(0.5, 0.95))) |> 
  bind_rows() |> 
  filter(str_detect(name, "est")) |>  #keep only the mean prediction
  separate_wider_delim(cols = name, delim = " ", names = c("id", "level", NA)) |>  #split messy 'name' column
  mutate(id_orig = as.vector(str_match(id, "[0-9]+")),
         .after = id) |>  #add column for animal ID
  st_sf(crs = 'epsg:32736')




### Viz map of estimates

# Plot separately by ID
ggplot() +
  geom_sf(data = africa |> 
            st_transform(crs = 32736) |> 
            dplyr::select(-level)) +
  geom_sf(data = akde_sf, aes(color = level), fill = NA, linewidth = 0.5) +
  scale_color_brewer(palette = 'Set1') +
  theme_bw() +
  coord_sf(xlim = st_bbox(akde_sf)[c(1,3)],
           ylim = st_bbox(akde_sf)[c(2,4)]) +
  facet_wrap(~id)


# Plot all tracks per UD level
ggplot() +
  geom_sf(data = africa |> 
            st_transform(crs = 32736) |> 
            dplyr::select(-level)) +
  geom_path(data = dat3 |> 
              rename(id_orig = id, id = burst_id),  #match colnames to akde_sf
            aes(x, y, color = factor(id_orig), group = id), linewidth = 0.25, alpha = 0.5) +
  geom_sf(data = akde_sf, aes(color = factor(id_orig)), fill = NA, linewidth = 1) +
  scale_color_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = st_bbox(akde_sf)[c(1,3)],
           ylim = st_bbox(akde_sf)[c(2,4)]) +
  facet_wrap(~ level, ncol = 2)


# Plotted all tracks and UD levels together
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = dat3 |> 
              rename(id_orig = id, id = burst_id), aes(lon, lat, group = id, color = factor(id_orig)),
            alpha = 0.5, linewidth = 0.25) +
  # Plot 95% isopleths
  geom_sf(data = akde_sf |> 
            st_transform(crs = 4326) |> 
            filter(level == "95%"),
          aes(color = factor(id_orig)), fill = NA, linewidth = 0.75) +
  # Plot 50% isopleths
  geom_sf(data = akde_sf |> 
            st_transform(crs = 4326) |> 
            filter(level == "50%"),
          aes(color = factor(id_orig), fill = factor(id_orig)), linewidth = 0.75, alpha = 0.7) +
  scale_color_brewer("ID", palette = "Dark2") +
  scale_fill_brewer("ID", palette = "Dark2") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))







###########################################
### Measure space use overlap among IDs ###
###########################################

# Quantify Bhattacharyya's Affinity among UDs (w/ uncertainty)
overlap(akde)

# Show only the average estimated overlap
overlap(akde)$CI[,,"est"]






####################################################
### Assess pairwise interactions between animals ###
####################################################

### Calculate conditional distribution of encounters (CDE)
#a spatially defined PDF that describes the long-term encounter location probabilities for movement within home ranges

# Calculate for each pair
ctmm_cde <- combn(names(akde), 2, simplify = FALSE) |>  #create list of unique combos
  map(function(pair) {
    id1 <- pair[1]
    id2 <- pair[2]
    
    cde(akde[c(id1, id2)]) 
  }) |> 
  (\(x) set_names(x, map_chr(x, ~{.x@info$identity})))()  #set pair names per element

map(ctmm_cde, summary)

# Viz regions of potential encounters
plot(dat.telem2[c('5605_5','6471_2')], UD = ctmm_cde$`5605_5 6471_2`, col = c('red','green'))
plot(dat.telem2[c('6471_12','6471_6')], UD = ctmm_cde$`6471_12 6471_6`, col = c('red','blue'))
plot(dat.telem2[c('6469_2','6469_6')], UD = ctmm_cde$`6469_2 6469_6`, col = c('green','blue'))







########################################################
### Generate population-level inference of space use ###
########################################################

# mean() estimates population range; however, doesn't extrapolate to rest of pop.
#assumption that the sample of distributions is the population of distributions - such as for combining summer and winter ranges
pop_mean <- mean(akde)

summary(pop_mean)
plot(pop_mean)
plot(dat.telem2[idx], col = rainbow(length(dat.telem2[idx])), add = T)



# pkde() estimates range for entire pop. (via extrapolation) by accounting for inter-individual variability
#for when you have a small sample of a larger population - such as for a herd or colony

tic()
pkde <- pkde(data = dat.telem2[idx], UD = akde, kernel = "individual", ref = "Gaussian", weights = TRUE,
             dt = dt1, population = 3)
toc()  #took 10 min

summary(pkde)
plot(dat.telem2[idx], UD = pkde, col = rainbow(length(dat.telem2[idx])))
#in this case, results actually seem smaller compared to mean(); likely related to use of "messy" bursts being used

# Another way to viz pkde results
pkde_rast <- terra::rast(raster::raster(pkde))  #convert to SpatRaster
pkde_contours <- terra::as.contour(pkde_rast, levels = c(0.50, 0.95)) |>  #extract contours
  st_as_sf() |>  #convert to sf object
  st_cast("POLYGON")

ggplot() + 
  tidyterra::geom_spatraster(data = pkde_rast) + 
  scale_fill_viridis_c(option = "inferno", direction = -1) + 
  geom_point(data = dat3, aes(x, y), color = "chartreuse", size = 0.1, alpha = 0.25) +
  geom_sf(data = pkde_contours, aes(color = factor(level)), fill = NA, linewidth = 1) +
  scale_color_brewer("Level", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = range(dat3$x),
           ylim = range(dat3$y))
#Despite errors from plot() about grid size, results look the same



# meta() averages home range areas (in square meters) to get the population average and coefficient of variation
pop_meta <- ctmm::meta(akde, level.UD = 0.95, sort = TRUE)
#need to include 'ctmm' namespace explicitly due to conflict w/ another pkg

pop_meta






################################################################
### Make predictions of regularized tracks from fitted CTMMs ###
################################################################

### Predict from CTMMs
tic()
dat_pred <- map2(.x = ctmm_fit_best,
                 .y = dat.telem2[idx],
                 .f = ~ctmm::predict(object = .x, data = .y, dt = 3600, complete = TRUE))
toc()  # takes 5 s to run

# Convert from list to df
dat_pred_df <- dat_pred |> 
  map(data.frame) |> 
  bind_rows(.id = "id") |> 
  mutate(id_orig = as.vector(str_match(id, "[0-9]+")),
         .after = id) |> 
  dplyr::select(id, id_orig, timestamp, longitude, latitude, x, y, vx, vy)


# Viz regularized tracks
ggplot() +
  geom_point(data = dat_pred_df, aes(x, y, color = factor(id_orig)), size = 0.1) +
  scale_color_brewer("ID", palette = "Dark2") +
  coord_equal() +
  theme_bw()




### Simulate from CTMMs
tic()
sim_5605_1 <- ctmm::simulate(object = ctmm_fit_best$`5605_1`, nsim = 5, seed = 2026,
                           data = dat.telem2$`5605_1`, dt = 3600, complete = TRUE)
toc()  #took 0.25 sec

sim_5605_1_df <- sim_5605_1 |> 
  map(data.frame) |> 
  bind_rows(.id = "sim")


# Viz simulated tracks for burst 5605_1
ggplot() +
  geom_path(data = sim_5605_1_df, aes(x, y, color = factor(sim)), linewidth = 0.1, alpha = 0.5) +
  scale_color_brewer("Simulation #", palette = "Set1") +
  guides(color = guide_legend(override.aes = list(linewidth = 1, alpha = 1))) +
  coord_equal() +
  theme_bw()
#they're all identical in this case





#####################################
### Perform workflow in Shiny app ###
#####################################

#pak::pak("ctmm-initiative/ctmmweb")
ctmmweb::app(dat.telem2)





##########################################
#### Export datasets for easy loading ####
##########################################

save(akde, file = "processed_data/Session_2/AKDE_fits.RData")  #fitted AKDE
save(akde_sf, file = "processed_data/Session_2/AKDE_contours.RData")  #contours
save(ctmm_fit_best, file = "processed_data/Session_2/CTMM_fits.RData")  #fitted CTMMs
save(pkde, file = "processed_data/Session_2/PKDE_fits.RData")  #fitted PKDE
