
### Fit (integrated) step-selection functions (i)SSF ###

library(tidyverse)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(sf)
library(tictoc)
library(amt)

source("R/utils.R")



#-- While fitting any SSFs, it is a good idea to use observations recorded at the primary time interval (but especially so when fitting an iSSF). There are multiple ways to approach this, such as 1) using built-in functions from {amt} to find the steps at some defined primary time interval and specify some tolerance if not exact, 2) use best-fitting regularized tracks from any of the continuous-time models (via aniMotum, crawl, or ctmm), or 3) use the multiple imputations per track to account for location error and irregular time steps. --#

## For simplicity, we'll use the best-fitting track from the SSM (separated into bursts per ID)



###################
#### Load data ####
###################

# Load tracks
dat <- read_csv("processed_data/Session_3/regularized_tracks.csv")

glimpse(dat)  #coords (x,y) are in km
summary(dat)  #no missing values


# Change burst ID to numeric value
dat2 <- dat |> 
  rename(burst = id,
         id = id_orig) |> 
  mutate(burst = as.numeric(sub(".*_", "", burst)))



###########################################################
### Convert to {amt} 'track' object and prepare for SSF ###
###########################################################

# Check time intervals
dat2 |> 
  group_by(interaction(id, burst)) |>  #need to use interaction() since bursts are no longer unique
  mutate(dt = c(as.numeric(diff(date, units = "hour")), NA)) |> 
  pull(dt) |> 
  summary()  #looks good; all at 1 h intervals


# Convert directly to steps using pre-existing burst column
track_steps <- dat2 |> 
  make_track(.x = x, .y = y, .t = date, id = id, burst_ = burst, crs = 32736) |>  #"burst_" col must have trailing underscore to match required {amt} syntax
  nest(data = -id) |> 
  ## This is what you would do if resampling time intervals
  # mutate(steps = map(data, ~{
  #   .x |> track_resample(rate = hours(1), tolerance = minutes(15)) |> steps_by_burst()
  # })) |> 
  mutate(steps = map(data, steps_by_burst)) |>  #create step metrics for our pre-defined bursts; otherwise just need the steps() function
  select(id, steps) |>  #leave out "data" col
  unnest(cols = steps)

summary(track_steps)


# Viz distribution of step lengths and turning angles
track_steps |> 
  pivot_longer(cols = c(sl_,ta_), names_to = "metric", values_to = "values") |> 
  
  ggplot() +
  geom_density(aes(values, color = factor(id)), alpha = 0.5) +
  theme_bw() +
  facet_wrap(~metric, scales = "free")
#all nearly identicial across IDs




# Generate 10 random steps per observed step (same ratio as for RSF)

#As with RSF (and maybe even moreso), these points should be thought of as an approach to approximate the integral of nearby habitat on the landscape. See Michelot et al (2024) "Understanding step selection analysis through numerical integration" for a more detailed explanation

tic()
track_steps_presabs <- track_steps |> 
  random_steps(n_control = 10)  #assumes Gamma and Von Mises distribs for SL and TA, respectively, by default
toc()  #took 12 sec

track_steps_presabs
# We now have columns 'case_' and 'step_id_', referring to whether the point was observed (TRUE) or randomly sampled (FALSE) and the 'strata' or step number per ID/burst
#'step_id_' must start at 3 because turning angles require at least 3 consecutive observations






##########################
### Extract covariates ###
##########################

### Load covars
dist2water <- rast("rasters/dist2water_rsf.tif")
dist2pop <- rast("rasters/dist2pop_rsf.tif")
ndvi <- rast("rasters/ndvi_rsf.nc")

# Convert rasters to same proj as steps (i.e., units must be in km)
dist2water_proj <- project(dist2water, "+proj=utm +zone=36 +ellps=WGS84 +units=km +no_defs +south",
                           threads = 10, by_util = TRUE)
dist2pop_proj <- project(dist2pop, "+proj=utm +zone=36 +ellps=WGS84 +units=km +no_defs +south",
                         threads = 10, by_util = TRUE)
ndvi_proj <- project(ndvi, "+proj=utm +zone=36 +ellps=WGS84 +units=km +no_defs +south",
                     threads = 10, by_util = TRUE)

# Create SpatRaster stack for static covars
static_covars <- c(dist2water_proj, dist2pop_proj)
names(static_covars) <- c("dist2water","dist2pop")


### For static covars
track_steps_covars <- track_steps_presabs |> 
  # Extract static covars
  extract_covariates(covariates = static_covars, where = "end") |>  #extract covar values at "end" of each step
  # Extract dynamic covars (we specify that we're looking for nearest time w/in 1 week)
  extract_covariates_var_time(covariates = ndvi_proj, when = "any", max_time = weeks(x = 1),
                              name_covar = "ndvi", where = "end") |> 
  # Z-scale all covars
  mutate(across(.cols = c(dist2pop, dist2water, ndvi),
                .fns = ~{scale(.x) |>  #transform w/ scale() function
                    c()},
                .names = "{.col}_s"))  #create new cols for scaled covars





###############
### Fit SSF ###
###############

# Fit simple SSF
fit_ssf1 <- track_steps_covars |> 
  fit_clogit(case_ ~ dist2pop_s + dist2water_s + ndvi_s + strata(step_id_), model = TRUE)  #fit_ssf() could also be used

summary(fit_ssf1)
#these results differ from the RSF
#they suggest that (at the "step" scale) male elephants choose to be at greater distances from human settlements, closer to water sources, and in areas of high NDVI

# Viz estimated coeffs
ssf_coefs <- broom::tidy(fit_ssf1$model, exponentiate = TRUE, conf.int = TRUE) |> 
  filter(term %in% c("dist2water_s", "dist2pop_s", "ndvi_s"))

ggplot(ssf_coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3, color = "#2c7fb8") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkred") +
  labs(
    title = "Simple SSF Parameter Estimates",
    x = "Relative Selection Strength (RSS per 1 SD increase)",
    y = NULL
  ) +
  theme_minimal(base_size = 14)



# Viz RSS
covars <- c("dist2water_s", "dist2pop_s", "ndvi_s")

ssf_preds <- predict_ssf_margeff(fit = fit_ssf1, focal_covars = covars, data = track_steps_covars)

# Plot RSS
ggplot(ssf_preds, aes(x = x_natural, y = rss)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = rss_lwr, ymax = rss_upr), fill = "#2c7fb8", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "Simple SSF Marginal Effects",
    subtitle = "Relative Selection Strength (RSS) relative to average available step endpoint",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.placement = "outside", strip.text = element_text(face = "bold"))






################
### Fit iSSF ###
################

track_steps_covars2 <- track_steps_covars |> 
  mutate(log_sl_ = log(sl_),      # Required for updated gamma step length parameters
         cos_ta_ = cos(ta_)       # Required for updated von Mises turn angle parameters
         )

### Fit simple integrated SSF
fit_issf1 <- track_steps_covars2 |> 
  fit_issf(case_ ~ dist2pop_s + dist2water_s + ndvi_s +  #env covars
             sl_ + log_sl_ + cos_ta_ +  #movement covars
             strata(step_id_),  #strata for each step
           model = TRUE)  #fit_clogit could also be used, but with less post-fitting functionality

summary(fit_issf1)
#these results are essentially the same as for SSF (but won't always be)



### Fit iSSF w/ interactions
fit_issf2 <- track_steps_covars2 |> 
  fit_issf(case_ ~ dist2pop_s + dist2water_s + ndvi_s +  #env covars
             sl_ + log_sl_ + cos_ta_ +  #movement covars
             sl_:dist2pop_s + sl_:dist2water_s +  #interactions of SL w/ dist covars
             strata(step_id_),  #strata for each step
           model = TRUE)  #fit_clogit could also be used, but with less post-fitting functionality

summary(fit_issf2)
#these results slightly differ now, w/ a signif interaction between dist2water and SL




# Viz estimated coeffs
issf_coefs <- broom::tidy(fit_issf2$model, exponentiate = TRUE, conf.int = TRUE)

ggplot(issf_coefs, aes(x = estimate, y = term)) +
  geom_point(size = 3, color = "#2c7fb8") +
  geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "darkred") +
  labs(
    title = "iSSF Parameter Estimates",
    x = "Relative Selection Strength (RSS per 1 SD increase)",
    y = NULL
  ) +
  theme_minimal(base_size = 14)




# Viz RSS
covars <- c("dist2water_s", "dist2pop_s", "ndvi_s")

issf_preds <- predict_ssf_margeff(fit = fit_issf2, focal_covars = covars, data = track_steps_covars2)

# Plot RSS
ggplot(issf_preds, aes(x = x_natural, y = rss)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = rss_lwr, ymax = rss_upr), fill = "#2c7fb8", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "iSSF Marginal Relative Selection Strength (RSS)",
    subtitle = "Evaluated at mean step length with interactions incorporated",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.placement = "outside", strip.text = element_text(face = "bold"))




### Explore relationship between dist2water and SL

# Extract tentative distribution parameters using sl_distr()
tentative_sl <- sl_distr(fit_issf2)
alpha_0 <- tentative_sl$params$shape
scale_0 <- tentative_sl$params$scale
beta_0  <- 1 / scale_0  # Rate parameter

# Extract model coeffs
b_sl       <- coef(fit_issf2)["sl_"]
b_log_sl   <- coef(fit_issf2)["log_sl_"]
b_sl_water <- coef(fit_issf2)["dist2water_s:sl_"]

# Calculate expected mean SL across dist2water gradient
water_grid_s <- seq(min(track_steps_covars2$dist2water_s), max(track_steps_covars2$dist2water_s), length.out = 100)

mean_sl_df <- data.frame(dist2water_s = water_grid_s) |> 
  mutate(
    dist2water_m  = (dist2water_s * sd(track_steps_covars2$dist2water)) + mean(track_steps_covars2$dist2water),
    beta_sl_eff   = b_sl + (b_sl_water * dist2water_s),
    alpha_updated = alpha_0 + b_log_sl,
    beta_updated  = beta_0 - beta_sl_eff,
    # Expected mean step length: E[X] = shape / rate
    mean_sl       = alpha_updated / beta_updated
  )

# Plot expected mean SL over dist2water grad
ggplot(mean_sl_df, aes(x = dist2water_m / 1000, y = mean_sl)) +
  geom_line(color = "#2c7fb8", linewidth = 1.2) +
  labs(
    title = "Expected Step Length Conditional on Distance to Water",
    subtitle = "Calculated from updated selection-free Gamma parameters",
    x = "Distance to Water (km)",
    y = "Expected Mean Step Length (km)"
  ) +
  theme_minimal(base_size = 14)
#steps are ~30% longer when elephants are at water vs 57 km away





##########################################################
### Simulate from iSSF to generate spatial predictions ###
##########################################################

# Create mean NDVI layer for prediction
ndvi_proj_mean <- mean(ndvi_proj, na.rm = TRUE)
names(ndvi_proj_mean) <- "ndvi"

# Combine all rasters together
rast_stack <- c(static_covars, ndvi_proj_mean)
rast_stack <- rast_stack[[c(2,1,3)]]  #reorder to match coeff order in iSSF formula

# Scale rasters
rast_stack_s <- rast_stack
names(rast_stack_s) <- paste0(names(rast_stack_s), "_s")

rast_stack_s$dist2water_s <- (rast_stack_s$dist2water_s - mean(track_steps_covars2$dist2water, na.rm = TRUE)) / sd(track_steps_covars2$dist2water, na.rm = TRUE)
rast_stack_s$dist2pop_s <- (rast_stack_s$dist2pop_s - mean(track_steps_covars2$dist2pop, na.rm = TRUE)) / sd(track_steps_covars2$dist2pop, na.rm = TRUE)
rast_stack_s$ndvi_s <- (rast_stack_s$ndvi_s - mean(track_steps_covars2$ndvi, na.rm = TRUE)) / sd(track_steps_covars2$ndvi, na.rm = TRUE)



set.seed(2026)

# Get unique animal IDs
animal_ids <- unique(dat2$id)

# Iterate through each ID, calculate its kernel, and run 10 simulations each
#Example from Johannes Signer shown here: https://github.com/jmsigner/amt/issues/93
tic()
sim_paths_all <- map_dfr(animal_ids, function(animal_id) {  #map across IDs
  
  # Extract starting location for current animal (although we could randomly sample from study area)
  start_pt <- dat2 |> 
    filter(id == animal_id) |> 
    make_track(.x = x, .y = y, .t = date, crs = 32736) |> 
    slice(1) |> 
    make_start()
  
  # Calculate redistribution kernel specific to this start point
  kernel <- redistribution_kernel(
    x     = fit_issf2, 
    start = start_pt,
    map   = rast_stack_s, 
    fun   = function(xy, map) {
      extract_covariates(xy, map, where = "end") |> 
        mutate(
          log_sl_ = log(sl_),
          cos_ta_ = cos(ta_)
        )
    }
  )
  
  # Generate 10 simulated paths of 200 steps for this animal (~10 days)
  #This will create a transient UD
  #For a steady-state UD, either increase `n` and/or choose many more random starting points
  map_dfr(1:10, function(sim_i) {  #Map across sim number
    simulate_path(kernel, n.steps = 200) |> 
      mutate(id = animal_id, sim = sim_i)
  })
  # replicate(10, simulate_path(kernel, n.steps = 200), simplify = FALSE)
  
}, .progress = TRUE)

toc()  #took 36.5 min


# Convert all pooled simulation locations (30 total paths) to sf points
sim_sf <- make_track(sim_paths_all, .x = x_, .y = y_, crs = 32736) |> 
  as_sf_points()

# Rasterize point counts across all animals and paths onto the landscape grid
#alternatively, the UD estimation could also be done w/ hr_kde() for contours of 'track' objects
ud_counts <- rasterize(
  x          = sim_sf, 
  y          = rast_stack_s[[1]], 
  fun        = "count", 
  background = 0
)

# Normalize cell values so total landscape occupancy sums to 1
ud_raster <- ud_counts / sum(values(ud_counts), na.rm = TRUE)

# Quick plot of combined population-level simulated UD
plot(ud_raster, main = "Combined Simulated Population UD (3 Animals, 10 Paths Each)")


# Compare SIMULATED tracks against predicted UD
ggplot() +
  geom_spatraster(data = ud_raster) +
  scale_fill_viridis_c("Est. UD") +
  geom_path(data = sim_paths_all, aes(x_, y_, group = interaction(id, sim)), alpha = 0.3, color = "white",
            linewidth = 0.2) +
  theme_bw(base_size = 14) +
  coord_sf(xlim = range(sim_paths_all$x_),
           ylim = range(sim_paths_all$y_))


# Compare OBSERVED tracks against predicted UD
ggplot() +
  geom_spatraster(data = ud_raster) +
  scale_fill_viridis_c("Est. UD") +
  geom_path(data = dat2 |> 
              slice(1:200, .by = id), aes(x, y, group = id), alpha = 0.5, color = "white",
            linewidth = 0.5) +
  theme_bw(base_size = 14) +
  coord_sf(xlim = range(sim_paths_all$x_),
           ylim = range(sim_paths_all$y_))
#simulated tracks do show different movement patterns, and therefore predicted UDs
#behavioral state likely plays a large role in habitat selection and movement process not captured by this model

# ggsave(filename = "website/images/iSSF_pred_UD.png", units = "in", width = 6, height = 4, dpi = 400)



##############
### Export ###
##############

# Save simulated paths
write_csv(sim_paths_all, "processed_data/Session_5/iSSD_sim_paths.csv")

# Save predicted UD layer
writeRaster(ud_raster, "processed_data/Session_5/iSSF_UD_200steps.tif")
