
### Estimate effect of environment on space use ###

library(tidyverse)
library(sf)
library(rnaturalearth)
library(terra)
library(tidyterra)
library(exactextractr)
library(brms)
library(rstan)
library(bayesplot)
library(tidybayes)

source("R/utils.R")

# Enable parallel processing for brms
options(mc.cores = parallel::detectCores())
rstan_options(auto_write = TRUE)



#################
### Load data ###
#################

# Load tracks
dat <- read_csv('processed_data/cleaned_tracks.csv')

# Add project coords and calc displacement per ID
dat <- dat |> 
  add_trans_coords(coords = c('lon','lat'), proj = 4326, new_proj = 32736) |> 
  group_by(id) |> 
  arrange(date, .by_group = TRUE) |> 
  mutate(disp = sqrt((x - x[1])^2 + (y - y[1])^2),  #calc distance from first location
         disp = disp / 1000) |>  #convert units from m to km
  ungroup()


# Load vector layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")


# Load AKDE contours (from segmented tracks)
#Note: since these tracks were segmented relatively subjectively, we may want to be cautious when interpreting the results
load("processed_data/AKDE_contours.RData")


# Load in environmental rasters
dist2water <- rast("rasters/dist2water_rsf.tif")  #at 440 m res
ndvi <- rast("rasters/ndvi_rsf.nc")  #at 440 m res
shrub_prop <- rast("rasters/shrub_prop_rsf.tif")  #at 440 m res




# Viz contours
ggplot() +
  geom_sf(data = africa |> st_transform(32736)) +
  geom_path(data = dat, aes(x, y, group = id, color = factor(id_orig)), linewidth = 0.2, alpha = 0.5) +
  geom_sf(data = akde_sf |> 
            filter(level == "95%"), aes(color = factor(id_orig)), fill = NA, linewidth = 1) +
  scale_color_brewer("ID", palette = "Dark2") +
  labs(title = "95% KDE") +
  theme_bw() +
  coord_sf(xlim = st_bbox(dat_id_kde_href2)[c("xmin","xmax")],
           ylim = st_bbox(dat_id_kde_href2)[c("ymin","ymax")])


# Let's explore how home range size differs in relation to distance to water, proportion shrub cover, and NDVI variability




##################################
### Segment tracks into bursts ###
##################################

### Split tracks exact same way as done for AKDE

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
dat2 <- dat |> 
  inner_join(burst_windows3,
             by = join_by(id, between(date, start, end))) |> 
  filter(phase == 'resident') |> 
  mutate(burst_id = paste(id, burst, sep = "_"),
         .after = id)






###########################
### Prep data for model ###
###########################

# Define start and end dates per burst
burst_windows_resident <- burst_windows3 |> 
  filter(phase == 'resident') |> 
  mutate(burst_id = paste(id, burst, sep = "_"),
         .after = id)



### Calculate NDVI variability for each bursts time window
#as SD across time, averaged across all pixels w/in contour

# Join metadata time windows into the 95% AKDE sf object
akde_95 <- akde_sf |>
  filter(level == "95%") |>
  inner_join(
    burst_windows_resident |> select(burst_id, start, end),
    by = c("id" = "burst_id")
  )

# Extract NDVI times
ndvi_times <- time(ndvi)

# Iterate through each burst UD to calculate variability and extract
akde_ndvi_summary <- akde_95 |>
  mutate(
    # Extract spatial mean of temporal standard deviation
    mean_temporal_sd = pmap_dbl(
      list(st_geometry(akde_95), start, end),
      function(geom, start_t, end_t) {
        
        # Identify raster layer indices matching the burst time window
        idx <- which(ndvi_times >= start_t & ndvi_times <= end_t)
        
        if (length(idx) == 0) return(NA)
        
        # Temporal subset of NDVI layers
        ndvi_sub <- ndvi[[idx]]
        
        # Pixel-wise temporal variability (Standard Deviation across time)
        temporal_sd <- app(ndvi_sub, fun = "sd", na.rm = TRUE)
        
        # Spatial extraction within polygon
        poly <- st_sf(geometry = st_sfc(geom, crs = st_crs(akde_95)))
        exact_extract(temporal_sd, poly, fun = "mean", progress = TRUE)
      }
    )
  )




# Create stack of static rasters
static_rast <- c(dist2water, shrub_prop)
names(static_rast) <- c("dist2water", "shrub_prop")


# Extract mean values of static covars per polygon
akde_rast_ext <- akde_ndvi_summary |> 
  cbind(exact_extract(static_rast, akde_ndvi_summary, fun = "mean", progress = TRUE)) |> 
  rename(ndvi_var = mean_temporal_sd,
         dist2water = mean.dist2water,
         shrub_prop = mean.shrub_prop)

# Scale all covars
akde_rast_ext <- akde_rast_ext |> 
  mutate(across(.cols = c(ndvi_var, dist2water, shrub_prop),
                .fns = ~{scale(.x) |>  #transform w/ scale() function
                    c()},
                .names = "{.col}_s"))  #create new cols for scaled covars


# Add col for area of each UD
akde_rast_ext2 <- akde_rast_ext |> 
  mutate(area = as.numeric(st_area(geometry)) / 1e6)  #convert to km^2





##############################
### Prior predictive check ###
##############################

# Inspect default/suggested priors by {brms}
get_prior(
  area ~ ndvi_var_s + dist2water_s + shrub_prop_s + (1 | id_orig),  #this is our model formula
  data = akde_rast_ext2,
  family = lognormal()  #using log Normal distrib since area must be strictly positive
)



# Define domain-informed, weakly informative priors (to regularize estimated params)
#Log link means coefficients represent proportional changes on log scale
log(mean(akde_rast_ext2$area))  # ~3500 km; 8.16 on log scale

priors <- c(
  prior(normal(8.16, 1.5), class = "Intercept"), # Baseline home range size
  prior(normal(0, 0.5), class = "b"),              # Fixed effects (on log scale)
  prior(exponential(2), class = "sd"),           # Between-individual variability
  prior(exponential(2), class = "sigma")        # Residual variability
)

# Fit PRIOR-ONLY model (Sampling from Regularizing Priors)
fit_prior <- brm(
  formula = area ~ ndvi_var_s + dist2water_s + shrub_prop_s + (1 | id_orig),
  data = akde_rast_ext2,
  family = lognormal(),
  prior = priors,   #supply pre-defined priors
  sample_prior = "only", # IGNORES DATA; draws only from prior distributions
  chains = 2,
  cores = 2,
  iter = 2000,
  seed = 2026
)


summary(fit_prior)  #values look pretty good


# Prior Predictive Check
# Evaluate if defined priors produce biologically plausible areas
pp_check(fit_prior, ndraws = 50) +
  scale_x_log10() +
  labs(title = "Prior Predictive Check (Regularizing Priors)",
       subtitle = "Are simulated home range areas biologically plausible?")
#yep, the data (black) matches relatively well w/ the draws from the posterior (blue)




##############################################
### Fit full model (sample from posterior) ###
##############################################

fit_hr_area <- brm(
  formula = area ~ ndvi_var_s + dist2water_s + shrub_prop_s + (1 | id_orig),
  data = akde_rast_ext2,
  family = lognormal(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 2000,
  warmup = 1000,
  seed = 2026,
  control = list(adapt_delta = 0.9) # Helps prevent divergent transitions
)
# took 5 sec



###########################################
### Check model convergence diagnostics ###
###########################################

# Check Numerical Diagnostics (Rhat < 1.01, Bulk/Tail ESS > 1000)
summary(fit_hr_area)
#all values look pretty good

# Visual Convergence - Traceplots (Fuzzy caterpillars indicate convergence)
plot(fit_hr_area)

# Rank Plots (Uniform distribution across chains indicates proper mixing)
mcmc_rank_overlay(fit_hr_area)





##################################
### Posterior predictive check ###
##################################

# Density Overlay (Does model simulated data match observed data distribution?)
pp_check(fit_hr_area, ndraws = 100) +
  scale_x_log10() +
  labs(title = "Posterior Predictive Check", subtitle = "Observed (y) vs. Model Simulated (yrep)")
#yes, looks pretty close!

# Summary Statistic Checks (e.g., checking if model captures maximum step length)
pp_check(fit_hr_area, type = "stat", stat = "max", ndraws = 100)
pp_check(fit_hr_area, type = "stat", stat = "median", ndraws = 100)
#model currently estimates max and median home range areas well






##############################################################
### Visualize the posterior and make ecological inferences ###
##############################################################

### Plotting Half-Eye Density Interval Distributions of Parameters

# On log (link) scale
fit_hr_area |> 
  gather_draws(b_ndvi_var_s, b_dist2water_s, b_shrub_prop_s) |> 
  
  ggplot(aes(x = .value, y = .variable, fill = .variable)) +
  stat_halfeye(point_interval = median_hdi, .width = c(0.89, 0.95)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkred") +
  labs(
    title = "Posterior Probability Distributions of Fixed Effects",
    x = "Effect Size (log-scale)",
    y = "Parameter"
  ) +
  theme_bw(base_size = 14) +
  theme(legend.position = "none")
#all covars show no effect






### Simple Conditional Marginal Effects Plots
conditional_effects(fit_hr_area, effects = c("ndvi_var_s", "dist2water_s", "shrub_prop_s"))



### Custom Conditional Marginal Effects Plots

# Generate prediction list
ce <- conditional_effects(
  fit_hr_area, 
  effects = c("ndvi_var_s", "dist2water_s", "shrub_prop_s"), 
  plot = FALSE
)

# Un-scale x-values and combine into a single tidy data.frame
effects_df <- map_dfr(names(ce), function(eff) {
  raw_var <- sub("_s$", "", eff)
  
  var_mean <- mean(akde_rast_ext2[[raw_var]], na.rm = TRUE)
  var_sd   <- sd(akde_rast_ext2[[raw_var]], na.rm = TRUE)
  
  ce[[eff]] |> 
    mutate(
      covariate = raw_var,
      x_natural = (.data[[eff]] * var_sd) + var_mean
    )
})

# Plot all marginal effects
ggplot(effects_df, aes(x = x_natural, y = estimate__)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = lower__, ymax = upper__), fill = "#2c7fb8", alpha = 0.2) +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "Bayesian Marginal Effects (Natural Scale)",
    y = "Expected Home Range Area (km^2)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside", 
    strip.text = element_text(face = "bold")
  )
#we see essentially no effect of any of these covariates
##these results suggest that differences in home range size aren't captured by these covariates, which may also differ by season, age, or other environmental factors at various scales

