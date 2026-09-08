
### Fit Bayesian resource selection function (RSF) using GLM ###
# Also including spatial random effect (as Gaussian Random Field) #

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(inlabru)
library(fmesher)
library(tictoc)
library(patchwork)

source("R/utils.R")



###################
#### Load data ####
###################

# Load tracks
dat <- read_csv("processed_data/Session_5/dat_presabs_covars.csv") |> 
  mutate(id = as.character(id))

glimpse(dat)
summary(dat)


# Load vector layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")



# Viz tracks
ggplot() +
  geom_sf(data = africa) +
  geom_sf(data = gl_pa, color = "black", linewidth = 0.5) +
  geom_path(data = dat |> 
              filter(obs == 1),  #only keep the observed locs
            aes(lon, lat, group = id, color = id), linewidth = 0.25, alpha = 0.5) +
  theme_bw(base_size = 14) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top") +
  labs(x = "Easting", y = "Northing") +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))





##############################
### Prep data for modeling ###
##############################

# For now, we're just going to use 3 covars (none of the proportions)

# Remove any obs w/ NAs for covars
dat2 <- dat |> 
  drop_na(dist2pop, dist2water, ndvi)

# Z-score transform covars (i.e., center and scale by SD)
dat2 <- dat2 |> 
  mutate(across(.cols = c(dist2pop, dist2water, ndvi),
                .fns = ~{scale(.x) |>  #transform w/ scale() function
                    c()},
                .names = "{.col}_s"))  #create new cols for scaled covars


# Define weights for fitting as inhomogeneous Poisson point process (IPP)
# Referred to as infinitely weight logistic regression (IWLR)
# Standard choice for IWLR RSF model: 'Used' points weight = 1, background points weight = 5000 (or other large constant weight)
dat2$wts <- ifelse(dat2$obs == 1, 1, 5000)





###########################
### Fit RSF as IWLR GLM ###
###########################

# Define components
cmp <- ~ Intercept(1) + dist2pop_s + dist2water_s + ndvi_s

# Logit link logistic regression
tic()
fit_rsf_linear <- bru(
  components = cmp,
  formula = obs ~ Intercept + dist2pop_s + dist2water_s + ndvi_s,
  data = dat2,
  family = "binomial",
  weights = dat2$wts,
  options = list(
    control.inla  = list(int.strategy = "eb", strategy = "adaptive"),
    control.compute = list(openmp.strategy = "default"),
    inla.mode = "experimental",
    verbose = TRUE,
    bru_verbose = 4
  )
)
toc()  #took 3 sec

summary(fit_rsf_linear)
# Identical to the version we fit before w/ glm()
# For IWLR, we ignore the intercept since this is directly tied to N_avail


### Predict marginal effects
covars <- c("dist2pop_s", "dist2water_s", "ndvi_s")

# Generate predictions for all covariates and combine
linear_margeff <- map_dfr(covars,
                          get_inla_margeff, fit = fit_rsf_linear, dat_orig = dat2)

# Plot faceted marginal curves
ggplot(linear_margeff, aes(x = x_natural, y = mean)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "#2c7fb8", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "Marginal Effects of RSF Covariates",
    subtitle = "Relative Selection Strength (RSS) with other covariates held at their mean",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.text = element_text(face = "bold")
  )






######################################
### Fit RSF as mixed-effects model ###
######################################

# Ensure 'id' is formatted correctly (factor or integer)
dat2$id <- factor(dat2$id)

# 2. Define components with fixed effects + random intercepts and slopes by 'id'
cmp_mixed <- ~ Intercept(1) +
  dist2pop_s + dist2water_s + ndvi_s +
  
  # Random Intercept by ID; use priors per Muff et al. (2020)
  id_int(id, model = "iid", hyper = list(theta = list(initial = log(1e-6), fixed = TRUE))) +
  
  # Random Slopes by ID (covariate passed via 'weights'); use priors per Muff et al. (2020)
  id_dist2pop(id, weights = dist2pop_s, model = "iid",
              hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  id_dist2water(id, weights = dist2water_s, model = "iid",
                hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05)))) +
  id_ndvi(id, weights = ndvi_s, model = "iid",
          hyper = list(theta = list(initial = log(1), fixed = FALSE, prior = "pc.prec", param = c(1,0.05))))


# Fit the mixed-effects RSF
tic()
fit_rsf_mixed <- bru(
  components = cmp_mixed,
  formula    = obs ~ Intercept + dist2pop_s + dist2water_s + ndvi_s +
    id_int + id_dist2pop + id_dist2water + id_ndvi,
  family     = "binomial",
  data       = dat2,
  weights = dat2$wts,
  options = list(
    control.inla  = list(int.strategy = "eb", strategy = "adaptive"),
    control.compute = list(openmp.strategy = "default"),
    # use priors per Muff et al. (2020)
    control.fixed = list(mean = 0,
                         prec = list(dist2pop_s = 1e-4, dist2water_s = 1e-4, ndvi_s = 1e-4)),
    inla.mode = "experimental",
    verbose = TRUE,
    bru_verbose = 4
  )
)
toc()  #took 7 sec

summary(fit_rsf_mixed)
#also similar to frequentist model



# Generate predictions for all covariates and combine
linear_mixed_margeff <- map_dfr(covars,
                          get_inla_margeff, fit = fit_rsf_mixed, dat_orig = dat2)

# Plot faceted marginal curves
ggplot(linear_mixed_margeff, aes(x = x_natural, y = mean)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "#2c7fb8", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "Marginal Effects of RSF Covariates",
    subtitle = "As mixed-effects model",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.text = element_text(face = "bold")
  )




#########################################
### Fit RSF w/ spatial random effects ###
#########################################

# Build 2D spatial mesh from coordinate bounds (ideally UTM coordinates in meters/kilometers)
mesh_2d <- fm_mesh_2d(
  loc      = as.matrix(dat2[dat2$obs == 1, c('x','y')]),  #only on observed locs
  max.edge = c(10000, 50000), # Inner and outer domain edge lengths (in meters)
  cutoff   = 5000            # Minimum distance between mesh nodes (in meters)
)

mesh_2d  #number of nodes should generally not be >1000, but definitely not >3000

# Viz 2D mesh
ggplot() +
  gg(mesh_2d) +
  geom_point(data = dat2 |> 
               filter(obs == 1), aes(x, y), shape = ".", size = 0.5) +
  theme_minimal()
#mesh looks good compared to the points
#ideally, don't want mesh to be too fine or too coarse
# See Dambly et al (2023) "Integrated species distribution models fitted in INLA are sensitive to mesh parameterisation" (https://doi.org/10.1111/ecog.06391) for suggestions when building a mesh



# Define Matérn Gaussian Process model (SPDE)
# Prior assumption: Pr(spatial range) > 10,000 m = 0.05; Pr(marginal SD) < 1 = 0.01
spde_spatial <- INLA::inla.spde2.pcmatern(
  mesh        = mesh_2d,
  prior.range = c(10000, 0.05),
  prior.sigma = c(1, 0.01)
)




# Define model components
cmp_grf <- ~ Intercept(1) +
  dist2pop_s + dist2water_s + ndvi_s +
  # 2D mesh for spatial random effect
  grf(cbind(x, y), model = spde_spatial)

# Fit RSF w spatial random effect
tic()
fit_rsf_grf <- bru(
  components = cmp_grf,
  formula    = obs ~ Intercept + dist2pop_s + dist2water_s + ndvi_s + grf,
  family     = "binomial",
  data       = dat2,
  weights = dat2$wts,
  options = list(
    control.inla  = list(int.strategy = "eb", strategy = "adaptive"),
    control.compute = list(openmp.strategy = "default"),
    inla.mode = "experimental",
    verbose = TRUE,
    bru_verbose = 4
  )
)
toc()  #took 3 min

summary(fit_rsf_grf)
#results are VERY different from other RSFs
# strong negative relationship w/ dist2water and positive relationships for dist2pop and ndvi



# Generate predictions for all covariates and combine
linear_grf_margeff <- map_dfr(covars,
                                get_inla_margeff, fit = fit_rsf_grf, dat_orig = dat2)

# Plot faceted marginal curves
ggplot(linear_grf_margeff, aes(x = x_natural, y = mean)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = q0.025, ymax = q0.975), fill = "#2c7fb8", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free", strip.position = "bottom") +
  labs(
    title = "Marginal Effects of RSF Covariates",
    subtitle = "Spatial random effect model",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(
    strip.placement = "outside",
    strip.text = element_text(face = "bold")
  )
#by accounting for spatial autocorrelation in the points, we now estimate that these male elephants prefer habitats that are at least 15 km from human settlements, 15 km or closer to water sources, and in areas where NDVI > 0.5


# Viz spatial prediction of GRF
pxl <- fm_pixels(mesh_2d, mask = TRUE)

# Extract coordinates and assign explicit 'x' and 'y' columns
coords <- sf::st_coordinates(pxl)
pxl$x <- coords[, 1]
pxl$y <- coords[, 2]

# Predict GRF on grid
grf_pred <- predict(fit_rsf_grf, pxl, ~ grf)



#mean
ggplot() +
  gg(data = grf_pred, aes(fill = mean), geom = "tile") +
  scale_fill_viridis_c("GRF") +
  theme_bw() +
  coord_sf(xlim = fm_bbox(mesh_2d)[1] |> unlist(),
           ylim = fm_bbox(mesh_2d)[2] |> unlist())

#SD
ggplot() +
  gg(data = grf_pred, aes(fill = sd), geom = "tile") +
  scale_fill_viridis_c("SD GRF", option = "rocket") +
  theme_bw() +
  coord_sf(xlim = fm_bbox(mesh_2d)[1] |> unlist(),
           ylim = fm_bbox(mesh_2d)[2] |> unlist())




#########################################
### Make spatial prediction from RSFs ###
#########################################

## Load in environmental rasters
dist2water <- rast("rasters/dist2water_rsf.tif")
dist2pop <- rast("rasters/dist2pop_rsf.tif")
ndvi <- rast("rasters/ndvi_rsf.nc")

# Z-scale transform rasters
dist2water_s <- (dist2water - mean(dat2$dist2water)) / sd(dat2$dist2water)
dist2pop_s <- (dist2pop - mean(dat2$dist2pop)) / sd(dat2$dist2pop)
ndvi_s <- (ndvi - mean(dat2$ndvi)) / sd(dat2$ndvi)

# Calc mean of scaled NDVI
ndvi_s_mean <- mean(ndvi_s, na.rm = TRUE)




## Predict relative selection intensity (w(x) = exp(x * beta); or relative abundance) in space
#doesn't include "Intercept" in prediction

# Construct raster stack matching covariate names (w/ mean NDVI)
#likewise, time-varying values of NDVI could be used to make dynamic predictions of habitat selection
cov_stack <- c(dist2water_s, dist2pop_s, ndvi_s_mean)
names(cov_stack) <- c("dist2water_s", "dist2pop_s", "ndvi_s")

#Convert to data.frame
covs_df <- as.data.frame(cov_stack, xy = TRUE) |> 
  drop_na(dist2water_s:ndvi_s)


# Predict relative intensity (i.e., relative abundance)
#Also generally a good idea to show the uncertainty in predictions too, although this should be done w/o intercept
rsf_map_linear <- predict(fit_rsf_linear, newdata = covs_df,
                          formula = ~dist2water_s + dist2pop_s + ndvi_s,
                          num.threads = 10, n.samples = 50)
rsf_map_mixed <- predict(fit_rsf_mixed, newdata = covs_df,
                          formula = ~dist2water_s + dist2pop_s + ndvi_s,
                          num.threads = 10, n.samples = 50)
rsf_map_grf <- predict(fit_rsf_grf, newdata = covs_df,
                       formula = ~dist2water_s + dist2pop_s + ndvi_s,  #leave off GRF component
                       num.threads = 10, n.samples = 50) 





# Combine predictions
rsf_preds <- list(GLM = rsf_map_linear |> mutate(model = "GLM"),
                  GLMM = rsf_map_mixed |> mutate(model = "GLMM"),
                  GLM_GRF = rsf_map_grf |> mutate(model = "GLM + GRF"))



# Viz mapped predictions of RSS (mean)
rsf_mean_maps <- map(rsf_preds,
                ~{
                  ggplot() +
                    geom_raster(data = .x, aes(x, y, fill = mean)) +
                    scale_fill_viridis_c("log(Relative Intensity)", na.value = "transparent") +
                    geom_sf(data = africa |> 
                              st_transform(32736), fill = NA, color = "white", lwd = 1) +
                    geom_path(data = dat2 |> 
                                filter(obs == 1),  #only keep the observed locs
                              aes(x, y, group = id), color = "white", lwd = 0.25, alpha = 0.5) +
                    theme_bw(base_size = 14) +
                    theme(strip.text = element_text(face = "bold")) +
                    labs(x = "Easting", y = "Northing") +
                    coord_sf(xlim = range(rsf_map_linear$x),
                             ylim = range(rsf_map_linear$y),
                             expand = FALSE) +
                    facet_wrap(~ model)
                })


# Create composite plot
rsf_mean_maps$GLM + rsf_mean_maps$GLMM + rsf_mean_maps$GLM_GRF +
  plot_layout(ncol = 2)
#Each map slightly different, but GLM + GRF model has by far the strongest effect size


# Viz mapped prediction uncertainty of RSS (SE)
rsf_sd_maps <- map(rsf_preds,
                     ~{
                       ggplot() +
                         geom_raster(data = .x, aes(x, y, fill = sd)) +
                         scale_fill_viridis_c("SD log(Rel. Int.)", option = "rocket", na.value = "transparent") +
                         geom_sf(data = africa |> 
                                   st_transform(32736), fill = NA, color = "white", lwd = 1) +
                         geom_path(data = dat2 |> 
                                     filter(obs == 1),  #only keep the observed locs
                                   aes(x, y, group = id), color = "white", lwd = 0.25, alpha = 0.5) +
                         theme_bw(base_size = 14) +
                         theme(strip.text = element_text(face = "bold")) +
                         labs(x = "Easting", y = "Northing") +
                         coord_sf(xlim = range(rsf_map_linear$x),
                                  ylim = range(rsf_map_linear$y),
                                  expand = FALSE) +
                         facet_wrap(~ model)
                     })


# Create composite plot
rsf_sd_maps$GLM + rsf_sd_maps$GLMM + rsf_sd_maps$GLM_GRF +
  plot_layout(ncol = 2)
#Different scales, but patterns are quite similar


#NOTE: The GLM and GLMM RSFs were identical to that fitted in frequentist framework w/ {glmmTMB}. But with {inlabru}, it's a little easier to work with the posterior distributions and to extend this to non-linear models. The estimates for the GLM + GRF (model with a spatial random effect) appeared to be quite different from the other two, but the relationships here made sense and also accounted for the spatial autocorrelation in the hourly telemetry data. Therefore, I'd expect the results from this model to be better than the other two. However, this model could be extended further by also accounting for inter-individual variability and possibly non-linear relationships for each of the covariates.

