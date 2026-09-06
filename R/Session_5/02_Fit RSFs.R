
### Fit resource selection function (RSF) using GLM and GAM ###

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(mgcv)
library(gratia)
library(tictoc)
library(glmmTMB)

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

# Logit link logistic regression
rsf_linear <- glm(obs ~ dist2pop_s + dist2water_s + ndvi_s,
  data = dat2,
  family = binomial(link = "logit"),
  weights = wts
)

summary(rsf_linear)
# For IWLR, we ignore the intercept since this is directly tied to N_avail


# Calculate Relative Selection Strengths (Exponential of Coefficients)
# Coefficients > 1 indicate positive selection; < 1 indicate avoidance
selection_ratios <- exp(coef(rsf_linear)[-1])
selection_ratios
#so we see some slight positive selection for higher NDVI and distances greater from water, and a negative relationship w/ distance to human settlements (i.e., closer to humans than expected)
#however, it's hard to fully understand these patterns (and the estimated uncertainty) from these numbers alone; let's plot these marginal effects


### Predict marginal effects
covars <- c("dist2pop_s", "dist2water_s", "ndvi_s")

rsf_preds <- predict_rsf_margeff(fit = rsf_linear, focal_covars = covars, intercept = TRUE, data  = dat2)

# Plot all marginal effects
ggplot(rsf_preds, aes(x = x_natural, y = rss)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = rss_lwr, ymax = rss_upr), fill = "#2c7fb8", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "Marginal Effects of RSF Covariates",
    subtitle = "Relative Selection Strength (RSS) with other covariates held at their mean",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.placement = "outside", strip.text = element_text(face = "bold"))





###########################
### Fit RSF as IWLR GAM ###
###########################

#we'll use a set of 5 knots (k) for defining the relative complexity of our smoothed curves
#we'll also use cubic regression splines instead of the default thin plate splines

tic()
rsf_gam <- gam(obs ~ s(dist2pop_s, k = 5, bs = "cr") + s(dist2water_s, k = 5, bs = "cr") + 
    s(ndvi_s, k = 5, bs = "cr"),
  data = dat2,
  family = binomial(link = "logit"),
  weights = wts,
  method = "REML"
)
toc()  #took 14 sec

summary(rsf_gam)
#estimated degrees of freedom (edf) are all ~4, so we have nonlinear relationships
#all p-values for smooth terms are significant, so they represent *statistically significant* relationships
#low R^2 and Dev. explained, so probably not the most informative model

# Check goodness-of-fit (although typically harder for logistic regression)
gratia::appraise(rsf_gam)



### Visualize Partial Non-linear Effects

# base R plotting
plot(rsf_gam, pages = 1, all.terms = TRUE)

#nicer functions for ggplot2 from {gratia}
draw(rsf_gam) &
  geom_hline(yintercept = 0, linewidth = 0.5, linetype = "dashed") &
  theme_bw()

# Custom marginal effects plot
gam_preds <- predict_gam_margeff(fit = rsf_gam, data = dat2)

ggplot(gam_preds, aes(x = x_natural, y = rss)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = rss_lwr, ymax = rss_upr), fill = "#2c7fb8", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "Non-Linear Marginal Effects of RSF Covariates (GAM)",
    subtitle = "Relative Selection Strength (RSS) centered on average availability",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.placement = "outside", strip.text = element_text(face = "bold"))
#relationships are highly non-linear; indicating that GLM would likely obscure true relationship
#these plots suggest that (when other covars held at their mean value) 1) elephants are more likely to select habitat 1-30 km from human settlements, 2) elephants prefer habitat 10-40 km from water sources, and 3) elephants are more likely to select habitat where vegetation has a moderately high NDVI value (0.5-0.75)




#########################################
### Fit RSF as mixed-effects IWLR GLM ###
#########################################

#This example follows recommendations from Muff et al. (2020) "Accounting for individual-specific variation in habitat-selection studies: Efficient estimation of mixed-effects models using Bayesian or frequentist computation" (https://doi.org/10.1111/1365-2656.13087)

#Essentially, we'll be fitting the RSF GLM again, but now allowing the intercept AND slope to vary by individual. This is much better practice than applying random intercepts alone. Another option is the 2-step approach recommended by John Fieberg (e.g., Fieberg et al. 2009; https://doi.org/10.1111/j.1365-2664.2009.01692.x) where each individual has a separate RSF fitted, and then these are combined in a principled way to achieve population inference; this is a good option if not wanting to rely on assumptions typically used for random effects.

# This example implements this hierarchical approach in {glmmTMB}, which will slightly differ from using INLA/inlabru or a purely Bayesian model


# Set up (but don't yet fit) model
rsf_linear_h.tmp <- glmmTMB(obs ~ dist2pop_s + dist2water_s + ndvi_s +  #fixed effects
                          (1|id) + (0 + dist2pop_s|id) + (0 + dist2water_s|id) + (0 + ndvi_s|id),  #varying effects
                        family = binomial(link = "logit"),
                        data = dat2, 
                        doFit = FALSE,
                        weights = wts)

# Fix SD of first random term (`(1|id)`; varying intercept) to 1e3 (i.e., 1000), which corresponds to variance of 1e6
#must be on log scale
rsf_linear_h.tmp$parameters$theta[1] <- log(1e3)

# Tell glmmTMB to leave the first param "theta[1]" as fixed, but estimate all others
rsf_linear_h.tmp$mapArg <- list(theta = factor(c(NA, 1:3)))  #vector must be length of fixed effects


### Fit the hierarchical model ###
tic()
rsf_linear_h <- fitTMB(rsf_linear_h.tmp)
toc()  #took 45 sec

summary(rsf_linear_h)
#results show some general differences compared to simple RSF




### Extract Population-Level & Individual Selection Ratios ###

# Extract pop-level fixed effects
pop_coefs <- fixef(rsf_linear_h)$cond[-1]  #exclude intercept

# Calculate population relative selection strengths (excluding intercept)
selection_ratios_pop <- exp(pop_coefs)
selection_ratios_pop
#relatively similar to simple linear RSF, but effect of dist2pop is actually stronger

# Extract individual-level coefficients (Fixed + Random Slopes per ID)
indiv_coefs <- coef(rsf_linear_h)$cond$id[,-1]  #remove intercept
indiv_selection_ratios <- exp(indiv_coefs)
indiv_selection_ratios
#plenty of inter-individual variability per covar


### Predict marginal effects
covars <- c("dist2pop_s", "dist2water_s", "ndvi_s")

# Predict population-level marginal effects for RSS w/ 95% CI 
rsf_preds <- predict_rsf_margeff(fit = rsf_linear_h, focal_covars = covars, intercept = TRUE, data = dat2,
                                 level = "population")

# Plot all marginal effects
ggplot(rsf_preds, aes(x = x_natural, y = rss)) +
  geom_line(color = "#2c7fb8", linewidth = 1) +
  geom_ribbon(aes(ymin = rss_lwr, ymax = rss_upr), fill = "#2c7fb8", alpha = 0.2) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "Pop.-Level Marginal Effects of Covariates for Mixed-Effects RSF",
    subtitle = "Relative Selection Strength (RSS) with other covariates held at their mean",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.placement = "outside", strip.text = element_text(face = "bold"))




# Predict individual-level marginal effects for RSS w/ 95% CI 
rsf_preds_id <- predict_rsf_margeff(fit = rsf_linear_h, focal_covars = covars, intercept = TRUE, data = dat2,
                                 level = "individual")

# Plot all marginal effects
ggplot(rsf_preds_id, aes(x = x_natural, y = rss, group = id, color = factor(id))) +
  geom_line(linewidth = 1) +
  scale_color_brewer("ID", palette = "Set1") +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkred") +
  facet_wrap(~ covariate, scales = "free_x", strip.position = "bottom") +
  labs(
    title = "ID-Level Marginal Effects of Covariates for Mixed-Effects RSF",
    subtitle = "Relative Selection Strength (RSS) with other covariates held at their mean",
    y = "Relative Selection Strength (RSS)",
    x = NULL
  ) +
  theme_minimal(base_size = 14) +
  theme(strip.placement = "outside", strip.text = element_text(face = "bold"))






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

# Construct raster stack matching covariate names (w/ mean NDVI)
#likewise, time-varying values of NDVI could be used to make dynamic predictions of habitat selection
cov_stack <- c(dist2water_s, dist2pop_s, ndvi_s_mean)
names(cov_stack) <- c("dist2water_s", "dist2pop_s", "ndvi_s")


# Predict log(relative intensity) (i.e., log relative abundance)
#Also generally a good idea to show the uncertainty in predictions too, although this should be done w/ intercept removed
rsf_map_linear <- predict(cov_stack, rsf_linear, fun = predict_rsf_raster, type = "all")
rsf_map_gam <- predict(cov_stack, rsf_gam, fun = predict_gam_rsf, type = "all")
rsf_map_linear_h <- predict(cov_stack, rsf_linear_h, fun = predict_rsf_raster, type = "all") 

# Exponentiate log-scale models


# Map predictions compared to tracks
rast_preds <- c(rsf_map_linear$rss, rsf_map_gam$rss, rsf_map_linear_h$rss)
names(rast_preds) <- c("GLM", "GAM", "GLMM")

rast_pred_se <- c(rsf_map_linear$se_link, rsf_map_gam$se_link, rsf_map_linear_h$se_link)
names(rast_pred_se) <- c("GLM", "GAM", "GLMM")


# Viz mapped predictions of RSS (mean)
ggplot() +
  geom_spatraster(data = rast_preds) +
  scale_fill_viridis_c("Relative Intensity", na.value = "transparent") +
  geom_sf(data = africa, fill = NA, color = "white", lwd = 1) +
  geom_path(data = dat2 |> 
              filter(obs == 1),  #only keep the observed locs
            aes(x, y, group = id), color = "white", lwd = 0.25, alpha = 0.5) +
  theme_bw(base_size = 14) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top") +
  guides(fill = guide_colorbar(
    title.position = "top",         # Moves title above the bar for more space
    barwidth = unit(10, "cm"),      # Lengthens the colorbar (adjust as needed)
    barheight = unit(0.5, "cm")     # Adjusts the thickness
  )) +
  labs(x = "Easting", y = "Northing") +
  coord_sf(xlim = ext(rast_preds)[1:2],
           ylim = ext(rast_preds)[3:4],
           expand = FALSE) +
  facet_wrap(~ lyr)
#GLM seems to be an oversimplification, whereas GAM seems to predict these tracks better


# Viz mapped prediction uncertainty of RSS (SE)
ggplot() +
  geom_spatraster(data = rast_pred_se) +
  scale_fill_viridis_c("SD of Relative Intensity", option = "rocket", na.value = "transparent") +
  geom_sf(data = africa, fill = NA, color = "white", lwd = 1) +
  geom_path(data = dat2 |> 
              filter(obs == 1),  #only keep the observed locs
            aes(x, y, group = id), color = "white", lwd = 0.25, alpha = 0.5) +
  theme_bw(base_size = 14) +
  theme(strip.text = element_text(face = "bold"),
        legend.position = "top") +
  guides(fill = guide_colorbar(
    title.position = "top",         # Moves title above the bar for more space
    barwidth = unit(10, "cm"),      # Lengthens the colorbar (adjust as needed)
    barheight = unit(0.5, "cm")     # Adjusts the thickness
  )) +
  labs(x = "Easting", y = "Northing") +
  coord_sf(xlim = ext(rast_preds)[1:2],
           ylim = ext(rast_preds)[3:4],
           expand = FALSE) +
  facet_wrap(~ lyr)
#GLMM shows much greater uncertainty since it actually accounts for ID variability


#NOTE: The GAM results are still suspect given that the marginal effects suggested that staying close to human settlements was the strongest effect of the different covars. Marginal effects also suggested that these 3 elephants were more likely to choose areas that were 15-40 km from water (likely related to coarsening of water layer). Also of mention is that these estimates are VERY precise, owing to that fact that it treats all points as independent observations, which is certainly not the case (due to autocorrelation w/in tracks). Also expected to be variability among individuals that will likely drive greater uncertainty and lower effect sizes, which was somewhat accounted for in the mixed-effects model. While not shown here, it's also possible to account for varying slopes in GAMs as well.

