
### Fit resource selection function (RSF) using GLM and GAM ###

library(tidyverse)
library(sf)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(mgcv)
library(gratia)
library(tictoc)



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


# Predict relative intensity (i.e., relative abundance)
rsf_map_linear <- predict(cov_stack, rsf_linear, type = "response")
rsf_map_gam <- predict(cov_stack, rsf_gam, type = "response")


# Map predictions compared to tracks
rast_preds <- c(rsf_map_linear, rsf_map_gam)
names(rast_preds) <- c("GLM", "GAM")


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


#NOTE: The GAM results are still suspect given that the marginal effects suggested that staying close to human settlements was the strongest effect of the different covars. Marginal effects also suggested that these 3 elephants were more likely to choose areas that were 15-40 km from water (likely related to coarsening of water layer). Also of mention is that these estimates are VERY precise, owing to that fact that it treats all points as independent observations, which is certainly not the case (due to autocorrelation w/in tracks). Also expected to be variability among individuals that will likely drive greater uncertainty and lower effect sizes.

