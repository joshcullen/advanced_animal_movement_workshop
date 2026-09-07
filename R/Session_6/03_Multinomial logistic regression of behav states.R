
### Estimate effect of LULC on probability of behavioral states ###

library(tidyverse)
library(brms)
library(bayesplot)
library(rnaturalearth)
library(sf)
library(terra)
library(tidyterra)

# Enable parallel processing for brms
options(mc.cores = parallel::detectCores())
rstan::rstan_options(auto_write = TRUE)



#################
### Load data ###
#################

# Load tracks (from 3-state HMM w/ multiple imputation and covariates)
dat <- read_csv("processed_data/Session_3/HMM_3state_MultImp_covar.csv") |> 
  mutate(ID = as.character(ID))

summary(dat)
glimpse(dat)


# Load vector layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")


# Load rasters
lulc <- rast("rasters/lulc.tif")
tree_prop <- rast("rasters/tree_prop_50m.tif")
shrub_prop <- rast("rasters/shrub_prop_50m.tif")





#############################
### Extract static covars ###
#############################

# Check that rasters are in same CRS as tracks (UTM Zone 36S; EPSG:32736)
tree_prop; shrub_prop; lulc
# lc prop layers both look good (stored at 10 m res)
# general lulc layer needs to be projected

# Reproject LULC
lulc_proj <- project(lulc, tree_prop, threads = 10, by_util = TRUE)

# Creat SpatRaster stack
lulc_rast <- rast(c(tree = tree_prop, shrub = shrub_prop, lulc = lulc_proj))

# Extract
dat2 <- dat |> 
  cbind(extract(lulc_rast, dat[,c("x","y")], ID = FALSE))



# Explore whether behavioral states ("state_vit") or predictors are missing
summary(dat2)
table(dat2$state_vit, useNA = "ifany")
#6218 points w/ missing LULC


# Explore how many rows have 'tree' and 'shrub' cols that sum to 1
#this will impact parameter identifiability
dat2 |> 
  mutate(lc_sum = tree + shrub) |> 
  filter(lc_sum == 1) |> 
  nrow()
# > 15k rows; looks like we should only keep one of these
#alternatives would be to increase buffer radius during data prep (which would potentially result in fewer probs that sum to 1) or create dummy variables per class




###########################
### Prep data for model ###
###########################

# Data for use in LC proportion models
dat3 <- dat2 |> 
  mutate(
    # Set 'Encamped' as the baseline reference category
    #Use Viterbi estimates (none are unclassified); although you could use the other est. if desired
    state_vit = factor(state_vit, levels = c("Encamped", "Exploratory", "Transit"))#,
    # Ensure discrete LULC is a factor (e.g., Tree, Shrub, Other)
    # lulc_class = factor(lulc_class, levels = c("Other", "Tree", "Shrub"))
  )

# Data for use in LC dummy var models
dat4 <- dat2 |> 
  drop_na(lulc) |> 
  mutate(
  # Set 'Encamped' as the baseline reference category
  #Use Viterbi estimates (none are unclassified); although you could use the other est. if desired
  state_vit = factor(state_vit, levels = c("Encamped", "Exploratory", "Transit")),
  # Ensure discrete LULC is a factor (e.g., Tree, Shrub, Other); treat "Other" as reference
  lulc = factor(lulc, levels = c("(Other)","Tree cover","Shrubland","Grassland","Cropland","Herbaceous wetland"))
)



######################
### Specify priors ###
######################

# Inspect default/suggested priors by {brms}
get_prior(
  state_vit ~ tree,  #this is our model formula
  data = dat3,
  family = categorical(link = "logit")  #multinomial logistic regression
)


# In categorical models, class "b" applies across all non-reference linear predictors 
# (muExploratory and muTransit)
priors <- c(
  # Intercept priors for non-reference states
  prior(normal(0, 1.5), class = "Intercept", dpar = "muExploratory"),
  prior(normal(0, 1.5), class = "Intercept", dpar = "muTransit"),
  
  # Slope priors for proportion tree cover
  prior(normal(0, 1), class = "b", dpar = "muExploratory"),
  prior(normal(0, 1), class = "b", dpar = "muTransit")
)



# Inspect default/suggested priors by {brms}
get_prior(
  state_vit ~ lulc,  #this is our model formula
  data = dat4,
  family = categorical(link = "logit")  #multinomial logistic regression
)
#use same as other models


##################
### Fit models ###
##################

# Let's fit a few different models (either w/ only tree or shrub proportions, or with dummy-coded LULC classes)

### Tree proportion
fit_tree_prop <- brm(
  formula = state_vit ~ tree,
  data    = dat3,
  family  = categorical(link = "logit"),
  prior   = priors,
  chains  = 4,
  cores   = 4,
  iter    = 2000,
  seed    = 2026
)
# took 12 sec

summary(fit_tree_prop)  #R_hat and ESS look good
plot(fit_tree_prop)  #traceplots and histograms look good
mcmc_rank_overlay(fit_tree_prop, regex_pars = "^b_")  #trace rank plots look good




### Shrub proportion
fit_shrub_prop <- brm(
  formula = state_vit ~ shrub,
  data    = dat3,
  family  = categorical(link = "logit"),
  prior   = priors,
  chains  = 4,
  cores   = 4,
  iter    = 2000,
  seed    = 2026
)
# took 12 sec

summary(fit_shrub_prop)  #R_hat and ESS look good
plot(fit_shrub_prop)  #traceplots and histograms look good
mcmc_rank_overlay(fit_shrub_prop, regex_pars = "^b_")  #trace rank plots look good




### Shrub proportion
fit_lulc <- brm(
  formula = state_vit ~ lulc,
  data    = dat4,
  family  = categorical(link = "logit"),
  prior   = priors,
  chains  = 4,
  cores   = 4,
  iter    = 2000,
  seed    = 2026
)
# took 20 sec

summary(fit_lulc)  #R_hat and ESS look good
plot(fit_lulc)  #traceplots and histograms look good
mcmc_rank_overlay(fit_lulc, regex_pars = "^b_")  #trace rank plots look good





###########################################
### Perform posterior predictive checks ###
###########################################

# Tree proportion model
pp_check(fit_tree_prop, ndraws = 100) +
  labs(title = "Tree Prop Posterior Predictive Check", subtitle = "Observed (y) vs. Model Simulated (yrep)")
#Fits very tightly!

# Shrub proportion model
pp_check(fit_shrub_prop, ndraws = 100) +
  labs(title = "Shrub Prop Posterior Predictive Check", subtitle = "Observed (y) vs. Model Simulated (yrep)")
#Fits very tightly!

# LULC dummy var model
pp_check(fit_lulc, ndraws = 100) +
  labs(title = "LULC Posterior Predictive Check", subtitle = "Observed (y) vs. Model Simulated (yrep)")
#Fits very tightly!




########################################
### Viz conditional marginal effects ###
########################################

# Plot tree prop model
plot(conditional_effects(fit_tree_prop, categorical = TRUE, plot = FALSE))[[1]] + 
  labs(
    title = "Behavioral State Probability by Proportion Tree Cover",
    subtitle = "(w/in 50 km buffer)",
    x     = "Proportion Tree Cover",
    y     = "Predicted Probability",
    color  = "Behavioral State",
    fill = "Behavioral State"
  ) +
  theme_minimal(base_size = 14)
#Increased tree cover associated w/ increaing Encamped prob, but decreasing Exploratory and Transit probs
#Encamped has highest prob across all proportions


# Plot shrub prop model
plot(conditional_effects(fit_shrub_prop, categorical = TRUE, plot = FALSE))[[1]] + 
  labs(
    title = "Behavioral State Probability by Proportion Shrub Cover",
    subtitle = "(w/in 50 km buffer)",
    x     = "Proportion Shrub Cover",
    y     = "Predicted Probability",
    color  = "Behavioral State",
    fill = "Behavioral State"
  ) +
  theme_minimal(base_size = 14)
#Increased shrub cover associated w/ decreasing Encamped prob, but increasing Exploratory and Transit probs
#Encamped has highest prob across all proportions, however


# Plot LULC dummy var model
plot(conditional_effects(fit_lulc, categorical = TRUE, plot = FALSE))[[1]] + 
  labs(
    title = "Behavioral State Probability by Land Cover Class",
    subtitle = "'Encamped' is reference state; 'Other' is reference LC class",
    x     = "Land Cover Class",
    y     = "Predicted Probability",
    color  = "Behavioral State",
    fill = "Behavioral State"
  ) +
  ylim(0, 0.65) +
  theme_minimal(base_size = 14)
#Encamped dominates in all but "Herbaceous wetlands" (where Exploratory more likely)





################################
### Make spatial predictions ###
################################

#Rename layer to match model formula
names(tree_prop) <- "tree"


### Subset region for computational purposes (otherwise will take too long and/or crash RStudio session)
# Extract full extent boundaries
e <- ext(tree_prop)

# Calculate center coordinates
x_center <- (e$xmin + e$xmax) / 2
y_center <- (e$ymin + e$ymax) / 2

# Define a 10 km x 10 km bounding box around the center (5 km buffer on each side)
buffer_m <- 5000  #5000 m

crop_extent <- ext(
  x_center - buffer_m, 
  x_center + buffer_m, 
  y_center - buffer_m, 
  y_center + buffer_m
)

# Crop the raster to the central region
tree_prop_sub <- crop(tree_prop, crop_extent)


# Aggregates cells by factor of 10 (just for demonstration; 100 km res now)
tree_toy <- aggregate(tree_prop_sub, fact = 10, fun = mean)



# Custom prediction function that also restricts MCMC draws (provided to terra::predict())
predict_brms_fast <- function(model, data, ...) {
  df_chunk <- as.data.frame(data)
  complete_rows <- complete.cases(df_chunk)
  
  cat_levels <- levels(model$data[[as.character(formula(model)$formula[[2]])]])
  out <- matrix(NA_real_, nrow = nrow(df_chunk), ncol = length(cat_levels))
  colnames(out) <- cat_levels
  
  if (any(complete_rows)) {
    # Capping ndraws avoids memory limits
    prob_array <- fitted(
      model, 
      newdata = df_chunk[complete_rows, , drop = FALSE], 
      re_formula = NA,
      ndraws = 50 
    )
    out[complete_rows, ] <- prob_array[, "Estimate", ]
  }
  return(out)
}




# Predict on toy grid 
state_maps_toy <- predict(tree_toy, fit_multi_prop, fun = predict_brms_fast)


ggplot() +
  geom_spatraster(data = state_maps_toy) +
  scale_fill_viridis_c("Pr(State)") +
  theme_bw(base_size = 14) +
  facet_wrap(~lyr, ncol = 2)
