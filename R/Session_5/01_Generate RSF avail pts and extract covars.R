
### Generate "available" points and extract environ covars ###

library(tidyverse)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(sf)
library(tictoc)

source("R/utils.R")



###################
#### Load data ####
###################

# Load tracks
dat <- read_csv("processed_data/Session_1/cleaned_tracks.csv") |> 
  mutate(id = as.character(id),
         obs = 1)  #denotes "used" habitat for model

glimpse(dat)
summary(dat)


# Load vector layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")
pop_polygons <- st_read("rasters/pop_polygons.fgb")


# Load rasters
# dem <- rast("rasters/dem.tif")
lulc <- rast("rasters/lulc.tif")
dist2pop <- rast("rasters/dist2pop.tif")
ndvi <- rast("rasters/ndvi.nc")
water <- st_read("rasters/water.fgb")
tree_prop <- rast("rasters/tree_prop_50m.tif")
shrub_prop <- rast("rasters/shrub_prop_50m.tif")


# Viz spatial extent of tracks
ggplot() +
  geom_sf(data = africa) +
  geom_path(data = dat, aes(lon, lat, group = id), color = "black", linewidth = 0.15, alpha = 0.5) +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))
#while we could use an MCP or KDE estimate to define where to sample habitat, let's just stick with a bounding box





#############################################
### Generate background/quadrature points ###
#############################################

# While many studies treat these points as a representation of the habitat "available" to an animal, it's better to think of them as a numerical integration technique
# So essentially we're trying to define all the different types of habitat in the vicinity of the observed points to better estimate the species-habitat associations. And the more points (or "better" selection of background points) the greater the estimated coefficient precision by the model

# In this example, we'll just use a single set of background points at a 10:1 (available:used) ratio. But it's highly suggested that a few different ratios are selected (e.g., 10:1, 30:1, 50:1); once the estimated coefficients have stabilized, that means that the model has more-or-less reached the asymptotic estimate


# Define bounding box for tracks
bbox <- st_bbox(c(xmin = min(dat$lon), ymin = min(dat$lat), xmax = max(dat$lon),  #using range of coords
                          ymax = max(dat$lat)) + c(-0.5, -0.5, 0.5, 0.5), crs = 4326) |>  #buffered by 0.5° on each side
  st_as_sfc() |> 
  st_as_sf()


# Viz tracks compared to new bbox
ggplot() +
  geom_sf(data = africa) +
  geom_point(data = dat, aes(lon, lat), color = "black", size = 0.15, alpha = 0.5) +
  geom_sf(data = bbox, color = "blue", fill = NA) +
  theme_bw() +
  coord_sf(xlim = st_bbox(bbox)[c("xmin","xmax")],
           ylim = st_bbox(bbox)[c("ymin","ymax")])



# Generate available points and randomly assign to date per ID (for dynamic NDVI layer)
# ratio of 1:10 (used:available)

n_avail <- 10 * nrow(dat)

set.seed(2026)
id_idx <- rep(1:nrow(dat), 10)

avail_pts <- st_sample(bbox, size = n_avail) |>  #random background sampling from bbox
  data.frame() |>
  mutate(obs_id = sample(id_idx),  #randomly sample a row from 'dat' (10 times)
         id = dat$id[obs_id],  #select ID based on sampled row order from 'dat'
         date = dat$date[obs_id],  #pull date associated with observation; this could be randomly sampled, however
         obs = 0,  #denotes 'available' points
         lon = unlist(map(geometry, 1)),
         lat = unlist(map(geometry, 2)),
         .before = geometry) |>
  select(-c(geometry, obs_id))


# Viz available point spread by month
ggplot() +
  geom_point(data = avail_pts, aes(lon, lat, color = month(date)), size = 0.05, alpha = 0.5) +
  geom_point(data = dat, aes(lon, lat), size = 0.2) +
  geom_sf(data = bbox, color = "red", fill = NA, linewidth = 0.5) +
  theme_bw() +
  facet_wrap(~ month(date)) +
  coord_sf(xlim = st_bbox(bbox)[c("xmin","xmax")],
           ylim = st_bbox(bbox)[c("ymin","ymax")])

# Viz example of available point for ID 5605 (strictly related to dates used)
ggplot() +
  geom_point(data = avail_pts |> 
               filter(id == 5605), aes(lon, lat), size = 0.05, alpha = 0.5) +
  geom_point(data = dat, aes(lon, lat, color = month(date)), size = 0.3) +
  scale_color_viridis_c("Month") +
  geom_sf(data = bbox, color = "red", fill = NA, linewidth = 0.5) +
  theme_bw() +
  coord_sf(xlim = st_bbox(bbox)[c(1,3)],
           ylim = st_bbox(bbox)[c(2,4)])



### Combine used and available pts
dat2 <- rbind(dat |> select(-temp),
              avail_pts |> tibble() |> select(id, date, lat, lon, obs))




##########################################################
### Check that covariates have same CRS and resolution ###
##########################################################

#Before any extraction can be done, we need to make sure all layers have the same CRS so that we're getting correct values
#We also need to make sure all rasters are at the same spatial resolution. This is particularly important if wanting to generate mapped predictions of the fitted RSF, but still important if not b/c it's difficult to compare an effect for one variable at a fine resolution that possibly had high variability within a single cell of a coarser resolution raster variable


### Make sure all layers are in UTM Zone 36S proj

dist2pop  #WGS84
ndvi  #WGS84
shrub_prop  #UTM
tree_prop  #UTM
water  #WGS84

# Project NDVI and dist2pop
ndvi_proj <- project(ndvi, crs(tree_prop), threads = 10, use_gdal = TRUE)
dist2pop_proj <- project(dist2pop, crs(tree_prop), threads = 10, use_gdal = TRUE)



### Check that all rasters are at coarsest (i.e., limiting) variable res
res(dist2pop_proj)  #439 m
res(ndvi_proj)  #257 m
res(tree_prop)  #10 m
res(shrub_prop)  #10 m
# resample all to same res as 'dist2pop'


# Resample rasters to coarser resolution
ndvi_proj2 <- resample(ndvi_proj, dist2pop_proj, threads = 10, by_util = TRUE)
tree_prop2 <- resample(tree_prop, dist2pop_proj, threads = 10, by_util = TRUE)
shrub_prop2 <- resample(shrub_prop, dist2pop_proj, threads = 10, by_util = TRUE)




### Create dist2water raster

# Rasterize the water layer (and convert to UTM)
tic()
water_binary <- rasterize(
  x = st_transform(water, 32736), 
  y = tree_prop,          # rasterizes over 10 m grid
  field = 1,              # Assign all water polygons a value of 1
  touches = TRUE,         # ANY physical overlap turns the pixel to 1
  background = NA,        # Set everything else to NA
  fun = "max"             # If multiple water polygons overlap a cell, max(1,1) = 1
)
toc()  #took 38 sec
plot(water_binary, col = "blue")

# Coarsen water raster to match other rasters
water_binary2 <- resample(water_binary, dist2pop_proj, method = "near", threads = 10, by_util = TRUE)
plot(water_binary2, col = "blue")

# Compute distance to water raster
#this approach is much more computationally efficient for large datasets (as opposed to calculating nearest distance to vector layer for >100k points)
dist2water <- distance(water_binary2)




### Create stack of all static rasters

static_covars <- c(dist2pop_proj, tree_prop2, shrub_prop2, dist2water)
names(static_covars) <- c("dist2pop", "tree_prop", "shrub_prop", "dist2water")
plot(static_covars)




##########################
### Extract covariates ###
##########################

# Add UTM coords
dat3 <- dat2 |> 
  add_trans_coords(coords = c('lon','lat'), proj = 4326, new_proj = 32736)


### Extract static covars ###
dat4 <- dat3 |> 
  cbind(extract(static_covars, dat3[,c("x","y")], ID = FALSE))


# Need to create "Other" LULC class that catches everything else (these all need to sum to 1)
dat4 <- dat4 |> 
  mutate(other_prop = 1 - c(tree_prop + shrub_prop),
         .after = shrub_prop)



### Extract time-matched dynamic covars ###

# Ensure track timestamps and raster times are comparable both "Date" vectors
track_dates <- as_date(dat4$date) 
ndvi_dates <- as_date(time(ndvi_proj2))

# Match each GPS point to the index of the nearest NDVI layer
dat4$ndvi_lyr_idx <- sapply(track_dates, function(d) {
  which.min(abs(ndvi_dates - d))
})

# Initialize an empty column to hold the extracted values
dat4$ndvi <- NA

# Extract data iteratively (by time slice)
# We only loop through layers that actually have corresponding tracking points
unique_layers <- unique(dat4$ndvi_lyr_idx)


tic()
for (i in unique_layers) {
  # Find all track points associated with this specific time slice
  pt_idx <- which(dat4$ndvi_lyr_idx == i)
  
  # Extract values from just this single raster layer for these specific points
  # ID = FALSE prevents terra from returning a column of ID numbers
  extracted_vals <- extract(
    ndvi_proj2[[i]],
    dat4[pt_idx, c('x','y')],
    ID = FALSE
  )
  
  # Assign the extracted values back to the correct rows in the sf object
  dat4$ndvi[pt_idx] <- extracted_vals[[1]]
  message(paste0("Extracted NDVI for Layer ", i))
}
toc()  # took 1.5 min

# Clean up the temporary index column using dplyr
dat4 <- dat4 |> 
  select(-ndvi_lyr_idx)




########################################################
### Summarize and explore patterns in extracted data ###
########################################################

summary(dat4)
#NDVI is only covar w/ NAs (4348)
#possible to do some imputation, but we'll just remove these records


# Viz density distributions per covar (separately between 'used' and 'available')
dat4 |> 
  pivot_longer(cols = dist2pop:ndvi, names_to = "covar", values_to = "value") |> 
  mutate(presabs = ifelse(obs == 0, "Available", "Used")) |> 
  
  ggplot() +
  geom_density(aes(value, color = presabs)) +
  geom_rug(aes(value, color = presabs)) +
  scale_color_brewer("", palette = "Set1") +
  labs(x = "Value", y = "Density") +
  theme_bw() +
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text = element_text(size = 10, face = "bold"),
        plot.title = element_text(size = 14, face = "bold"),
        legend.position = "top") +
  facet_wrap(~ covar, scales = "free")
#Looks like we can see separation between the 'used' and 'available' distributions for just about all covars
#This means that these covariates should be informative during model fitting
#Also, there don't appear to be any outliers, which is good; in some cases, it makes sense to filter out extreme values that would otherwise negatively effect the ability of the model to estimate the maximum likelihood (or find the posterior distribution)





###########################################
### Export data and transformed rasters ###
###########################################

write_csv(dat4, "processed_data/Session_5/dat_presabs_covars.csv")


writeRaster(dist2pop_proj, "rasters/dist2pop_rsf.tif")
writeRaster(dist2water, "rasters/dist2water_rsf.tif")
writeRaster(tree_prop2, "rasters/tree_prop_rsf.tif")
writeRaster(shrub_prop2, "rasters/shrub_prop_rsf.tif")
writeCDF(ndvi_proj2, "rasters/ndvi_rsf.nc", varname = "NDVI", timename = "time", overwrite = TRUE)  #700 MB
