
### Extract environmental covariates for modeling ###

library(tidyverse)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(sf)
library(tictoc)
library(osmdata)
library(exactextractr)




###################
#### Load data ####
###################

# Load tracks
dat <- read_csv("processed_data/Session_1/cleaned_tracks.csv")

glimpse(dat)
summary(dat)


# Load vector layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")
pop_polygons <- st_read("rasters/pop_polygons.fgb")


# Load rasters
dem <- rast("rasters/dem.tif")
lulc <- rast("rasters/lulc.tif")
dist2pop <- rast("rasters/dist2pop.tif")
ndvi <- rast("rasters/ndvi.nc")




######################################
### Get vector water data from OSM ###
######################################

# Define region of interest
gl_bbox <- st_bbox(c(xmin = min(dat$lon), ymin = min(dat$lat), xmax = max(dat$lon), ymax = max(dat$lat))) + 
  c(-1, -1, 1, 1)

# Query water features w/in ROI
gl_water <- opq(bbox = gl_bbox) |>
  add_osm_feature(key = "water", key_exact = FALSE, value_exact = FALSE) |>
  osmdata_sf()

# Separate by class type
water_multipolygons <- gl_water$osm_multipolygons
water_polygons <- gl_water$osm_polygons

# Merge OSM layers together
water_polygons2 <- water_multipolygons[["geometry"]] |> 
  st_cast("POLYGON")
water_osm <- c(st_geometry(water_polygons), st_geometry(water_polygons2))


# Viz water
ggplot() +
  geom_sf(data = africa, color = "black", linewidth = 1) +
  geom_sf(data = water_osm, color = 'grey30', fill = "lightblue2") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))





####################################
### Calculate derived covariates ###
####################################

### From DEM
slope <- terrain(dem, v = "slope")  #slope (in degrees)
aspect <- terrain(dem, v = "aspect")  #aspect (in degrees)
tri <- terrain(dem, v = "TRI")  #terrain ruggedness index (TRI)


### From OSM water

## Calc distance to water

# Reproject each layer to UTM Zone 36S
dat_sf <- dat |> 
  st_as_sf(coords = c('lon','lat'), crs = 4326, remove = FALSE) |> 
  st_transform(32736)

water_osm_proj <- water_osm |> 
  st_transform(32736)

# Find index of nearest water feature for each GPS point
nearest_idx <- st_nearest_feature(dat_sf, water_osm_proj)

# Calc distance (in meters) to water and assign to relocations
dat_sf$dist2water <- as.numeric(  #convert to numeric value
  st_distance(dat_sf, water_osm_proj[nearest_idx, ], by_element = TRUE)  #calc distance
  )

# Viz density distrib of distances
ggplot(dat_sf) +
  geom_density(aes(dist2water / 1000), fill = "lightblue3") +  #show in km!
  theme_minimal()
#looks like they're often close (< 20 km), but sometimes pretty far (> 40 km)




### From LULC

# Print IDs and LULC classes
levels(lulc)

ggplot() +
  geom_spatraster(data = lulc, use_coltab = TRUE) +
  scale_fill_coltab(name = "LULC", data = lulc) +
  geom_sf(data = africa, fill = NA, linewidth = 0.5, color = "black") +
  geom_path(data = dat, aes(lon, lat, group = id), color = "black", linewidth = 0.15, alpha = 0.5) +
  theme_bw() +
  coord_sf(xlim = ext(lulc)[1:2],
           ylim = ext(lulc)[3:4],
           expand = FALSE)
#seems like 'Tree cover' and 'Shrubland' are often used
#but let's explore this a bit closer via GeoLibre (for interactivity)

plotly::ggplotly(
  ggplot() +
    geom_spatraster(data = lulc, use_coltab = TRUE) +
    scale_fill_coltab(name = "LULC", data = lulc) +
    geom_sf(data = africa, fill = NA, linewidth = 0.5, color = "black") +
    geom_path(data = dat, aes(lon, lat, group = id), color = "black", linewidth = 0.15, alpha = 0.5) +
    theme_bw() +
    coord_sf(xlim = ext(lulc)[1:2],
             ylim = ext(lulc)[3:4],
             expand = FALSE)
)


## Generate stand-alone layers per each LULC class

# Water (including permanent water and herbaceous wetlands)
water_lc <- lulc  #copy lulc layer
water_lc[!water_lc %in% c(80,90)] <- NA  #convert all non-water to NA
plot(water_lc)

# Tree cover
tree_lc <- lulc
tree_lc[tree_lc != 10] <- NA
plot(tree_lc)

# Shrub cover
shrub_lc <- lulc
shrub_lc[shrub_lc != 20] <- NA
plot(shrub_lc)


# Aggregate (i.e., coarsen) from 10m to 100m (fact = 10) to prevent crashing computer
# 'max' ensures that if ANY 10m sub-pixel was water, the 100m cell remains water
# water_lc_100m <- aggregate(water_lc, fact = 10, fun = "max", na.rm = TRUE)
# tree_lc_100m <- aggregate(tree_lc, fact = 10, fun = "max", na.rm = TRUE)
# shrub_lc_100m <- aggregate(shrub_lc, fact = 10, fun = "max", na.rm = TRUE)





#############################################################
### Merge OSM and LULC water data together (data fusion)  ###
#############################################################

# Viz each layer
par(mfrow = c(1,2))
plot(water_osm, col = "blue", border = NA, main = "OSM")

plot(water_lc, main = "LULC")
par(mfrow = c(1,1))



# Vectorize the LULC water pixels
# 'values = FALSE' ignores cell values and vectorizes only non-NA pixels
tic()
lulc_water_poly <- as.polygons(water_lc, values = FALSE, dissolve = TRUE)
toc()  #took 35 sec

# Convert from SpatVector to sf object
lulc_water_sf <- st_as_sf(lulc_water_poly) |> 
  st_cast("POLYGON")

# Filter out tiny noise patches (e.g., single 10m pixels or artifacts < 500 m²)
water_area_m2 <- as.numeric(st_area(lulc_water_sf))
lulc_water_sf_sub <- lulc_water_sf[water_area_m2 >= 500,]

# Combine LULC geometries with your existing OSM vector geometries
# c(st_geometry(...)) strips attribute metadata and merges spatial shapes cleanly
water <- c(st_geometry(lulc_water_sf_sub), water_osm)


# Dissolve all polygons together (to remove overlapping duplicates) and then break apart again
tic()
water_sub <- water |> 
  st_union() |> 
  st_cast("POLYGON")
toc()  #took 1 min

plot(water_sub, col = "blue", border = NA, main = "OSM + LULC")


# Calc new distance to water measure
water_proj <- st_transform(water_sub, 32736)

tic()
nearest_idx <- st_nearest_feature(dat_sf, water_proj)
dat_sf$dist2water_2 <- as.numeric(
  st_distance(dat_sf, water_proj[nearest_idx, ], by_element = TRUE)
)
toc()  #took 1.5 min


# Compare density distrib of distances between data sources (OSM vs combined)
ggplot(dat_sf) +
  geom_density(aes(dist2water / 1000, fill = "OSM"), alpha = 0.4) +  #show in km!
  geom_density(aes(dist2water_2 / 1000, fill = "OSM + LULC"), alpha = 0.4) +  #show in km!
  scale_fill_manual("Data Source", values = c("lightblue", "darkblue")) +
  labs(x = "Distance to water (km)", y = "Density") +
  theme_minimal()
#We see that merging both datasets together produces distance measures that are much closer to water compared to using OSM alone
#Both are static layers of water and therefore are a simplified representation (i.e., snapshot) of water on the landscape, but likely represent locations with at least seasonal if not semi-permanent water over time







#########################################
### Extract covariates at exact point ###
#########################################

### Static layers ###

# Combine all DEM-related products
dem_prods <- c(dem, slope, aspect, tri)
names(dem_prods)[1] <- "elev"
plot(dem_prods)


# Extract from a SpatRaster stack
dat_sf2 <- dat_sf |> 
  cbind(extract(dem_prods, dat_sf, ID = FALSE))  #auto-projects our 'sf' object if different CRS from raster

# Extract from a single SpatRaster layer
dat_sf2$dist2pop <- extract(dist2pop, dat_sf2, ID = FALSE) |> 
  unlist() |> 
  as.numeric()

dat_sf2$lulc <- extract(lulc, dat_sf2, ID = FALSE) |> 
  unlist()





### Dynamic layer ###

#While more common in marine studies, dynamic variables (e.g., vegetation indices, water extent, burned area, duman development, etc) play major roles on when and where animals move. While some regions that are relatively unchanging and stable over years and seasons may be well-suited to use of a static LULC layer, there are often times that dynamic layers are necessary for ecological analyses

#Since the layers we want to extract vary over space for a specific time (or time window), we need to extract time-matched values in space

#This is more straightforward if done at a monthly scale or daily scale, but for products that have repeat satellite visits every 8, 10, or even 16 days, we need to take a slightly more complicated approach to time-matched extraction of covariates since very few points will likely fall on these specific dates of Earth observation (assuming the products aren't some composite over that window)

# The 16-day MODIS NDVI product we downloaded has dates now stored for the midpoint of each time window. So we will match up locations based on the closest date and then extract them for those particular layers


# Ensure track timestamps and raster times are comparable both "Date" vectors
track_dates <- as_date(dat_sf2$date) 
ndvi_dates <- as_date(time(ndvi))

# Match each GPS point to the index of the nearest NDVI layer
dat_sf2$ndvi_lyr_idx <- sapply(track_dates, function(d) {
  which.min(abs(ndvi_dates - d))
})

# Initialize an empty column to hold the extracted values
dat_sf2$ndvi <- NA

# Extract data iteratively (by time slice)
# We only loop through layers that actually have corresponding tracking points
unique_layers <- unique(dat_sf2$ndvi_lyr_idx)


tic()
for (i in unique_layers) {
  # Find all track points associated with this specific time slice
  pt_idx <- which(dat_sf2$ndvi_lyr_idx == i)

  # Extract values from just this single raster layer for these specific points
  # ID = FALSE prevents terra from returning a column of ID numbers
  extracted_vals <- extract(
    ndvi[[i]],
    dat_sf2[pt_idx, ],
    ID = FALSE
  )

  # Assign the extracted values back to the correct rows in the sf object
  dat_sf2$ndvi[pt_idx] <- extracted_vals[[1]]
  message(paste0("Extracted NDVI for Layer ", i))
}
toc()  # took 1.5 sec

# Clean up the temporary index column using dplyr
dat_sf2 <- dat_sf2 |> 
  select(-ndvi_lyr_idx)




##########################################
### Extract covariates within a buffer ###
##########################################

#this may be preferred when accounting for tag location error and/or an animal's perceived habitat in its vicinity


### Extract proportion cover of Trees and Shrubland w/in 50 m of location

# Convert to binary layers
tree_bin <- ifel(!is.na(tree_lc), 1, 0)
shrub_bin <- ifel(!is.na(shrub_lc), 1, 0)

# Stack layers
lc_stack <- c(tree_bin, shrub_bin) |> 
  project(crs(dat_sf2), res = 10, method = "near", threads = 10, use_gdal = TRUE)  #reproject to match tracks
names(lc_stack) <- c('tree','shrub')
plot(lc_stack)
lc_stack
  
# Create polygons of buffered points
pts_buff <- dat_sf2 |> 
  st_buffer(dist = 50)  #units already in meters for EPSG:32736

# Use function from {exactextractr} for fast extraction of rasters by polygons
props <- exact_extract(lc_stack, pts_buff, fun = "mean")  #mean gives us proportion coverage per class
names(props) <- c("prop_tree","prop_shrub")

dat_sf2 <- cbind(dat_sf2, props)



### Calculate spatial layers of 'tree' and 'shrub' proportions using ~50 m buffer distance

# Calc focal props (each cell is 10 m, so we need 11 cells to cover full diameter to equate 50 m buffer)

tic()
tree_prop <- focal(lc_stack[["tree"]], w = 11, fun = "mean", na.rm = TRUE, cores = 10)
toc()  #took 20 sec

tic()
shrub_prop <- focal(lc_stack[["shrub"]], w = 11, fun = "mean", na.rm = TRUE, cores = 10)
toc()  #took 20 sec




##############
### Export ###
##############

# Save new hybrid water vector layer (OSM + LULC)
st_write(water_sub, "rasters/water.fgb")

# Save layers for proportion of 'tree' and 'shrub' cover w/in 50 m buffer
writeRaster(tree_prop, "rasters/tree_prop_50m.tif", overwrite = TRUE)  #1.7 GB
writeRaster(shrub_prop, "rasters/shrub_prop_50m.tif", overwrite = TRUE)  #1.7 GB

# Convert from sf object to data.frame (including UTM coords; x,y)
dat_out <- dat_sf2 |> 
  mutate(x = st_coordinates(geometry)[,"X"],
         y = st_coordinates(geometry)[,"Y"],
         .after = lon) |> 
  st_drop_geometry()


write_csv(dat_out, "processed_data/Session_4/tracks_covars.csv")
