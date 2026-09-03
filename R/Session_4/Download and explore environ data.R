
### Explore and process remote sensing data ###

library(tidyverse)
library(terra)
library(tidyterra)
library(rnaturalearth)
library(sf)
library(tictoc)
library(rsi)
library(elevatr)
library(osmdata)
library(rstac)
library(ncdf4)

## Set GDAL settings for optimized cloud data streaming

# Prevent GDAL from scanning remote directories
setGDALconfig("GDAL_DISABLE_READDIR_ON_OPEN", "EMPTY_DIR")
# Bundle requests
setGDALconfig("GDAL_HTTP_MERGE_CONSECUTIVE_READS", "YES")
setGDALconfig("GDAL_HTTP_MULTIPLEX", "YES")  
setGDALconfig("GDAL_HTTP_VERSION", "2")
# Enable caching for the Virtual File System
setGDALconfig("VSI_CACHE", "TRUE")  



###################
#### Load data ####
###################

# Load tracks
dat <- read_csv("processed_data/cleaned_tracks.csv")

glimpse(dat)
summary(dat)


# Load spatial layers
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")
gl_pa <- st_read("raw_data/gltfca_protectedAreasDetailed.shp")





############################################
### Explore mapping political boundaries ###
############################################

# {rgeoboundaries} another good resource, especially for lower level admin bounds
#remotes::install_github("wmgeolab/rgeoboundaries")

# Full countries
ggplot() +
  geom_sf(data = africa |> 
            filter(name %in% c('South Africa','Mozambique'))) +
  # geom_sf_text(data = africa |> 
  #                filter(name %in% c('South Africa','Mozambique')),
  #              aes(label = name)) +  #add labels
  theme_bw() +
  coord_sf(xlim = c(15,42),
           ylim = c(-35,-10))


# Lower-level admin bounds
za_mz <- ne_states(country = c('South Africa','Mozambique'), returnclass = "sf")

ggplot() +
  geom_sf(data = za_mz) +
  # geom_sf_text(data = za_mz, aes(label = name), size = 3) +  #add labels
  theme_bw() +
  coord_sf(xlim = c(15,42),
           ylim = c(-35,-10))




########################
### Process DEM data ###
########################

# Extract elevation at points only
tic()
dat2 <- dat |> 
  rename(x = lon, y = lat) |>  #needs to have cols named 'x' and 'y'
  data.frame() |>  #can't be a tibble
  get_elev_point(prj = "epsg:4326", src = "aws", ncpu = 5, z = 10)
toc()  #took 16 sec

summary(dat2$elevation)  #elevation ranges from 23 to 367 m asl


# Download DEM layer based on extent of locs
tic()
dem <- dat |> 
  rename(x = lon, y = lat) |>  #needs to have cols named 'x' and 'y'
  data.frame() |>  #can't be a tibble
  get_elev_raster(z = 10, prj = "epsg:4326", src = "aws", expand = 0.5,
                  clip = "bbox", ncpu = 5)
toc()  #took 37 sec

dem
dem <- rast(dem)  #convert to terra SpatRaster
plot(dem)
plot(africa$geometry, col = NA, border = "black", add = TRUE)
plot(gl_pa$geometry, col = NA, border = "black", lwd = 0.5, add = TRUE)


# Create plot w/ ggplot2
ggplot() +
  geom_spatraster(data = dem) +
  geom_sf(data = africa, color = "black", linewidth = 0.5, fill = NA) +
  geom_sf(data = gl_pa, color = "black", linewidth = 0.25, fill = NA) +
  # scale_fill_terrain_c() +
  scale_fill_hypso_c() +
  theme_bw() +
  coord_sf(xlim = ext(dem)[1:2],
           ylim = ext(dem)[3:4],
           expand = FALSE)


### Calculate slope, aspect, and terrain ruggedness index (TRI)

slope <- terrain(dem, v = "slope")
asp <- terrain(dem, v = "aspect")
tri <- terrain(dem, v = "TRI")  #a number of related metrics also available

plot(slope)
plot(asp)
plot(tri)




######################################
### Get vector water data from OSM ###
######################################

# Explore available features
available_features()

# See what's available for 'water'
available_tags(feature = 'water')

# Define region of interest
gl_bbox <- st_bbox(c(xmin = min(dat$lon), ymin = min(dat$lat), xmax = max(dat$lon), ymax = max(dat$lat))) + 
  c(-1, -1, 1, 1)

# Query water features w/in ROI
gl_water <- opq(bbox = gl_bbox) |>
  add_osm_feature(key = "water") |>
  osmdata_sf()

# Separate by class type
water_multipolygons <- gl_water$osm_multipolygons
water_polygons <- gl_water$osm_polygons

# Viz water
ggplot() +
  geom_sf(data = africa, color = "black", linewidth = 1) +
  geom_sf(data = water_multipolygons, color = 'grey30', fill = "lightblue2") +
  geom_sf(data = water_polygons, color = 'grey30', fill = "lightblue2") +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), linewidth = 0.15) +
  scale_color_brewer("ID", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = range(dat$lon),
           ylim = range(dat$lat))

# Compare to satellite map
dat |> 
  rename(x = lon, y = lat) |> 
  bayesmove::shiny_tracks(4326)





######################################################
### Stream in Human Footprint Index from the cloud ###
######################################################

# Values range from 0-50, where 0 is least impacted and 50 is most
# Pulling in latest global layer (from 2021)
# Data citation: F. Gassert, O. Venter, J. E.M. Watson, S.P. Brumby, J.C. Mazzariello, S.C. Atkinson, S. Hyde. 2023: "Global 100m Terrestrial Human Footprint (HFP-100) v1.2"
# Paper citation: Gassert F., Venter O., Watson J.E.M., Brumby S.P., Mazzariello J.C., Atkinson S.C. and Hyde S., An Operational Approach to Near Real Time Global High Resolution Mapping of the Terrestrial Human Footprint. Front. Remote Sens. 4:1130896 doi: 10.3389/frsen.2023.1130896 (2023). https://www.frontiersin.org/articles/10.3389/frsen.2023.1130896/full

# Stream in Human Footprint Index from cloud
hfp <- rast("https://data.source.coop/vizzuality/hfp-100/hfp_2021_100m_v1-2_cog.tif", vsi = TRUE)

# Define bounding box for data
bbox <- ext(min(dat$lon), max(dat$lon), min(dat$lat), max(dat$lat)) |>  #xmin, xmax, ymin, ymax
  vect(crs = "epsg:4326") |>  #define vector layer
  project(crs(hfp)) |>  #reproject to match Mollweide proj used by HFP
  ext() |>  #convert back to extent
  extend(10000)  #extend bbox in each direction by 10 km (i.e., 10,000 m)

# Crop HFP to match data spatial extent and project back to WGS84
hfp_gl <- crop(hfp, bbox)
hfp_gl <- project(hfp_gl, "epsg:4326")  #reproject back into WGS84
hfp_gl <- hfp_gl / 1000  #convert to original 0-50 scale

hfp_gl
plot(hfp_gl)


# Viz tracks and HFP
ggplot() +
  geom_spatraster(data = hfp_gl) +
  scale_fill_viridis_c("HFP") +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id))) +
  scale_color_brewer("ID", palette = "Set1") +
  theme_bw()





################################################################
### Stream in human population density (2020) from the cloud ###
################################################################

# Modeled residential population per 100 m grid cell for 2020, derived from the JRC Global Human Settlement Layer GHS-POP R2023A (but in WGS84 proj). 
# Values are census estimates from 2020, representing the number of ppl per grid cell
# Data taken from https://source.coop/cboettig/population w/ more details on original data at https://source.coop/nlebovits/ghsl


### Access data ###

# Stream-in 2020 dataset from cloud
pop <- rast("https://data.source.coop/cboettig/population/raw/ghs-pop-2020-cog.tif", vsi = TRUE)  #100 m res

# Define regional bbox on which to crop global raster
gl_bbox <- ext(c(xmin = min(dat$lon), xmax = max(dat$lon), ymin = min(dat$lat), ymax = max(dat$lat))) + 
  rep(1, 4)  #expand by 1° in every direction

# Crop raster
pop_gl <- crop(pop, gl_bbox)

pop_gl
plot(pop_gl)

ggplot() +
  geom_spatraster(data = pop_gl) +
  scale_fill_viridis_c("Human Pop.\nDensity") +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), linewidth = 0.15) +
  scale_color_brewer("ID", palette = "Set1") +
  theme_bw()


# Find pixels where pop > 15 ppl per 100 m
#essentially making subjective decision to identify higher-use areas; may want to adjust what this threshold is in a senstivity analysis

pop_gl_high <- pop_gl  #store original raster in new object
pop_gl_high[pop_gl_high < 15] <- NA  #change all small values (i.e., < 15 ppl) to NA

ggplot() +
  geom_sf(data = africa, fill = NA, linewidth = 0.5) +
  geom_spatraster(data = pop_gl_high) +
  scale_fill_viridis_c("Human Pop.\nDensity", na.value = "transparent") +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), linewidth = 0.15) +
  scale_color_brewer("ID", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = ext(pop_gl_high)[1:2],
           ylim = ext(pop_gl_high)[3:4])




### Calculate distance to human settlement ###

# Lower raster resolution (to speed up processing; from 100 to 500 m)
pop_high_500m <- aggregate(pop_gl_high, fact = 5, fun = "sum", na.rm = TRUE, cores = 5)

pop_high_500m
plot(pop_high_500m)


# Calculate distance surface to nearest non-NA pixel (i.e., distance to human settlement)
tic()
dist2pop <- distance(pop_high_500m)  #measured in meters (when using WGS84)
toc()  #5 sec; takes 30 min at 100 m res to calc distances

minmax(dist2pop)  #distance ranges from 0 to 104 km
plot(dist2pop)
points(as.points(pop_high_500m), col = "white", cex = 0.1)





### Creating polygon from population raster ###

## Binary method
bin_pop <- ifel(!is.na(pop_high_500m), 1, NA)  #need to convert to binary raster first
pop_high_500m_poly <- as.polygons(bin_pop)

# Can also convert to 'sf' POLYGON object
pop_high_500m_sf <- st_as_sf(pop_high_500m_poly) |> 
  st_cast("POLYGON")  #convert from MULTIPOLYGON so that each shape has its own record
plot(pop_high_500m_sf, border = NA)

ggplot() +
  geom_sf(data = za_mz) +
  geom_sf(data = pop_high_500m_sf, fill = "black", color = NA) +
  scale_fill_continuous("Human Pop. Density") +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), linewidth = 0.15, alpha = 0.5) +
  scale_color_brewer("ID", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = st_bbox(pop_high_500m_sf)[c("xmin","xmax")],
           ylim = st_bbox(pop_high_500m_sf)[c("ymin","ymax")])



## Aggregation method

# Group contiguous non-NA cells into unique patch IDs
pop_patches <- patches(pop_high_500m, directions = 8)  #using all 8 neighbors

# Aggregate raster values for each patch ID
patch_stats <- zonal(pop_high_500m, pop_patches, fun = "sum", na.rm = TRUE)

# Convert the patch raster into dissolved polygons
pop_polys <- as.polygons(pop_patches)

# Join the calculated stats back to the polygon attribute table
pop_polys <- merge(pop_polys, patch_stats, by = "patches") |> 
  st_as_sf()  #convert to 'sf'

ggplot() +
  geom_sf(data = za_mz) +
  geom_sf(data = pop_polys, aes(fill = `ghs-pop-2020-cog`), color = NA) +
  scale_fill_viridis_c("Human Pop. Density") +
  geom_path(data = dat, aes(lon, lat, group = id, color = factor(id)), linewidth = 0.15, alpha = 0.5) +
  scale_color_brewer("ID", palette = "Set1") +
  theme_bw() +
  coord_sf(xlim = st_bbox(pop_polys)[c("xmin","xmax")],
           ylim = st_bbox(pop_polys)[c("ymin","ymax")])






###############################################
### Access land use/land cover (LU/LC) data ###
###############################################

# Define study area as an sf object (EPSG:4326)
aoi <- st_as_sf(st_as_sfc(st_bbox(gl_bbox)), crs = 4326)

# Stream in ESA WorldCover with STAC parameters
tic()
worldcover_file <- get_stac_data(
  aoi = aoi,
  start_date = "2021-01-01",
  end_date = "2021-12-31",
  stac_source = "https://planetarycomputer.microsoft.com/api/stac/v1",
  collection = "esa-worldcover",
  asset_names = "map",
  sign_function = rsi::sign_planetary_computer,
  output_filename = tempfile(fileext = ".tif")
)
toc()  #took 4.5 min

# Load results and convert to classes
worldcover <- terra::rast(worldcover_file) |> 
  as.factor()

worldcover
names(worldcover) <- "lulc"  #change layer name
plot(worldcover)


# Define a Raster Attribute Table (RAT)
coltab(worldcover)  #print color table
wc_rat <- data.frame(
  ID = c(10, 20, 30, 40, 50, 60, 80, 90, 95),
  lulc = c("Tree cover", "Shrubland", "Grassland", "Cropland", 
           "Built-up", "Bare / sparse vegetation", 
           "Permanent water bodies", "Herbaceous wetland", "Mangroves")
)

# Set factor levels permanently on the SpatRaster pointer
levels(worldcover) <- wc_rat
levels(worldcover)


ggplot() +
  geom_spatraster(data = worldcover, use_coltab = TRUE) +
  scale_fill_coltab(name = "LULC", data = worldcover) +
  geom_sf(data = africa, fill = NA, linewidth = 0.5, color = "black") +
  geom_path(data = dat, aes(lon, lat, group = id), color = "black", linewidth = 0.15, alpha = 0.5) +
  theme_bw() +
  coord_sf(xlim = ext(worldcover)[1:2],
           ylim = ext(worldcover)[3:4],
           expand = FALSE)





##################################
### Access Sentinel 2A imagery ###
##################################

sm_aoi <- st_bbox(c(xmin = 32, ymin = -23, xmax = 32.25, ymax = -22.75), crs = 4326) |> 
  st_as_sfc() |> 
  st_as_sf() |> 
  st_transform(32736) # Needs to be in project CRS; reproject to meters in UTM

# Access Sentinel 2A imagery
tic()
s2a <- get_sentinel2_imagery(
  aoi = sm_aoi, 
  start_date = "2025-08-01",
  end_date = "2025-08-15",
  output_filename = tempfile(fileext = ".tif")
)
toc()  #took 2.5 min

# Viz image
s2a |> 
  rast(lyrs = c("R", "G", "B")) |>
  stretch() |>
  plotRGB()






##############################################
### Access MODIS NDVI data (16-day, 250 m) ###
##############################################

# Define bounding box for the STAC query (WGS84)
# bbox_wgs84 <- c(xmin = 32, ymin = -23, xmax = 32.25, ymax = -22.75)

# Define datetime range of tracking data (in format for STAC API)
date_range <- range(dat$date, na.rm = TRUE) |> 
  format("%Y-%m-%dT%H:%M:%SZ", tz = "UTC") |> 
  paste(collapse = "/")

# Verify the output
date_range
# "2022-03-17T06:53:00Z/2025-07-04T09:09:00Z"

# Query Microsoft Planetary Computer STAC for MODIS NDVI
tic()
items <- stac("https://planetarycomputer.microsoft.com/api/stac/v1") |>
  stac_search(
    collections = "modis-13Q1-061",
    bbox = st_bbox(gl_bbox),
    datetime = date_range
  ) |>
  get_request() |>
  items_fetch() |>  # Automatically loops through all API pages (otherwise limited to 250 max)
  items_sign(sign_planetary_computer()) # Appends the required SAS tokens
toc()  #took 4 sec

# Extract the specific NDVI remote URLs
ndvi_urls <- sapply(items$features, function(x) x$assets$`250m_16_days_NDVI`$href)  #304 layers

# Extract start and end datetimes directly from the STAC metadata
start_dates <- as_date(sapply(items$features, function(x) x$properties$start_datetime))
end_dates <- as_date(sapply(items$features, function(x) x$properties$end_datetime))

# Split URLs by date and check how many tiles per date
urls_by_date <- split(ndvi_urls, start_dates)  #148 dates
lengths(urls_by_date)  #nearly all have 2 tiles


# Build virtual raster dataset (VRT) to mosaic tiles per date
tic()
ndvi_vrt <- map(urls_by_date, ~{
  # Force streaming by prepending the VSI curl prefix to the HTTPS links
  vsi_urls <- paste0("/vsicurl/", .x)
  
  # Stitch the streaming URLs into a spatial mosaic
  vrt(vsi_urls)
  })
toc()  #took 1.5 min

ndvi_vrt$`2022-03-06`
plot(ndvi_vrt$`2022-03-06`)


# Create a directory to hold the cropped slices
temp_dir <- file.path(tempdir(), "modis_crops")
dir.create(temp_dir, showWarnings = FALSE)

# Iterate through each time step's VRT and crop it locally
tic()
cropped_paths <- map(seq_along(ndvi_vrt), function(i) {
  
  # Generate a unique file path for this time step
  out_file <- file.path(temp_dir, paste0("crop_", i, ".tif"))
  
  # Crop just this one layer and write directly to disk
  crop(
    ndvi_vrt[[i]], 
    gl_bbox |> 
      vect(crs = "epsg:4326") |>
      project(crs(ndvi_vrt$`2022-03-06`)), 
    filename = out_file, 
    overwrite = TRUE,
    wopt = list(gdal = c("COMPRESS=LZW"))
  )
  
  # Print progress (great for workshops!)
  message(sprintf("Cropped layer %d of %d", i, length(ndvi_vrt)))
  
  return(out_file)
})
toc()  #took 2.5 min


# Create single time series stack
ndvi <- rast(unlist(cropped_paths))

ndvi  #148 layers
plot(ndvi[[1]])


# Apply the scale factor
#Reported scale factor (0.0001) on NASA product website and in metadata is wrong
#Refer to forum thread: https://forum.earthdata.nasa.gov/viewtopic.php?t=5837&sid=b964068e3b5992875d02368abad8bfbc
ndvi2 <- ndvi * 1e-8  #apply correct scale factor
ndvi2

# Reproject raster to WGS84
tic()
ndvi3 <- project(ndvi2, "epsg:4326", threads = 10, use_gdal = TRUE, by_util = TRUE)
toc()  # took 20 sec



# Calculate the temporal midpoint
unique_starts <- unique(start_dates)
unique_ends <- unique(end_dates)
mid_dates <- unique_starts + (unique_ends - unique_starts) / 2

# Apply the calculated midpoints to your raster layer names
names(ndvi3) <- paste0("NDVI_", format(mid_dates, "%Y_%m_%d"))
time(ndvi3) <- mid_dates

### Explore time-varying NDVI layers ###
ndvi3
plot(ndvi3[[1:4]])

# Map latest 9 layers
ggplot() +
  geom_spatraster(data = ndvi3[[1:9]]) +
  scale_fill_viridis_c("NDVI", na.value = "transparent") +
  theme_bw() +
  facet_wrap(~ lyr)

# Compare difference between first and last layer of time series (2025-07-04 vs 2022-03-14)
ggplot() +
  geom_spatraster(data = ndvi3[[1]] - ndvi3[[nlyr(ndvi3)]]) +
  scale_fill_gradient2("Diff in NDVI", na.value = "transparent") +
  theme_bw()

# Calculate mean of all layers
ggplot() +
  geom_spatraster(data = mean(ndvi3, na.rm = TRUE)) +
  scale_fill_viridis_c("Mean NDVI", na.value = "transparent") +
  theme_bw()

# Calculate mean value per layer and plot over time
ndvi3 |> 
  global(fun = "mean", na.rm = TRUE) |> 
  data.frame() |> 
  mutate(date = time(ndvi3)) |> 
  
  ggplot() +
  geom_line(aes(date, mean)) +
  geom_point(aes(date, mean)) +
  theme_bw()


### Calculate monthly mean value (based on midpoint date)

# Create grouping index
ymon_idx <- format(mid_dates, "%Y_%m")

# Apply the mean function across the groups using tapp()
monthly_ndvi <- tapp(
  ndvi3, 
  index = ymon_idx, 
  fun = mean, 
  na.rm = TRUE  #ensures any masked pixels don't ruin calc
)

# Create 'clean' names
names(monthly_ndvi) <- paste0("NDVI_Avg_", str_remove(names(monthly_ndvi), "X"))

monthly_ndvi
plot(monthly_ndvi)

ggplot() +
  geom_spatraster(data = monthly_ndvi[[1:9]]) +
  scale_fill_viridis_c("NDVI", na.value = "transparent") +
  theme_bw() +
  facet_wrap(~ lyr)





############################
### Export raster layers ###
############################
#for extraction and model prediction

# DEM
writeRaster(dem, "rasters/dem.tif", overwrite = TRUE)

# Distance to (high) human pops
writeRaster(dist2pop, "rasters/dist2pop.tif", overwrite = TRUE)

# Polygons of high human pop areas
st_write(pop_polys, "rasters/pop_polygons.fgb", overwrite = TRUE)

# LULC
writeRaster(worldcover, "rasters/lulc.tif", overwrite = TRUE)  #350 MB

# NDVI
ndvi_sorted <- sort(ndvi3)  #first sort raster in chronological order
writeCDF(ndvi_sorted, "rasters/ndvi.nc", varname = "NDVI", timename = "time", overwrite = TRUE)  #2.6 GB
