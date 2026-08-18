
# Function to add transformed coordinates

add_trans_coords <- function(data, coords, proj, new_proj) {
  ## data = a data.frame with coordinates to be re-projected
  ## coords = char vector; The names of the coordinates to be tranformed
  ## proj = either an EPSG code or proj string defining the projection of `coords`
  ## new_proj = an EPSG code or proj string defining the projection of the new coords to be added
  
  tmp <- sf::st_as_sf(data, coords = coords, crs = proj, remove = FALSE)
  
  if (sum(coords %in% c('x','X','y','Y')) > 0) {
    tmp |> 
      sf::st_transform(4326) |> 
      mutate(lon = unlist(map(geometry, 1)),
             lat = unlist(map(geometry, 2))) |>  #add cols from new CRS
      sf::st_drop_geometry()
  } else if (sum(coords %in% c('lon','long','longitude','Lon','Long','Longitude',
                               'lat','latitude','Lat','Latitude')) > 0) {
    tmp |> 
      sf::st_transform(new_proj) |> 
      mutate(x = unlist(map(geometry, 1)),
             y = unlist(map(geometry, 2))) |>  #add cols from new CRS
      sf::st_drop_geometry()
  } else {
    stop("Coordinate names must either be x/y or lon/lat.")
  }
  
}

#-----------------------------------

# Function to create common SpatRaster grid across IDs

create_shared_grid <- function(data, rasters, crs, res, ext = 0.3) {
  ## data = data.frame with coordinates labeled 'x' and 'y'
  ## rasters = a list object containing the SpatRasters of interest; list or raster layers should be named by ID for clarity
  ## crs = projection to be used for raster; could be EPSG code or proj4 string; For example, UTM Zone 36S could be specified as "EPSG:32736"
  ## res = numeric. Value of spatial resolution in units of specified CRS
  ## ext = the proportion by which to extend the bounding box around the data
  
  dat_bbox <- terra::ext(min(data$x), max(data$x), min(data$y), max(data$y))
  
  # Extend bbox
  dat_bbox2 <- dat_bbox * (1 + ext)
  
  # Create shared template grid
  r_template <- terra::rast(dat_bbox2, resolution = res, crs = crs)
  
  # Resample to template
  r <- lapply(rasters, \(x) terra::resample(x, r_template, method = "mean", threads = TRUE))
  
  # Replace NAs with 0 (areas where the animal didn't go should have 0 probability)
  r2 <- lapply(r, \(x) terra::subst(x, from = NA, to = 0)) |> 
    terra::rast()  #create stack
  
  # Add names
  if(!is.null(names(rasters))) {
    names(r2) <- names(rasters)
  } else if (lapply(rasters, names) |> unlist() |> unique() |> length() > 1) {
    names(r2) <- lapply(rasters, names) |> unlist()
  }
  
  return(r2)
}

#-----------------------------------

# Function to calculate space use overlap indices from SpatRasters
#works on whatever raster layers provided (e.g., full surface, 95% UD, 50% UD, etc)

calc_ud_overlap <- function(r, index, feature = NULL) {
  ## r = a SpatRaster object where the layers represent different UDs on shared grid
  ## index = character. "vi" (Volume of Intersection), "ba" (Bhattacharyya's Affinity), or "feature" (Overlap with spatial feature).
  ## feature = an sf or SpatVector polygon (required if index == "feature").
  
  
  # Normalize the SpatRaster so each layer sums to 1
  layer_sums <- terra::global(r, "sum", na.rm = TRUE)$sum
  r_norm <- r / layer_sums
  
  
  # ---------------------------------------------------------
  # Feature Overlap (Probability Mass inside a Polygon)
  # ---------------------------------------------------------
  if (index == 'feature') {
    if (is.null(feature)) {
      stop("You must provide an 'sf' or 'SpatVector' object to the 'feature' argument.")
    }
    
    # Convert sf to terra's SpatVector format
    if (inherits(feature, "sf")) {
      feature <- terra::vect(feature)
    }
    
    # Ensure CRS matches between raster and feature
    if (terra::crs(feature) != terra::crs(r_norm)) {
      feature <- terra::project(feature, terra::crs(r_norm))
    }
    
    # Mask the raster to the feature (sets all cells outside the polygon to NA)
    r_masked <- terra::mask(r_norm, feature)
    
    # Sum the remaining density values inside the polygon
    overlap_sums <- terra::global(r_masked, "sum", na.rm = TRUE)$sum
    
    overlap <- tibble(
      id = names(r_norm),
      overlap = overlap_sums
    )
    
    # ---------------------------------------------------------
    # Pairwise Overlap (VI or BA)
    # ---------------------------------------------------------
  } else if (index %in% c('vi', 'ba')) {
    
    # Extract spatial data into a list of numeric vectors for fast pairwise math
    vals_list <- as.list(r_norm) |> 
      map(~ as.numeric(terra::values(.x))) |> 
      set_names(names(r_norm))
    
    overlap <- combn(names(vals_list), 2, simplify = FALSE) |> 
      map_dfr(function(pair) {
        id1 <- pair[1]
        id2 <- pair[2]
        
        if (index == 'vi') {
          # Calculate Volume of Intersection
          val <- sum(pmin(vals_list[[id1]], vals_list[[id2]], na.rm = TRUE))
        } else if (index == 'ba') {
          # Calculate Bhattacharyya's Affinity
          val <- sum(sqrt(vals_list[[id1]] * vals_list[[id2]]), na.rm = TRUE)
        }
        
        tibble(
          from = id1,
          to = id2,
          overlap = val
        )
      })
    
  } else {
    stop("Index must be 'vi', 'ba', or 'feature'.")
  }
  
  return(overlap)
}
