
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