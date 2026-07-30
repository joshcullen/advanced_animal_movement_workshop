
### Clean and explore data using {amt} ###

library(tidyverse)
library(amt)
library(sf)


#################
### Load data ###
#################

# Load built-in dataset of fishers
data("amt_fisher")


# Convert tibble to amt 'track' object
trk <- amt_fisher |> 
  # rbind(amt_fisher[1,]) |> 
  make_track(x_, y_, t_, crs = 5070, all_cols = TRUE, check_duplicates = TRUE, verbose = TRUE)

# Check data completeness and format
trk
summary(trk)  #coordinates reported in NAD83 proj (based on help docs)
table(trk$name)  #4 tracks




####################################
### Perform exploratory analyses ###
####################################

# Viz tracks
ggplot(trk, aes(x_, y_, color = name, group = name)) +
  geom_path() +
  theme_bw() +
  coord_equal()

# Nest by ID (for subsequent ID-specific analyses)
trk2 <- trk |>
  nest(data = -"name")  #converts to compact object, w/ 1 row per ID

trk2  #we now have a tibble where all data besides ID is now a list
trk2$data[[1]]  #print first 10 obs for Leroy



### Calculate NSD per ID and explore time series

# Calc NSD per ID
trk3 <- trk2 |> 
  mutate(nsd = map(data, ~add_nsd(.x)))

trk3
trk3$nsd[[1]]  #We can see NSD added to rest of data now

# Viz NSD over time
trk3 |> 
  select(name, nsd) |> 
  unnest(cols = nsd) |>
  
  ggplot() +
  geom_path(aes(t_, nsd_, color = name, group = name)) +
  labs(x = "Time", y = "Net Squared Displacement (m^2)") +
  theme_bw(base_size = 14) +
  facet_wrap(~ name, scales = "free_x")




### Calculate metrics of movement patterns

trk3 |> 
  mutate(straightness = map(data, straightness) |> unlist(),
         tot_dist = map(data, tot_dist) |> unlist(),
         sinuosity = map(data, sinuosity) |> unlist())
# Leroy moved the straightest and farthest (> 2 km)




# Remove first day of data that could produce a 'capture effect'
trk4 <- trk2 |> 
  mutate(data2 = map(data, ~remove_capture_effect(.x, start = period(1, "day"))))

# Check for Leroy
trk4$data[[1]]  #starts on 2009-02-11
trk4$data2[[1]]  #starts on 2009-02-13


# Traditional removal of first day of data
trk |> 
  mutate(date = as_date(t_)) |> 
  group_by(name) |> 
  filter(date > first(date)) |> 
  ungroup()

#Or exactly 24 hrs after first transmission
trk |> 
  group_by(name) |> 
  filter(t_ > first(t_) + days(1)) |> 
  ungroup()



### Summarize sampling intervals

# Summarize across IDs
summarize_sampling_rate_many(amt_fisher, "name")
#looks like the schedule changes occasionally

# Viz time intervals over tracking duration
trk5 <- trk4 |> 
  mutate(data3 = map(data2, ~{
    .x |> 
      mutate(dt = c(summarize_sampling_rate(.x, summarize = FALSE), NA))
  })
  ) |> 
  select(name, data3) |> 
  unnest(cols = data3)

trk5  #dt units are shown in minutes


ggplot(trk5) +
  geom_point(aes(t_, dt, color = name)) +
  labs(x = "Time", y = "Time Interval (min)") +
  theme_bw(base_size = 14) +
  facet_wrap(~ name, scales = "free_x")
#occasional outliers present per ID





############################################
### Example for converting to sf objects ###
############################################

# Convert to sf objects
trk_sf_pt <- as_sf_points(trk)  #convert to sf POINT object
trk_sf_pt

trk_sf_line <- trk4 |> 
  mutate(sf_lines = map(data, as_sf_lines)) |>  #convert to sf LINESTRING object
  pull(sf_lines) |>  #pull the 'sf_lines' column only
  set_names(trk4$name) |>  #provide fisher names
  bind_rows(.id = "name") |>  #bind into single data.frame
  rename(geometry = 2)  #give proper name for spatial column
trk_sf_line

# Generate bounding box
bbox <- bbox(trk)
bbox

# Viz spatial objects
ggplot() +
  geom_sf(data = trk_sf_line, linewidth = 0.5) +
  geom_sf(data = trk_sf_pt, aes(color = name), size = 0.5) +
  geom_sf(data = bbox, fill = NA, color = "firebrick", linewidth = 1) +
  theme_bw(base_size = 14)



### Calculate track-level spatial objects

# Calculate bbox and centroid per ID
trk6 <- trk4 |> 
  mutate(bbox = map(data, amt::bbox),
         centroid = map(data, amt::centroid))

# Merge each into single object
id_bbox <- trk6 |> 
  pull(bbox) |> 
  set_names(trk6$name) |> 
  map(sf::st_sf) |> 
  bind_rows(.id = 'name') |> 
  rename(geometry = 2)

cent <- trk6 |> 
  select(name, centroid) |> 
  mutate(centroid = map(centroid, ~{.x |> t() |> as.data.frame()})) |> 
  unnest(cols = centroid)

ggplot() +
  geom_sf(data = trk_sf_line, linewidth = 0.1, alpha = 0.3) +
  geom_sf(data = id_bbox, aes(color = name), linewidth = 0.5, fill = NA) +
  geom_point(data = cent, aes(x_, y_, color = name), size = 3) +
  theme_bw(base_size = 14)
#it may be simpler to do this w/ basic functions where you derive these yourself

