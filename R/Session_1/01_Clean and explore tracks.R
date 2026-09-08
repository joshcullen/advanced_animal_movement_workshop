
### Clean and explore data ###

library(tidyverse)
library(readxl)
library(plotly)
library(bayesmove)
# install.packages("rnaturalearthhires", repos = "https://ropensci.r-universe.dev", type = "source")
library(rnaturalearth)
library(sf)
library(scico)


###################
#### Load data ####
###################

dat <- read_excel("raw_data/GL_bulls_example.xlsx")

# Explore data summaries
glimpse(dat)
summary(dat)
n_distinct(dat$`Elephant Id`)

## Structure of dataset ##
#Elephant Id: numeric animal ID (assigned to tag)
#Time Stamp: datetime of each relocation
#Date: the date of each relocation (in datetime format via read_excel())
#Time: the time of each relocation (in datetime format via read_excel())
#Latitude: latitude of relocation
#Longitude: longitude of relocation
#Temperature: the temperature reading from the on-board sensor
#Accelerometer: values associated w/ the x/y/z dims of the tri-axial accelerometer


# Inspect columns of pertinent variables
table(dat$`Elephant Id`, useNA = "ifany")

#Now, we'll need to make 'better' column names and convert temperature to an integer



##################################################################################
#### Rename columns, convert temperature to 'integer' class, and sort by date ####
##################################################################################

dat2 <- dat |> 
  select(-c(Date, Time, Accelerometer)) |>  #remove these columns since not needed
  rename(id = `Elephant Id`,
         date = `Time Stamp`,
         lon = Longitude,
         lat = Latitude,
         temp = Temperature) |> 
  mutate(temp = str_match(temp, pattern = "[0-9]+") |>  #only keep numbers
           as.integer()) |> 
  arrange(id, date)  #make sure data is properly sorted



#########################################
#### Check for and remove duplicates ####
#########################################

# Check for any duplicate records
dat2 |> 
  group_by(id, date) |>  #define groups as a timestamp per ID
  filter(n() > 1) |>  #only keep records where >1 occur (i.e., duplicates)
  ungroup()
#4 pairs total; 2 for ID 5605, 1 for ID 6469, 1 for ID 6471

# Remove duplicate records
dat3 <- dat2 |>
  group_by(id, date) |>
  distinct(date, .keep_all = TRUE) |>  #only keep one copy for each timestamp per ID
  ungroup()



#############################################
#### Check basic summary stats on tracks ####
#############################################

dat3 |> 
  summarize(.by = id,
            start = first(date),
            end = last(date),
            n = n(),
            median_dt = round(median(difftime(dat2$date, lag(dat2$date), units = "hours"),
                               na.rm = TRUE), 2),
            min_temp = min(temp, na.rm = TRUE),
            max_temp = max(temp, na.rm = TRUE)) |> 
  mutate(days = round(end - start, 1),
         .after = end)



########################################################
#### Visualize tracks to determine if any anomalies ####
########################################################

# Change 'id' to character so treated as discrete values (not continuous)
dat3$id <- as.character(dat3$id)

ggplot(dat3, aes(lon, lat, color = id)) +
  geom_path(aes(group = id), linewidth = 0.25) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  coord_equal()
#no visible outliers


# Let's add some land layers to give better context
africa <- ne_countries(scale = 10, continent = c("Africa"), returnclass = "sf")

ggplot() +
  geom_sf(data = africa) +
  geom_path(data = dat3, aes(lon, lat, group = id, color = id), linewidth = 0.25) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  theme(panel.grid = element_blank()) +
  coord_sf(xlim = c(min(dat3$lon) - 1, max(dat3$lon) + 1),
           ylim = c(min(dat3$lat) - 1, max(dat3$lat) + 1))


# Let's make this map interactive so we can see a little better
plotly::ggplotly(
  ggplot() +
    geom_sf(data = africa) +
    geom_path(data = dat3, aes(lon, lat, group = id, color = id,
                               text = paste(  #create custom tooltip on hover
                                 "ID:", id,
                                 "\nLon:", lon,
                                 "\nLat", lat,
                                 "\nDate:", date)),
              linewidth = 0.25) +
    scale_color_brewer(palette = "Set1") +
    theme_bw() +
    theme(panel.grid = element_blank()) +
    coord_sf(xlim = c(min(dat3$lon) - 1, max(dat3$lon) + 1),
             ylim = c(min(dat3$lat) - 1, max(dat3$lat) + 1)),
  tooltip = "text"
)
#there appear to be some limitations to movement (e.g., fences)


# Let's use a Shiny app from {bayesmove} to explore a little further
dat3 |>
  rename(x = lon, y = lat) |>
  bayesmove::shiny_tracks(epsg = 4326)



################################################
#### Explore how locations change over time ####
################################################

# Inspect time series plots of lat and long
dat3 |> 
  pivot_longer(cols = c(lon, lat), names_to = "coords", values_to = "values") |> 
  
  ggplot() +
  geom_line(aes(date, values, color = id)) +
  scale_color_brewer(palette = "Set1") +
  theme_bw() +
  facet_grid(coords ~ id, scales = "free")





############################################
#### Explore temperature data over time ####
############################################

# By date
dat3 |> 
  
  ggplot() +
  geom_point(aes(date, temp, color = temp), alpha = 0.05) +
  scale_color_viridis_c(option = "inferno") +
  labs(x = "Date", y = "Temperature (°C)") +
  theme_bw(base_size = 14) +
  facet_wrap(~id, ncol = 1)


# By hour
dat3 |> 
  
  ggplot() +
  geom_point(aes(hour(date), temp, color = month(date)), alpha = 0.05) +
  scico::scale_color_scico("Month", palette = "romaO") +  #circular palette
  labs(x = "Hour", y = "Temperature (°C)") +
  theme_bw(base_size = 14) +
  facet_wrap(~id, ncol = 1)




#############################
#### Export cleaned data ####
#############################

write_csv(dat3, "processed_data/Session_1/cleaned_tracks.csv")
