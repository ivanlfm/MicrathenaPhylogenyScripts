setwd("c:/a/maps")

# Load necessary libraries
library(ggplot2)
library(sf)
library(terra)

# Load data
data <- read.csv("Micrathena_all_distributions.csv")

# Create a list of unique species
species_list <- unique(data$Species)

# Convert the data frame to sf object
species_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)

# Load country borders shapefile
country_borders <- st_read("world-administrative-boundaries.shp")
country_borders_cropped <- st_crop(country_borders, c(xmin= min(data$Longitude)-10, ymin = min(data$Latitude)-10, xmax = max(data$Longitude)+10, ymax = max(data$Latitude)+10))

# Loop through each species
for (species in species_list) {
  # Subset data for current species
  species_data <- subset(species_sf, Species == species & Type != "inaturalist")
  dubious <- subset(species_sf, Species == species & Type == "inaturalist")
  # Plot the distribution of the species on the background map
  ggplot() +
    geom_sf(data = country_borders_cropped, aes(fill=NULL), color = "black", fill = "light green")  +
	geom_sf(data = dubious, color = "black", size = 2.7) +
	geom_sf(data = dubious, color = "red", size = 2.5) +
	geom_sf(data = species_data, color = "white", size = 2.2) +
	geom_sf(data = species_data, color = "black", size = 2) +
	ggtitle(species)
  # Save the plot as a PNG image
  ggsave(paste(species, ".png"))
}