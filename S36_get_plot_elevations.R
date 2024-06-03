#assigns elevation to points based on lat/long and and elevation raster
#generates a graph with altitudes for each species
#generates individual maps for each species, highlighting high-elevation records
#generates phylip matrices
#thanks chatGPT!

# Load required packages
library(sp)
library(raster)
library(ggplot2)
library(dplyr)
library(sf)
library(terra)

#elevation threshold to consider a species a lowland/mountain taxa
elev_threshold <- 1200 #Following Pérez-Escobar et al. 2022, p.10, https://doi.org/10.1016/j.tplants.2021.09.010

# Read coordinate data in raster file with altitude, and shapefile with world borders
setwd("c:/a/maps")
coords <- read.csv("Micrathena_all_distributions.csv")
raster_file <- raster("raster/wc2.1_30s_elev.tif")
country_borders <- st_read("world-administrative-boundaries.shp")

##################
###GET ELEVATION
##################

#gets coordinates
coordinates(coords) <- c("Longitude", "Latitude")
proj4string(coords) <- CRS("+proj=longlat +datum=WGS84")

# Extract altitude values for each point
altitudes <- extract(raster_file, coords)

#corrects negative and NA altitudes to 0
altitudes[altitudes < 0] <- 0
altitudes[is.na(altitudes)] <- 0

# Combine coordinates and altitude into new dataframe
data <- data.frame(coords, altitude = altitudes)
data$optional <- NULL

#creates column indicating if record is high-elevation
for(i in 1:nrow(data)){
	if(data$altitude[i] < elev_threshold) {
		data[i, 7] <- "low_elevation"
	} else {
		data[i, 7] <- "high_elevation"
	}
}

colnames(data)[colnames(data) == 'V7'] <- "type_elevation"

# Write out new dataframe to CSV file with elevation
write.csv(data, "output_elevation.csv", row.names = FALSE)

##################
###PLOT ALTITUDES
##################

# Calculate mean altitude for each species
median_altitudes <- data %>% 
  group_by(Species) %>% 
  summarize(median_altitude = median(altitude)) %>% 
  arrange(median_altitude)

# Create plot with all species ordered by mean altitude; marks median and 25/75 quantiles
plot <- ggplot(data, aes(x = factor(Species, levels = median_altitudes$Species), y = altitude)) +
  geom_point(aes(color = Species), alpha = 0.1, size = 3) +
  geom_hline(yintercept = elev_threshold, linetype = "dashed", color = "gray") +
  #stat_summary(fun = function(altitude) {quantile(altitude,0.25)}, geom = "point", alpha = 0.5, size = 1, color = "black", shape = 23, fill = "red") +
  #stat_summary(fun = function(altitude) {quantile(altitude,0.75)}, geom = "point", alpha = 0.5, size = 1, color = "black", shape = 23, fill = "red") +
  #stat_summary(fun = mean, geom = "point", size = 3, color = "black", shape = 21, fill = "red") +
  stat_summary(fun = median, geom = "point", size = 3, color = "black", shape = 21, fill = "gold") +
  labs(x = "Species", y = "Altitude", title = "Altitude by Species") +
  theme_minimal() +
  scale_color_manual(values = rep("black", length(unique(data$Species)))) +
  coord_flip() +
  guides(color = "none")

# Save plot as PDF and JPEG
ggsave("all_species_plot.pdf", plot, width = 8.27, height = 11.69, units = "in", device = "pdf")
ggsave("all_species_plot.jpeg", plot, width = 8.27, height = 11.69, units = "in", device = "jpeg")

##########################################
###GENERATE PHYLIP FILES FOR BIOGEOBEARS
##########################################

# Calculate median altitude for each species
median_altitudes <- data %>% 
  group_by(Species) %>% 
  summarize(median_altitudes = median(altitude))
  median_altitudes <- as.data.frame(median_altitudes)

phylip_matrix <- matrix(data=0, nrow=length(median_altitudes$Species), ncol = 3)
phylip_matrix <- as.data.frame(phylip_matrix)
names(phylip_matrix) <- c("species", "lowland", "mountain")

for(i in 1:length(median_altitudes$Species)){
	phylip_matrix [i,1] <- median_altitudes$Species[i]
	if(median_altitudes$median_altitudes[i] < elev_threshold) phylip_matrix [i,2] <- 1
	if(median_altitudes$median_altitudes[i] >= elev_threshold) phylip_matrix [i,3] <- 1
}

write.csv(phylip_matrix, "phylip_threshold.csv", row.names = FALSE)

# Calculate quantiles of altitude for each species
quantiles_altitudes_low <- data %>% 
  group_by(Species) %>% 
  summarize(quants_low = quantile(altitude, probs = 0.25))
  quantiles_altitudes_low <- as.data.frame(quantiles_altitudes_low)
  
  quantiles_altitudes_high <- data %>% 
  group_by(Species) %>% 
  summarize(quants_high = quantile(altitude, probs = 0.75))
  quantiles_altitudes_high <- as.data.frame(quantiles_altitudes_high)
  
  quantiles_altitudes <- cbind(quantiles_altitudes_low, quantiles_altitudes_high$quants_high)
  names(quantiles_altitudes) <- c("species", "quant25", "quant75")
  
phylip_matrix_quantiles <- matrix(data=0, nrow=length(quantiles_altitudes$species), ncol = 3)
phylip_matrix_quantiles <- as.data.frame(phylip_matrix_quantiles)
names(phylip_matrix_quantiles) <- c("species", "lowland", "mountain")

for(i in 1:length(quantiles_altitudes$species)){
	phylip_matrix_quantiles [i,1] <- quantiles_altitudes$species[i]
	if(quantiles_altitudes$quant25[i] < elev_threshold) phylip_matrix_quantiles [i,2] <- 1
	if(quantiles_altitudes$quant75[i] >= elev_threshold) phylip_matrix_quantiles [i,3] <- 1
}

write.csv(phylip_matrix_quantiles, "phylip_quantiles.csv", row.names = FALSE)
  
###############
###PLOTS MAPS
###############

# Create a list of unique species
species_list <- unique(data$Species)

# Convert the data frame to sf object
species_sf <- st_as_sf(data, coords = c("Longitude", "Latitude"), crs = 4326)

#crops world shapefile
country_borders_cropped <- st_crop(country_borders, c(xmin= min(data$Longitude)-10, ymin = min(data$Latitude)-10, xmax = max(data$Longitude)+10, ymax = max(data$Latitude)+10))

if (!file.exists("species_maps")) dir.create("species_maps/")
if (!file.exists("species_maps_color_ramp")) dir.create("species_maps_color_ramp/")

# Loop through each species
for (species in species_list) {
  # Subset data for current species
  species_data_low <- subset(species_sf, Species == species & type_elevation == "low_elevation")
  species_data_high <- subset(species_sf, Species == species & type_elevation == "high_elevation")
  # Plot the distribution of the species on the background map
  ggplot() +
    geom_sf(data = country_borders_cropped, aes(fill=NULL), color = "black", fill = "light green")  +
	geom_sf(data = species_data_low, shape = 21, color = "white", size = 2, fill = "black") +
	geom_sf(data = species_data_high, shape = 21, color = "white", size = 2, fill = "red") +
	ggtitle(species)
  # Save the plot as a PNG image
  ggsave(paste("species_maps/", species, ".png"), width = 4, height = 4, units = "in", device = "png")
  }
  
alt_range <- range(data$altitude)
 for (species in species_list) {
  # Subset data for current species
  species_data <- subset(species_sf, Species == species)
  # Plot the distribution of the species on the background map
  ggplot() +
    geom_sf(data = country_borders_cropped, aes(fill=NULL), color = "black", fill = "light green")  +
	geom_sf(data = species_data, aes(fill = altitude), shape = 21, color = "white", size = 2) +
	scale_fill_gradient(name = "Elevation", low = "black", high = "red", limits = alt_range,breaks = pretty(alt_range, n = 5)) +
	 theme(legend.position = "right", legend.key.width = unit(0.2, "cm"),
        legend.title = element_text(size = 6, face = "bold"), 
        legend.text = element_text(size = 5)) +
	ggtitle(species)
  # Save the plot as a PNG image
  ggsave(paste("species_maps_color_ramp/", species, ".png"), width = 4, height = 4, units = "in", device = "png")
  }