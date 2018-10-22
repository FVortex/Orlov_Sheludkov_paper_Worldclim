# Maps for the paper BIOCLIMATIC DATA OPTIMIZATION FOR SPATIAL DISTRIBUTION MODELS
# Author: Alexander Sheludkov
# Date: 7 August 2018

library(sp)
library(raster)
library(rgdal)
library(rgeos)
library(dplyr)
library(ggplot2)
library(sf)
library(broom)

# ===============
# 1. Reading data
# ===============

# Bioclimatic clusters
load("data/initial_full_bioclimatic_data.Rda")
load("data/reduced_full_bioclimatic_data.Rda")
initial_bioclim <- initial_full_bioclimatic_data %>% 
  as.data.frame() %>%
  select(Longitude, Latitude, Cluster)
reduced_bioclim <- reduced_full_bioclimatic_data %>% 
  as.data.frame() %>%
  select(Longitude, Latitude, Cluster)

# Actual localities
load("data/actual_localities.RData")

# Predictions
load("data/initial_dataset_prediction.rda")
load("data/reduced_dataset_prediction.rda")


# =================
# 2. Pre-processing
# =================

# ========
# 2.1. CRS
# ========

# Define CRS
WGS84 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")

# ==========================
# 2.2. Bioclim data clusters
# ==========================

# Create SpatialPointsDataFrame objects
initial_clusters <-
  SpatialPointsDataFrame(c(initial_bioclim[,1:2]), data = initial_bioclim, 
                         proj4string = WGS84)
reduced_clusters <- 
  SpatialPointsDataFrame(c(reduced_bioclim[,1:2]), data = reduced_bioclim, 
                         proj4string = WGS84)

# =======================
# 2.1. Actual localitites
# =======================

# Create data frame
colnames(actual_localities) <- c("Latitude","Longitude")
# Change the order of the columns
actual_localities %>% select(Longitude, Latitude) -> actual_localities

# Create SpatialPoints object
actual_points <- SpatialPoints(actual_localities, proj4string = WGS84)

# ========================
# 2.2. Initial predictions
# ========================

# Firstly, transform the data into tidy data frames

# Scrap models' names 
model_names <- initial %>% names()

# Join the data into single data frame
init_preds <- data_frame()
for(i in 1:length(model_names)){
  initial[[i]] %>% 
    mutate(model = model_names[i]) %>%
    select(Longitude, Latitude, model) %>% 
    bind_rows(init_preds) -> init_preds
}

# Change "model" to factor variable
init_preds %>% 
  mutate(model = factor(model, ordered = T, levels = model_names)) ->
  init_preds

# ========================
# 2.3. Reduced predictions
# ========================

# Join the data into single data frame
reduced_preds <- data_frame()
for(i in 1:length(model_names)){
  reduced[[i]] %>% 
    mutate(model = model_names[i]) %>%
    select(Longitude, Latitude, model) %>% 
    bind_rows(reduced_preds) -> reduced_preds
}

# Change "model" to factor variable
reduced_preds %>% 
  mutate(model = factor(model, ordered = T, levels = model_names)) ->
  reduced_preds

# ================
# 2.4. Aggregation
# ================

# Aggregate model results by raster cells with 10 min resolution

# Сreate empty raster
cell_size <- 1/6 # 10 min
lon_min <- 32; lon_max <- 36.9; lat_min <- 44.2; lat_max <- 46.3
ncols <- ((lon_max - lon_min)/cell_size)+1; nrows <- ((lat_max - lat_min)/cell_size)+1
temp_raster <- raster(nrows=nrows, ncols=ncols, xmn=lon_min, xmx=lon_max, 
                      ymn=lat_min, ymx=lat_max, crs="+proj=longlat +datum=WGS84")

# 2.4.1. Initial predictions

# Create empty data_frame
init_preds_aggr <- data_frame()

# Fill data frame with aggregated data
for(i in 1:length(levels(init_preds$model))){
  print(i)
  
  # Filter df and select coordinates columns
  init_preds %>% 
    filter(model == levels(init_preds$model)[i]) %>% 
    select(Longitude, Latitude) -> temp_df
  # Rasterize 
  temp_raster <- rasterize(temp_df, temp_raster, fun = "count")
  temp_spdf <- as(temp_raster, "SpatialPixelsDataFrame")
  temp_spdf <- as.data.frame(temp_spdf)
  colnames(temp_spdf) <- c("value", "x", "y")
  # Add model variable
  temp_spdf %>% 
    mutate(model = levels(init_preds$model)[i]) -> temp_spdf
  # Add results to final data frame
  temp_spdf %>% bind_rows(init_preds_aggr) -> init_preds_aggr
}

# Change "model" to factor variable
init_preds_aggr %>% 
  mutate(model = factor(model, ordered = T, levels = levels(init_preds$model))) ->
  init_preds_aggr

# 2.4.2. Reduced predictions

# Create empty data_frame
reduced_preds_aggr <- data_frame()

for(i in 1:length(levels(reduced_preds$model))){
  print(i)
  # Filter df and select coordinates columns
  reduced_preds %>% 
    filter(model == levels(reduced_preds$model)[i]) %>% 
    select(Longitude, Latitude) -> temp_df
  # Rasterize 
  temp_raster <- rasterize(temp_df, temp_raster, fun = "count")
  temp_spdf <- as(temp_raster, "SpatialPixelsDataFrame")
  temp_spdf <- as.data.frame(temp_spdf)
  colnames(temp_spdf) <- c("value", "x", "y")
  # Add model variable
  temp_spdf %>% 
    mutate(model = levels(reduced_preds$model)[i]) -> temp_spdf
  # Add results to final data frame
  temp_spdf %>% bind_rows(reduced_preds_aggr) -> reduced_preds_aggr
}

# Change "model" to factor variable
reduced_preds_aggr %>% 
  mutate(model = factor(model, ordered = T, levels = levels(reduced_preds$model))) ->
  reduced_preds_aggr

# ===========
# 3. Basemape
# ===========

# Loading costline data from OSM
coastline <- readOGR("data/coastline_Crimea_OSM.geojson") %>% 
  aggregate()
# We can take a look
plot(coastline)

# Create bounding box and clip coastline
bbox <- extent(c(32, 36.9, 44.2, 46.3)) %>% as("SpatialPolygons")
crs(bbox) <- crs(coastline)
# Clip coastline by bbox
raster::intersect(coastline, bbox) -> coastline

# Take a look at the basemap
ggplot()+
  geom_sf(data = st_as_sf(coastline))+
  theme_minimal()

# ======================
# 4. Final Visualization
# ======================

# 4.1. Bioclimatic clusters

# Initial Grey
p_initial_bioclim <- ggplot()+
  geom_sf(data = st_as_sf(coastline))+
  geom_raster(data = tidy(initial_clusters), 
              aes(x = Longitude, y= Latitude, fill = as.factor(Cluster)), alpha = 0.8)+
  scale_fill_grey(name = "Bioclim\ncluster")+
  scale_y_continuous(name = element_blank())+
  scale_x_continuous(name = element_blank())+
  # ggtitle(label = "Initial")+
  theme_minimal(base_size = 14)

# Initial Colored
p_initial_bioclim <- ggplot()+
  geom_sf(data = st_as_sf(coastline))+
  geom_raster(data = data_frame(Longitude = coordinates(initial_clusters)[,1], 
              Latitude = coordinates(initial_clusters)[,2], Cluster = initial_clusters@data$Cluster), 
              aes(x = Longitude, y= Latitude, fill = as.factor(Cluster)), alpha = 0.8)+
  scale_fill_viridis_d(name = "Bioclim\ncluster")+
  # scale_fill_viridis_d(name = "Bioclim\ncluster", option = "B", direction = -1)+
  scale_y_continuous(name = element_blank())+
  scale_x_continuous(name = element_blank())+
  # ggtitle(label = "Initial")+
  theme_minimal(base_size = 14)

# Reduced Grey
p_reduced_bioclim <- ggplot()+
  geom_sf(data = st_as_sf(coastline))+
  geom_raster(data = tidy(reduced_clusters), 
              aes(x = Longitude, y= Latitude, fill = as.factor(Cluster)), alpha = 0.8)+
  scale_fill_grey(name = "Bioclim\ncluster")+
  # scale_fill_grey(name = "Bioclim\ncluster", start = 0.9, end = 0)+
  scale_y_continuous(name = element_blank())+
  scale_x_continuous(name = element_blank())+
  # ggtitle(label = "Reduced")+
  theme_minimal(base_size = 14)

# Reduced Colored
p_reduced_bioclim <- ggplot()+
  geom_sf(data = st_as_sf(coastline))+
  geom_raster(data = data_frame(Longitude = coordinates(reduced_clusters)[,1], 
                                Latitude = coordinates(reduced_clusters)[,2], Cluster = reduced_clusters@data$Cluster), 
              aes(x = Longitude, y= Latitude, fill = as.factor(Cluster)), alpha = 0.8)+
  # scale_fill_viridis_d(name = "Bioclim\ncluster", option = "B", direction = -1)+
  scale_fill_viridis_d(name = "Bioclim\ncluster")+
  # scale_fill_grey(name = "Bioclim\ncluster", start = 0.9, end = 0)+
  scale_y_continuous(name = element_blank())+
  scale_x_continuous(name = element_blank())+
  # ggtitle(label = "Reduced")+
  theme_minimal(base_size = 14)

# Save the plots
ggsave(plot = p_initial_bioclim,
       filename = "initial_bioclim_clusters_3.svg", path = "plots/Colored/", 
       device ="svg", dpi = 1200)
ggsave(plot = p_reduced_bioclim,
       filename = "reduced_bioclim_clusters_3.svg", path = "plots/Colored", 
       device ="svg", dpi = 1200)

# 4.2. Actual localitites

p_actual_localities <- ggplot()+
  geom_sf(data = st_as_sf(coastline))+
  geom_point(data = actual_localities, aes(Longitude, Latitude), pch = 1, size = 0.8)+
  scale_y_continuous(name = element_blank())+
  scale_x_continuous(name = element_blank())+
  # ggtitle(label = "Actual localities")+
  ggtitle(label = "A")+
  theme_minimal(base_size = 10)+theme(plot.title = element_text(hjust = 0, size = 8))

# Save the plot
ggsave(plot = p_actual_localities,
       filename = "actual_localities_letters.svg", path = "plots/", 
       device ="svg", dpi = 1200, width = 3, height = 2)

# 4.3. Model predictions

make_letters <- c(Maxnet = "B", plsRglm = "C", LogitBoost = "D",
                  gbm = "E", mlpML = "F", nb = "G", rf = "H", svmRadial = "I")

# 4.3.1. Initial
p_initial_predictions <- ggplot() +  
  geom_sf(data = st_as_sf(coastline))+
  geom_tile(data=init_preds_aggr, aes(x=x, y=y, fill=value), alpha=0.8)+
  scale_fill_viridis_c(name = element_blank(), breaks = seq(1, 6, 1), guide = "legend")+
  # scale_fill_gradient(name = element_blank(),
  #                     low = "#D9D9D9", high = "#252525", guide = "legend",
  #                     limits = c(1, 6),  breaks = 1:6)+
  scale_x_continuous(name = element_blank(), breaks = c(32:36))+
  scale_y_continuous(name = element_blank())+
  # ggtitle(label = "Initial")+
  theme_minimal(base_size = 10)+
  facet_wrap(.~model, nrow = 2, labeller = labeller(model = make_letters))+
  theme(strip.text = element_text(hjust = 0))

# Save the plot
ggsave(plot = p_initial_predictions,
       filename = "initial_predictions.svg", path = "plots/colored/", 
       device ="svg", dpi = 1200, width = 9, height = 5, units = "in")

# Save the plots for every model results separately
for(i in 1:length(model_names)){
  print(i)
  data <- init_preds_aggr %>% 
    filter(model == model_names[i])
  tepm_plot<- ggplot()+
    geom_sf(data = st_as_sf(coastline))+
    geom_tile(data=data, aes(x=x, y=y, fill=value), alpha=0.8)+
    scale_fill_gradient(name = element_blank(),
                        low = "#D9D9D9", high = "#252525", guide = "legend",
                        limits = c(1, 6),  breaks = 1:6)+
    scale_x_continuous(name = element_blank())+
    scale_y_continuous(name = element_blank())+
    # ggtitle(label = "Initial")+
    theme_minimal(base_size = 14)
  # Save the plot
  ggsave(plot = tepm_plot,
         filename = paste0("initial_", model_names[i], ".svg"), path = "plots/", 
         device ="svg", dpi = 1200, width = 5, height = 3, units = "in")
}

# 4.3.2. Reduced
p_reduced_predictions <- ggplot() +  
  geom_sf(data = st_as_sf(coastline))+
  geom_tile(data=reduced_preds_aggr, aes(x=x, y=y, fill=value), alpha=0.8)+
  scale_fill_viridis_c(name = element_blank(), limits = c(1, 5), breaks = seq(1, 5, 1), guide = "legend")+
  # scale_fill_gradient(name = element_blank(),
  #                     low = "#D9D9D9", high = "#252525", guide = "legend",
  #                     limits = c(1, 6),  breaks = 1:6)+
  scale_x_continuous(name = element_blank(), breaks = c(32:36))+
  scale_y_continuous(name = element_blank())+
  # ggtitle(label = "Reduced")+
  theme_minimal()+
  facet_wrap(~model, nrow = 2, labeller = labeller(model = make_letters))+
  theme(strip.text = element_text(hjust = 0))

# Save the plot
ggsave(plot = p_reduced_predictions,
       filename = "reduced_predictions.svg", path = "plots/colored/", 
       device ="svg", dpi = 1200, width = 9, height = 5, units = "in")

# Save the plots for every model results separately
for(i in 1:length(model_names)){
  print(i)
  data <- reduced_preds_aggr %>% 
    filter(model == model_names[i])
  tepm_plot<- ggplot()+
    geom_sf(data = st_as_sf(coastline))+
    geom_tile(data=data, aes(x=x, y=y, fill=value), alpha=0.8)+
    scale_fill_gradient(name = element_blank(),
                        low = "#D9D9D9", high = "#252525", guide = "legend",
                        limits = c(1, 6),  breaks = 1:6)+
    scale_x_continuous(name = element_blank())+
    scale_y_continuous(name = element_blank())+
    # ggtitle(label = "Initial")+
    theme_minimal(base_size = 14)+
    facet_wrap(~model, nrow = 2)
  # Save the plot
  ggsave(plot = tepm_plot,
         filename = paste0("reduced_", model_names[i], ".svg"), path = "plots/", 
         device ="svg", dpi = 1200, width = 5, height = 3, units = "in")
}


# 4.3.3. Combined plots

library(gridExtra)

# Actual localities + initial predictions
initial_plot <- grid.arrange(grobs = list(p_actual_localities, p_initial_predictions), 
                             widths = c(1, 4))

ggsave(plot = initial_plot,
       filename = "Initial.svg", path = "plots/Colored/", 
       device ="svg", dpi = 1200, width = 10, height = 6, units = "in")


# Actual localities + reduced predictions
reduced_plot <- grid.arrange(grobs = list(p_actual_localities, p_reduced_predictions), 
                             widths = c(1, 4))

ggsave(plot = reduced_plot,
       filename = "Reduced.svg", path = "plots/Colored/", 
       device ="svg", dpi = 1200, width = 10, height = 6, units = "in")
