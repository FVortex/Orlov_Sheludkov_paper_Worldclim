###R script to the paper 'Bioclimatic Data Optimization for Spatial Distribution Models'

#Clearance to remove possible objects etc. from previous session
#to avoid interference
rm(list = ls())

###Extract climate data from WorldClim.org tiles for several locations and make data table

#download and unzip all relevant WorldClim geoTIFF files into a single directory.
setwd('wc2.0_2.5m_bio')
#check the folder content
dir()

#load packages: raster, rgdal, foreach, etc.
library(rgdal)
library(raster)
library(foreach) #used in the initial script that I used
library(ggplot2)
library(RColorBrewer)
library(gridExtra)
library(factoextra)
library(FactoMineR)
library(corrplot)
library(maps)
library(mapdata)
library(dendextend)

library(cluster)
library(caret)
library(maxnet) #Maxent implemented in the Maxent package as well;
#But this one doesn't require additional java script and evaluated entirely in R 
#Essentially this implementation is based on glmnet algorithm

#Read names of all .tif (geotif) files in directory into a list
#taken from http://stackoverflow.com/questions/5319839/read-multiple-csv-files-into-separate-data-frames
filenames <- grep('*.tif', list.files(), value = T)

#Load all geoTIFF files
for(i in filenames){
  filepath <- file.path(i)
  assign(i, raster(filepath))
}

#check that all files loaded properly by raster
#taken from http://stackoverflow.com/questions/15387727/use-object-names-as-list-names-in-r
list <- mget(filenames, envir=globalenv())

for(i in list){
  if (hasValues(i)==FALSE){
    print(i,"hasValues error")
  }
  if (inMemory(i)==TRUE){
    print(i, "inMemory error")
  }
  else{
    print("All checked out!")
  }
}

#to determine region of interest
crimea <- map(ylim=c(44.3, 46), xlim=c(32.5,36.6), col='gray90', fill=TRUE)  

#boundaries for which background data will be extracted
crimea_x_lims <- c(32.5, 36.6) 
crimea_y_lims <- c(44.4, 46)

abline(v = crimea_x_lims)
abline(h = crimea_y_lims)

#creating table of latitude and longitude for locations of interest
Latitude <- seq(from = crimea_y_lims[1], to = crimea_y_lims[2], by = 0.062) #approx 2'
Longitude <- seq(from = crimea_x_lims[1], to = crimea_x_lims[2], by = 0.062) #approx 2'
lat_long <- expand.grid(Latitude = Latitude, Longitude = Longitude) #making pairs of latitude and longitude
Pop <- as.factor(1:nrow(lat_long)) #populations ids
pop <- cbind(lat_long, Pop) #adding Ids
row.names(pop) <- pop$Pop
#head(pop)

###NOTE: WorldClim data, even at highest res, is averaged over 1 km2.
#If your location is too close to a coast (i.e. less than 1 km2),
#there will not be any information (nothing but NAs) in this data set.
#Therefore, it may be necessary to approximate some locations by moving
#them inland. I did this using Google Earth.

#load location coordinates as SpatialPoints
for(i in pop$Pop){
  assign(i, SpatialPoints(as.matrix(t(c(pop[i,2], pop[i,1])))))
}

#check that SpatialPoints load correctly from geoTIFFs
poplist <- mget(levels(pop$Pop), envir=globalenv())

tiffvector <- unlist(list)

#Optional quality check step. For smaller datasets, will tell you which population locations should be adjusted,
#in other words, which rows are all NA. See Note above, line 51. Or check after extracting data, see line 
#foreach(p=poplist, .combine='rbind') %:%
# foreach(t=tiffvector, .combine='cbind') %do%{
#  is.na(extract(t,p))
#} #may take a while

###make climate data table

# # # is a plug 

# # #climate <- foreach(p=poplist, .combine='rbind') %:%
# # #  foreach(t=tiffvector, .combine='cbind') %do%{
# # #  myValue<-extract(t, p)
# # #    } #may take a while

# # #save(climate, file = 'climate_complete_crimea_for_Orlov_Sheludkov.rda')
load('climate_complete_crimea_for_Orlov_Sheludkov.rda')
# # #str(climate)
# # #dim(climate)
# # #class(climate)

#tidying table
popnames <- sort(as.character(pop$Pop))
clim <- as.data.frame(climate, row.names=popnames)

#check for NAs
movepops <- clim[rowSums(is.na(clim)) == ncol(clim),] #from the initial script from the Internet
#or
table(complete.cases(clim))['TRUE'] #rows that are not NAs

head(clim)
# # #View(clim)
#how many points there are in total
unique(lapply(clim, length))

#matrix to data frame transformation; probably there's a better way
res <- c()
for (i in clim) {
  res <- cbind(res, i)
}

#setting bioclimatic variable names
raw <- 'BIO1 = Annual Mean Temperature, BIO2 = Mean Diurnal Range (Mean of monthly (max temp - min temp)), BIO3 = Isothermality (BIO2/BIO7) (* 100), BIO4 = Temperature Seasonality (standard deviation *100), BIO5 = Max Temperature of Warmest Month, BIO6 = Min Temperature of Coldest Month, BIO7 = Temperature Annual Range (BIO5-BIO6), BIO8 = Mean Temperature of Wettest Quarter, BIO9 = Mean Temperature of Driest Quarter, BIO10 = Mean Temperature of Warmest Quarter, BIO11 = Mean Temperature of Coldest Quarter, BIO12 = Annual Precipitation, BIO13 = Precipitation of Wettest Month, BIO14 = Precipitation of Driest Month, BIO15 = Precipitation Seasonality (Coefficient of Variation), BIO16 = Precipitation of Wettest Quarter, BIO17 = Precipitation of Driest Quarter, BIO18 = Precipitation of Warmest Quarter, BIO19 = Precipitation of Coldest Quarter'
#raw1 <- gsub('BIO.*? = ', '', raw)
bioclim_vars <- unlist(strsplit(raw, split = ', '))
bioclim_vars <- substr(bioclim_vars, start = 8, 100)

#to restore order
rownames(res) <- as.numeric(popnames)
res <- cbind(pop, res) #population Ids are added
colnames(res)[4:22] <- bioclim_vars

#removing NAs
res_no_nas <- res[complete.cases(res),]
res_no_nas_scaled <- scale(res_no_nas[,-c(1:3)], center = T, scale = T)
initial_res_no_nas_scaled <- res_no_nas_scaled
# # #save(initial_res_no_nas_scaled, file = 'initial_complete_climate_no_nas_scaled.Rda')

# # #View(res)
#plots just to check the data
#visualization of certain variables (colouring by isothermality)
#vector by which the color will be set
vectocol <- (min(sort(unique(na.omit(res_no_nas$`Isothermality (BIO2/BIO7) (* 100)`)))):max(sort(unique(na.omit(res_no_nas$`Isothermality (BIO2/BIO7) (* 100)`)))))
ggplot_Isothermality <- ggplot(data = res_no_nas, mapping = aes(x = res_no_nas$Longitude, 
                                 y = res_no_nas$Latitude, 
                                 lwd = 1, alpha = res_no_nas$`Annual Mean Temperature`, 
                                 color = res_no_nas$`Isothermality (BIO2/BIO7) (* 100)`))+geom_point(show.legend = F) + ggtitle(label = 'Isothermality (BIO 2/BIO 7)')+coord_fixed()

ggplot_Seasonality <- ggplot(data = res_no_nas, mapping = aes(x = res_no_nas$Longitude, 
                                 y = res_no_nas$Latitude, 
                                 lwd = 1, alpha = 1, 
                                 color = res_no_nas$`Temperature Seasonality (standard deviation *100)`))+geom_point(show.legend = F)+ ggtitle(label = 'Temperature Seasonality') + coord_fixed()
grid.arrange(ggplot_Isothermality, ggplot_Seasonality, ncol = 2)

#hierarchical clusterization etc.
hclusted <- hclust(dist(res_no_nas_scaled), method = 'ward.D2')

# # #svg('initial_hclust_no_numbers.svg', height = 4, width = 6)
par(mar = rep(3, 4))
plot(hclusted, labels = F, main = "Ward's method clusterization for\n complete bioclimate data", xaxt = 'n', xlab = '', col = 'grey60', sub ='', ylim = c(-10, 100))
rect.hclust(hclusted, k = 5, border =  1)
text(x = c(25, 125, 305, 440, 825), y = rep(-15,5),  as.roman(1:5), cex = 0.9)
# # #dev.off()

# cutting hclust object
cut_hclusted <- cutree(hclusted, k = 5)

# # #svg('ward_clustering_2_10_clusters.svg', height = 8, width = 8)
lwd= 4
par(mfrow = c(3,3),
    mar = rep(1.3,4))
for (i in 2:10){
  tmp <- cutree(hclusted, k = i)
  plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rainbow(length(unique(tmp)))[tmp], asp = 1, pch = 12, lwd = lwd, main = paste0('Worldclim data in ', i,' clusters (Ward method)'), ylim = crimea_y_lims, xlim = crimea_x_lims)
  lines(x = crimea$x, y = crimea$y)
}
# # #dev.off()

#greyscale ward's clusterization
simpheropole_coord <- c(44.952116, 34.102411)

lwd= 4
par(mfrow = c(3,3),
    mar = rep(2,4))
for (i in 2:10){
  tmp <- cutree(hclusted, k = i)
  plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rev(grey.colors(length(unique(tmp))))[tmp], asp = 1, pch = 12, lwd = lwd, main = paste0('Worldclim data in ', i,' clusters (Ward method)'), ylim = crimea_y_lims, xlim = crimea_x_lims)
  lines(x = crimea$x, y = crimea$y)
  
  points(y = simpheropole_coord[1], x = simpheropole_coord[2], pch = '\u2605', col = 'white', lwd = 1)
  points(y = simpheropole_coord[1]+0.025, x = simpheropole_coord[2]+0.025, pch = '\u2605', col = 'black', lwd = 1)
  
}

#another clusterization method - k-means
# # #svg('kmeans_clustering_2_10_clusters.svg', height = 8, width = 8)
lwd= 4
par(mfrow = c(3,3),
    mar = rep(1.3,4))
for (i in 2:10){
  tmp <- kmeans(res_no_nas_scaled, i)$cluster
  plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rainbow(length(unique(tmp)))[tmp], asp = 1, pch = 12, lwd = lwd, main = paste0('Worldclim data in ', i,' clusters (k-means method)'))
  lines(x = crimea$x, y = crimea$y)
}
# #dev.off()

# third clusterization technique - PAM (partition around medoids)
# # #svg('pam_clustering_2_10_clusters.svg', height = 8, width = 8)
lwd= 4
par(mfrow = c(3,3),
    mar = rep(1.3,4))
for (i in 2:10){
  tmp <- pam(res_no_nas_scaled, i)$cluster
  plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rainbow(length(unique(tmp)))[tmp], asp = 1, pch = 12, lwd = lwd, main = paste0('Worldclim data in ', i,' clusters (PAM method)'))
  lines(x = crimea$x, y = crimea$y)
}
# # #dev.off()


# three clustering methods plotted alongside - allows to select optimal method and number of clusters
# # #svg('3_methods_clustering_2_10_clusters.# #svg', height = 8, width = 24)
lwd= 4
par(mfrow = c(3,9),
    mar = rep(1.3,4))
for (i in 2:10){
  tmp <- cutree(hclusted, k = i)
  plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rainbow(length(unique(tmp)))[tmp], asp = 1, pch = 12, lwd = lwd, main = paste0(i,' clusters (Ward method)'), ylim = crimea_y_lims, xlim = crimea_x_lims)
  lines(x = crimea$x, y = crimea$y)
}
for (i in 2:10){
  tmp <- kmeans(res_no_nas_scaled, i)$cluster
  plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rainbow(length(unique(tmp)))[tmp], asp = 1, pch = 12, lwd = lwd, main = paste0(i,' clusters (k-means method)'))
  lines(x = crimea$x, y = crimea$y)
}
for (i in 2:10){
  tmp <- pam(res_no_nas_scaled, i)$cluster
  plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rainbow(length(unique(tmp)))[tmp], asp = 1, pch = 12, lwd = lwd, main = paste0(i,' clusters (PAM method)'))
  lines(x = crimea$x, y = crimea$y)
} 
# # #dev.off()

# hclust 5 clusters - greyscale and other representation
#cutting dendrogram, i.e. selecting threshold defining clusters
cut_hclusted <- cutree(hclusted, k = 5) #k came from visual dendrogram avealuation
#it also correpsonds to real climatic regions number
plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = grey.colors(length(unique(tmp)))[tmp], asp = 1, pch = 22, lwd = 9, main = "Ward's method clusterization", ylim = crimea_y_lims, xlim = crimea_x_lims, xlab = 'Longitude', ylab = 'Latitude')
lines(x = crimea$x, y = crimea$y)
legend('bottomright', legend = c(paste('Cluster', as.roman(1:5))), pt.bg =  grey.colors(5), pch = 22, bty = 'n')

#other way of plotting
plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rainbow(length(unique(cut_hclusted)))[cut_hclusted], asp = 1, pch = 20, xlab = 'Longitude', ylab = 'Latitude')
lines(x = crimea$x, y = crimea$y)

#trying different approach - scaling data prior to clusterization
#setting vector of colors
colvec <- brewer.pal(5, name = 'Dark2') 
kmeansed <- kmeans(scale(res_no_nas[,-(1:3)], center = T, scale = T), 3)
ggplot(data = res_no_nas, mapping = aes(y= res_no_nas$Latitude, x = res_no_nas$Longitude, lwd = res_no_nas$` Annual Precipitation`, alpha = res_no_nas$`Annual Mean Temperature`, color = colvec[kmeansed$cluster]))+geom_point()+scale_colour_identity()+ coord_fixed()

#PCA (principal component analysis), base R
#Compute PCA in R using prcomp()
res.pca <- prcomp(res_no_nas[,-c(1:3)], scale = TRUE) #no columns corresponding to coordinates and population names

#Visualize eigenvalues (scree plot). Show the percentage of variances explained by each principal component.
fviz_eig(res.pca)
#1 and 2 PCs explain over 70% of variance
#It is therefore reasonable to use 1, 2, 3, but not 4 and so forth PCs

#Graph of individuals. Individuals with a similar profile are grouped together.
colvec <- colorspace::rainbow_hcl(length(unique(cut_hclusted)))[cut_hclusted]
fviz_pca_ind(res.pca,axes = c(1,2),
                      col.var = "contrib", # Variables color
                      #col.ind = rainbow(length(unique(cut_hclusted)))[cut_hclusted], 
                      col.ind = colvec,
                      label = 'var',# Color by the quality of representation
                      #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                      repel = TRUE     # Avoid text overlapping
)


#Graph of variables. Positive correlated variables point to the same side of the plot. Negative correlated variables point to opposite sides of the graph.

fviz_pca_var(res.pca,
             label = 'var', # Individuals color
             col.var = "contrib", # Color by contributions to the PC
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE     # Avoid text overlapping
)

#Biplot of individuals and variables
#(couldn't set colors for dendrogram rigth)
# # #svg('initial_biplot.svg', height = 4, width = 6)
fviz_pca_biplot(res.pca, repel = TRUE,
                col.var = "contrib", # Variables color
                col.ind = 'grey75', 
                fill.ind = 'grey75',
                label = 'var',
                alpha_ind = 2# Individuals color
)+  scale_colour_gradient(low = "grey10", high = "black") +
  theme(legend.position = 'none')

#tripartite plot: bioclamate dataset exploratory analysis
#colors vector for dendrogram; on order to use unified colors
colvec_dendro <- c('#E495A5', '#BDAB66', '#65BC8C', '#55B8D0', '#C29DDE' )
plot(seq_along(colvec_dendro), col = colvec_dendro, lwd =100)
dendro_bioclimate <- hclusted %>% as.dendrogram %>% raise.dendrogram(10) %>% 
  color_branches(k = 5, col = colvec_dendro) %>% 
  set('labels_cex', 1e-9) %>% 
  ggplot() + ggtitle(label = "Ward's method clusterization")
map_bioclimate <- ggplot(data = res_no_nas, mapping = aes(y = res_no_nas$Latitude, 
                                                          x = res_no_nas$Longitude, 
                                                          lwd = 0.5, 
                                                          color = rainbow(length(unique(cut_hclusted)))[cut_hclusted]))+
  geom_point(show.legend = F) + ylab('Latitude') + xlab('Longitude') + 
  ggtitle(label = "Generic map")
pca_bioclimate <- fviz_pca_ind(res.pca,axes = c(1,2),
                               col.var = "contrib", # Variables color
                               #col.ind = rainbow(length(unique(cut_hclusted)))[cut_hclusted], 
                               col.ind = colvec,
                               label = 'var',# Color by the quality of representation
                               #gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
                               repel = TRUE     # Avoid text overlapping
)
grid.arrange(dendro_bioclimate,
             map_bioclimate,
             pca_bioclimate,
             ncol = 1)
# # #dev.off()

#correlation analysis (Pearson correlation coefficient)
#setting shorter names
colnames(res_no_nas) <- paste('BIO', -2:19) 
# # #svg('initial_corrplot.svg', height = 6, width = 6)
corrplot(cor((res_no_nas[,-c(1:3)]), method = 'p'), col= colorRampPalette(c("grey75", "grey15"))(10), tl.col = 1)
# # #dev.off()

# # # # # # PLANT DATA PREPARATION # # # # # # # # # # # # 
#reading in lacalities for multiple plants
species <- readxl::read_xlsx('/home/mikhail/Documents/lsh17/species-points_final_update.xlsx')[,1:3]
par (mar = c(15,3,2,2))
barplot(sort(table(species$Species), decreasing = T)[1:20], las = 2, main = 'Number of localities for species ')
#Pimpinella tragium is clearly stands out
P.tragium <- species[which(species$Species=='Pimpinella tragium Vill.'),]

#all plants representing species number of which equals ...
length(unique(species$Species))
#... are plotted together
map <- map(ylim=c(44.3, 46), xlim=c(32.5,36.6), col='gray90', fill = T)
plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rev(grey.colors(length(unique(tmp))))[tmp], asp = 1, pch = 12, lwd = lwd, main = 'All species localities', ylim = crimea_y_lims, xlim = crimea_x_lims)
lines(x = crimea$x, y = crimea$y)
points(x = species$Longitude, y = species$Latitude, col = as.numeric(as.factor(species$Species)), pch = 16, lwd = 0.2, ylim = crimea_y_lims, xlim = crimea_x_lims) #little trick to color different species with corresponding colors
#note: botanists are evading ridges

# all P. tragium actual localities on the map
# some points might be referenced to the sea by mistake
# there's no need in revoming them since geotiffs contain NAs in that points
plot(x = res_no_nas$Longitude, y = res_no_nas$Latitude, col = rev(grey.colors(length(unique(tmp))))[tmp], asp = 1, pch = 12, lwd = lwd, main = 'P. tragium actual localities', ylim = crimea_y_lims, xlim = crimea_x_lims)
lines(x = crimea$x, y = crimea$y)
points(y = P.tragium$Latitude, x = P.tragium$Longitude, pch = 16, lwd = 5)

#creating smaller clim object for P.tragium
#lattitude, longitude, and population Ids
P.tragium_lat_long<- as.data.frame(cbind(as.numeric(P.tragium$Latitude), as.numeric(P.tragium$Longitude)))
P.tragium_Pop <- as.factor(1:nrow(P.tragium_lat_long))
P.tragium_pop <- cbind((P.tragium_lat_long), P.tragium_Pop)
colnames(P.tragium_pop)[1:2] <- c('Latitude', 'Longitude')
row.names(P.tragium_pop) <- P.tragium_pop$P.tragium_Pop

head(P.tragium_pop)

#load location coordinates as SpatialPoints
for(i in P.tragium_pop$P.tragium_Pop){
  assign(i,SpatialPoints(as.matrix(t(c(P.tragium_pop[i,2], P.tragium_pop[i,1])))))
}

#check that SpatialPoints load correctly from geoTIFFs
P.tragium_poplist <- mget(levels(P.tragium_pop$P.tragium_Pop), envir=globalenv())
tiffvector <- unlist(list)

#make climate data table for P. tragium
# # #P.tragium_climate <- foreach(p=P.tragium_poplist, .combine='rbind') %:%
  # # #    foreach(t=tiffvector, .combine='cbind') %do%{
    # # #    myValue<-extract(t, p)
    # # #} #may take a while

# # #save(P.tragium_climate, file = 'P.tragium_climate_2.5min_to_paper.rda')
load('P.tragium_climate_2.5min_to_paper.rda')

# data tidying
P.tragium_popnames <- sort(as.character(P.tragium_pop$P.tragium_Pop))
# complete.cases function is applied to remove 1 row in P.tragium data that contains NAs
P.tragium_clim <- as.data.frame(P.tragium_climate[complete.cases(P.tragium_climate),], row.names=P.tragium_popnames)
colnames(P.tragium_clim) <- bioclim_vars

# matrix to data frame conversion
P.tragium_res <- c()
for (i in P.tragium_clim) {
  P.tragium_res <- cbind(P.tragium_res, i)
}
class(P.tragium_res); class(P.tragium_clim)

rownames(P.tragium_res) <- as.numeric(P.tragium_popnames[
  -which(!complete.cases(P.tragium_climate))]) # Id of the population data for which contains NAs is also removed
P.tragium_res <- cbind(P.tragium_pop[-which(!complete.cases(P.tragium_climate)),], P.tragium_res) #adding population Ids

colnames(P.tragium_res)[4:22] <- names(P.tragium_clim)

# # #View(P.tragium_res)
dim(P.tragium_res)
# # #save(P.tragium_res, file = 'P.tragium_to_ML.rda')

# # # # # # #ML PART# # # # # #
#to adjoin coordinates and population names with columns of SCALED bioclimatic variables for background
background_to_ML_scaled <- cbind(res_no_nas[,c(1:3)], res_no_nas_scaled)
colnames(background_to_ML_scaled)[1:3] <- c('Latitude', 'Longitude', 'Pop')

# # #background_to_ML_scaled[1:5,1:5]
#to adjoin coordinates and population names with columns of SCALED bioclimatic variables for the plant data
P.tragium_climate_scaled <- cbind(P.tragium_res[,1:3],
                                  scale(P.tragium_clim, center = T, scale = T))
# # #P.tragium_climate_scaled[1:5,1:5]

# column names are to be the same
colnames(background_to_ML_scaled)
colnames(P.tragium_climate_scaled) <- colnames(background_to_ML_scaled)

#preparing one data frame for training and testing
#130 background points are randomly selected and used for training
background_portion <- background_to_ML_scaled[sample(1:nrow(background_to_ML_scaled), 130),]
set.seed(999)
# df contains data for 130 P. tragium points ans 130 background points
df <- as.data.frame(rbind(P.tragium_climate_scaled, 
                          background_portion))
dim(df)

#separate data frame with coordinates which will be used for plotting prediction on maps
df_coords <- as.data.frame(rbind(as.matrix(P.tragium_climate_scaled[,1:2]), background_portion[,1:2]))
dim(df_coords)

#prediction outcomes are denoted; 1 means the presence and 0 maens 'the absence' (no information)
ps <- c(rep(1, nrow(P.tragium_clim)), rep(0, nrow(background_portion)))
length(ps)

#factor columns for outcomes
df$Class <- as.factor(ps)

#alternative to createDataPartition() way of creating subsets for training and testing
#data frames for outcome 0 or 1 only 
df_0s <- df[which(df$Class == 0),-c(1:3)]; dim(df_0s)
df_1s <- df[which(df$Class == 1),-c(1:3)]; dim(df_1s)                  

#separeta data frames for coordinates for outcome 0 or 1 only 
df_0s_crds <- df_coords[which(df$Class == 0),]; dim(df_0s_crds)
df_1s_crds <- df_coords[which(df$Class == 1),]; dim(df_1s_crds)                  

#how many points corresponding to outcome 1 to use for training?
nos <- sample(nrow(df_1s), 60)

#object for training and testing, correspondingly
Train <- rbind(df_0s[nos,],
               df_1s[nos,])
Test <- rbind(df_0s[-nos,],
              df_1s[-nos,])
#coordinates to Test object. There's no need in alike object for Train - we only plot predictions
Test_crds <- rbind(df_0s_crds[-nos,],
                   df_1s_crds[-nos,])

dim(Train)
dim(Test)
dim(Test_crds)
#new outcome vectors

p_train <- c(rep(0, nrow(df_0s[nos,])), rep(1, nrow(df_1s[nos,])))
p_test <- c(rep(0, nrow(df_0s[-nos,])), rep(1, nrow(df_1s[-nos,])))

dim(Train)
dim(Test)

# maxent model training
dim(Train)

Train_to_maxnet <- Train[,-19] #outcomes columns is excluded
colnames(Train_to_maxnet) <- NULL -> rownames(Train_to_maxnet)

#converting into format that peaky maxnet would use

mat <- as.matrix(unname(Train_to_maxnet))
test <- as.matrix(unname(Test[,-19]))
#model training
model <- maxnet(p = p_train, data = as.data.frame(mat),  maxnet.formula(p_train, as.data.frame(mat)))
plot(model, type="cloglog") #some diagnostics
class(model)
print(model)
summary(model)

#creating and modifying test object just for maxnet
test1 <- test
colnames(test1) <- NULL -> rownames(test1)
pred_maxnet <- predict(model, as.data.frame(test1), s = "lambda.min")

#selecting cutoff values
#in such a way that number of prediction is comparable to ones from caret models
which(pred_maxnet > 2)
#getting coordinates for prediction
#maxnet_coords_over05 <- res_no_nas[-c(1:130),1:3][which(pred_maxnet>0.5),1:3]

cutoff <- 2
maxnet_coords_over_cutoff <- as.data.frame(Test_crds[which(pred_maxnet>cutoff),])
colnames(maxnet_coords_over_cutoff) <- c('Latitude', 'Longitude')
#maxnet_coords_over3 <- testing_coords_nonas[which(pred_maxnet>cutoff),]

#specific prediction outcomes vector for maxent
predicted_p_test <- rep(0, length(p_test))
predicted_p_test[which(pred_maxnet>cutoff)] <- 1
#confusionMatrix() for performance assesment
confusionMatrix_maxent <- confusionMatrix(as.factor(predicted_p_test), as.factor(p_test))
#extracting main performance values
accuracy_maxent <- confusionMatrix_maxent$overall['Accuracy']
sensitivity_maxent <- confusionMatrix_maxent$byClass['Sensitivity']
specificity_maxent <- confusionMatrix_maxent$byClass['Specificity']

#training setting using basic caret function
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10
  #  savePredictions = T,
  ## Estimate class probabilities
  # classProbs = TRUE
  ## Evaluate performance using 
  ## the following function
  #summaryFunction = twoClassSummary
)

#caret train in a loop
#methods names
methods_to_loop <- c('plsRglm', 'LogitBoost','gbm', 'mlpML', 'nb', 'rf', 'svmRadial')

inner_accuracies <- c()
confusionMatrix_accuracies <- c()
sensitivities <- c()
specificities <- c()
predicted_coords <- list()

varimps <- list()

for (i in seq_along(methods_to_loop)) {
  tmpFit1 <- caret::train(Class ~ ., data = Train, #removing metadata columns, but no class column
                          method = methods_to_loop[i], 
                          trControl = fitControl,
                          #                      metric = 'ROC',
                          ## This last option is actually one
                          ## for gbm() that passes through
                          verbose = F,
                          preProcess = c('center','scale')
  )
  inner_accuracies <- c(inner_accuracies, mean(tmpFit1$results$Accuracy)) # one that come from Fit object itsels
  #prediction
  tmp2 <- predict(tmpFit1, newdata = Test)
  cm <- confusionMatrix(tmp2, as.factor(p_test))
  confusionMatrix_accuracies <- c(confusionMatrix_accuracies, cm$overall['Accuracy'])
  sensitivities <- c(sensitivities, cm$byClass['Sensitivity'])
  specificities <- c(specificities, cm$byClass['Specificity'])
  assign(paste0(methods_to_loop[i], '_predicted_coords'), Test_crds[which(tmp2 ==1),]) #each prediction as a separate object
  #variables importance extraction
  #varimps[i] <- list(varImp(tmpFit1))
}

#performance as barplots
barplot(inner_accuracies, names.arg = methods_to_loop) # no maxent inner accuracy is awailable
barplot(c(accuracy_maxent, confusionMatrix_accuracies), names.arg = c('Maxent', methods_to_loop))
barplot(c(sensitivity_maxent, sensitivities), names.arg = c('Maxent', methods_to_loop))
barplot(c(specificity_maxent, specificities), names.arg = c('Maxent', methods_to_loop))

#all the performance measures for initial dataset
initial_performance <- cbind(c(accuracy_maxent, confusionMatrix_accuracies),
                             c(sensitivity_maxent, sensitivities),
                             c(specificity_maxent, specificities))
colnames(initial_performance) <- c('Accuracy', 'Sensitivity', 'Specificity')
rownames(initial_performance) <- (c('Maxent', methods_to_loop))

#save(initial_performance, file = '/home/mikhail/Documents/lsh17/worldclim/wc2.0_2.5m_bio/initial_performance.rda')

#making a list with predicted coordinates
caret_predicted_coord <- grep('_predicted_coords', x = ls(), value = T)

#how many points were predicted?
sapply(caret_predicted_coord, function(x) {nrow(get(x[1]))})
pch = 11
lwd = 2

crimea <- map(ylim=c(44.3, 46), xlim=c(32.5,36.6), col='gray90', fill=TRUE)  

crimea_x_lims <- c(32.5, 36.6) 
crimea_y_lims <- c(44.4, 46)

par(mfrow = c(3,3),
    mar = rep(0.75, 4))
plot(x = crimea$x, y = crimea$y, xlim = c(32.5, 36.6), ylim = c(44.4, 46), xlab = 'Longitude', ylab = 'Latitude',type = 'l', asp = 1, main = 'Maxent')
points(x = P.tragium_res$Longitude, y = P.tragium_res$Latitude, col =1, pch = pch, lwd = lwd)
#points(x = maxnet_coords_over05$Longitude, y = maxnet_coords_over05$Latitude, col = adjustcolor("red", alpha=0.15),  pch = pch, lwd = lwd)
points(x = maxnet_coords_over_cutoff$Longitude, y = maxnet_coords_over_cutoff$Latitude, col = 2,  pch = pch, lwd = lwd)
#points(x = maxnet_coords_over2$Longitude, y = maxnet_coords_over2$Latitude, col = adjustcolor("red", alpha=1),  pch = pch, lwd = lwd)
## #dev.off()

for (i in caret_predicted_coord){
  plot(x = background_to_ML_scaled$Longitude, y = background_to_ML_scaled$Latitude, type = 'n', asp = 1, pch = 20, main = (gsub(pattern = 'predicted_coord_*', replacement = '', gsub(pattern = '*Fit1', replacement = '', i))), xlim = crimea_x_lims, ylim = crimea_y_lims)
  lines(x = crimea$x, y = crimea$y)
  points(x = P.tragium_res$Longitude, y = P.tragium_res$Latitude, col =1, pch = pch, lwd = lwd)
  #plot(x = crimea$x, y = crimea$y, xlim = c(32.5, 36.6), ylim = c(44.4, 46), xlab = 'Longitude', ylab = 'Latitude',type = 'l', asp = 1, main = 'Pimpinella tragium distribution')
  points(x = get(i)$Longitude, y = get(i)$Latitude, col = adjustcolor("red", alpha=1),  pch = pch, lwd = lwd)
  ## #dev.off()
}


#saving the file
initial <- list()
initial[1] <- list(maxnet_coords_over_cutoff)
for (i in seq_along(caret_predicted_coord)) {
  initial[i + 1] <- list(get(caret_predicted_coord[i]))
}

names(initial) <- c('Maxnet', methods_to_loop)
#save(initial, file = '/home/mikhail/Documents/lsh17/worldclim/wc2.0_2.5m_bio/initial_dataset_prediction.rda')


#performance barplots

par(mfrow = c(2,8),
    mar = c(2,2,2,.3))
barplot(initial_performance[1,], ylim = 0:1, col = grey.colors(3), las = 2)
abline(h = 0.8, lty = 2)
for (i in 2:nrow(initial_performance)){
  barplot(initial_performance[i,], ylim = 0:1, yaxt = 'n', col = grey.colors(3), las = 2)
  abline(h = 0.8, lty = 2)
}

# # # ML FOR REDUCED DATASET # # # # #
#variables that are less informative (highly correlated / less contributing to the variance) were determine previously
no_to_remove <- c(3, 12:15, 17) #which variables are present in the initial but absent in the reduced # 2 5 7
#removed bioclimatic variables
bioclim_vars[no_to_remove]

#alterniteve list of variables to remove
#+ mean diurnal range, precipitation of driest quarter, temperature annual range, isothermality

####
for (i in 3:8){

  
no_to_remove <- c(2, 3, 7, 12:15, 17)[1:i]
#reduced dataset itself
reduced_background_to_ML_scaled <- background_to_ML_scaled[,-(no_to_remove + 3)] #relative position would be bioclimateic parameter number +3
dim(reduced_background_to_ML_scaled)
colnames(reduced_background_to_ML_scaled)
# checking for variables removed
setdiff(x = colnames(background_to_ML_scaled), y = colnames(reduced_background_to_ML_scaled))
# # # # # # #ML PART# # # # # #
#to adjoin coordinates and population names with columns of SCALED bioclimatic variables for background
colnames(reduced_background_to_ML_scaled)[1:3] <- c('Latitude', 'Longitude', 'Pop')

dim(reduced_background_to_ML_scaled)
# # #reduced_background_to_ML_scaled[1:5,1:5]

#to adjoin coordinates and population names with columns of SCALED bioclimatic variables for the plant data
reduced_P.tragium_climate_scaled <- cbind(P.tragium_climate_scaled[,-(no_to_remove + 3)])#relative position would be bioclimateic parameter number +3
dim(reduced_P.tragium_climate_scaled)
# # #P.tragium_climate_scaled[1:5,1:5]

# column names are to be the same
colnames(reduced_background_to_ML_scaled)
colnames(reduced_P.tragium_climate_scaled) <- colnames(reduced_background_to_ML_scaled)

#preparing one data frame for training and testing
#130 background points are randomly selected and used for training
reduced_background_portion <- reduced_background_to_ML_scaled[sample(1:nrow(reduced_background_to_ML_scaled), 130),]
set.seed(999)
# df contains data for 130 P. tragium points ans 130 background points
reduced_df <- as.data.frame(rbind(reduced_P.tragium_climate_scaled[,-c(1:3)], #no coordinates and population Ids in this data frame...
                          reduced_background_portion[,-c(1:3)])) # ...since they are put into separate one
dim(reduced_df)

#separate data frame with coordinates which will be used for plotting prediction on maps
reduced_df_coords <- as.data.frame(rbind(as.matrix(reduced_P.tragium_climate_scaled[,1:2]), reduced_background_portion[,1:2]))
dim(reduced_df_coords)

#prediction outcomes are denoted; 1 means the presence and 0 maens 'the absence' (no information)
reduced_ps <- c(rep(1, nrow(reduced_P.tragium_climate_scaled)), rep(0, nrow(reduced_background_portion)))
length(reduced_ps)

#factor columns for outcomes
reduced_df$Class <- as.factor(reduced_ps)

#alternative to createDataPartition() way of creating subsets for training and testing
#data frames for outcome 0 or 1 only 
reduced_df_0s <- reduced_df[which(reduced_df$Class == 0),]; dim(reduced_df_0s)
reduced_df_1s <- reduced_df[which(reduced_df$Class == 1),]; dim(reduced_df_1s)                  

#separeta data frames for coordinates for outcome 0 or 1 only 
reduced_df_0s_crds <- reduced_df_coords[which(reduced_df$Class == 0),]; dim(reduced_df_0s_crds)
reduced_df_1s_crds <- reduced_df_coords[which(reduced_df$Class == 1),]; dim(reduced_df_1s_crds)                  

#how many points corresponding to outcome 1 to use for training?
nos <- sample(nrow(reduced_df_1s), 60)

#object for training and testing, correspondingly
reduced_Train <- rbind(reduced_df_0s[nos,],
                       reduced_df_1s[nos,])
reduced_Test <- rbind(reduced_df_0s[-nos,],
                      reduced_df_1s[-nos,])
#coordinates to Test object. There's no need in alike object for Train - we only plot predictions
reduced_Test_crds <- rbind(reduced_df_0s_crds[-nos,],
                           reduced_df_1s_crds[-nos,])

dim(reduced_Train)
dim(reduced_Test)
dim(reduced_Test_crds)
#new outcome vectors

reduced_p_train <- c(rep(0, nrow(reduced_df_0s[nos,])), rep(1, nrow(reduced_df_1s[nos,])))
reduced_p_test <- c(rep(0, nrow(reduced_df_0s[-nos,])), rep(1, nrow(reduced_df_1s[-nos,])))

dim(reduced_Train)
dim(reduced_Test)

# maxent model training
dim(reduced_Train)

reduced_Train_to_maxnet <- reduced_Train[,-ncol(reduced_Train)] #outcomes columns is excluded
colnames(reduced_Train_to_maxnet) <- NULL -> rownames(reduced_Train_to_maxnet)

#converting into format that peaky maxnet would use

reduced_mat <- as.matrix(unname(reduced_Train_to_maxnet))
reduced_test <- as.matrix(unname(reduced_Test[,-ncol(reduced_Test)]))
#model training
reduced_model <- maxnet(p = reduced_p_train, data = as.data.frame(reduced_mat),  maxnet.formula(reduced_p_train, as.data.frame(reduced_mat)))
plot(reduced_model, type="cloglog") #some diagnostics
class(reduced_model)
print(reduced_model)
summary(reduced_model)

#creating and modifying test object just for maxnet
reduced_test1 <- reduced_test
colnames(reduced_test1) <- NULL -> rownames(reduced_test1)
reduced_pred_maxnet <- predict(reduced_model, as.data.frame(reduced_test1), s = "lambda.min")

#selecting cutoff values
#in such a way that number of prediction is comparable to ones from caret models
which(reduced_pred_maxnet > 2)
#getting coordinates for prediction
#maxnet_coords_over05 <- res_no_nas[-c(1:130),1:3][which(pred_maxnet>0.5),1:3]

cutoff <- 2
reduced_maxnet_coords_over_cutoff <- as.data.frame(reduced_Test_crds[which(reduced_pred_maxnet>cutoff),])
colnames(reduced_maxnet_coords_over_cutoff) <- c('Latitude', 'Longitude')
#maxnet_coords_over3 <- testing_coords_nonas[which(pred_maxnet>cutoff),]

#specific prediction outcomes vector for maxent
reduced_predicted_p_test <- rep(0, length(reduced_p_test))
reduced_predicted_p_test[which(reduced_pred_maxnet>cutoff)] <- 1
#confusionMatrix() for performance assesment
reduced_confusionMatrix_maxent <- confusionMatrix(as.factor(reduced_predicted_p_test), as.factor(reduced_p_test))
#extracting main performance values
reduced_accuracy_maxent <- reduced_confusionMatrix_maxent$overall['Accuracy']
reduced_sensitivity_maxent <- reduced_confusionMatrix_maxent$byClass['Sensitivity']
reduced_specificity_maxent <- reduced_confusionMatrix_maxent$byClass['Specificity']

#training setting using basic caret function
fitControl <- trainControl(## 10-fold CV
  method = "repeatedcv",
  number = 10,
  ## repeated ten times
  repeats = 10
  #  savePredictions = T,
  ## Estimate class probabilities
  # classProbs = TRUE
  ## Evaluate performance using 
  ## the following function
  #summaryFunction = twoClassSummary
)

#caret train in a loop
#methods names
methods_to_loop <- c('plsRglm', 'LogitBoost','gbm', 'mlpML', 'nb', 'rf', 'svmRadial')

reduced_inner_accuracies <- c()
reduced_confusionMatrix_accuracies <- c()
reduced_sensitivities <- c()
reduced_specificities <- c()
reduced_predicted_coords <- list()

reduced_varimps <- list()

for (i in seq_along(methods_to_loop)) {
  tmpFit1 <- caret::train(Class ~ ., data = reduced_Train, #removing metadata columns, but no class column
                          method = methods_to_loop[i], 
                          trControl = fitControl,
                          #                      metric = 'ROC',
                          ## This last option is actually one
                          ## for gbm() that passes through
                          verbose = F,
                          preProcess = c('center','scale')
  )
  reduced_inner_accuracies <- c(reduced_inner_accuracies, mean(tmpFit1$results$Accuracy)) # one that come from Fit object itsels
  #prediction
  tmp2 <- predict(tmpFit1, newdata = reduced_Test)
  cm <- confusionMatrix(tmp2, as.factor(reduced_p_test))
  reduced_confusionMatrix_accuracies <- c(reduced_confusionMatrix_accuracies, cm$overall['Accuracy'])
  reduced_sensitivities <- c(reduced_sensitivities, cm$byClass['Sensitivity'])
  reduced_specificities <- c(reduced_specificities, cm$byClass['Specificity'])
  assign(paste0(methods_to_loop[i], '_reduced_predicted_coords'), reduced_Test_crds[which(tmp2 ==1),]) #each prediction as a separate object
  #variables importance extraction
  #reduced_varimps[i] <- list(varImp(tmpFit1))
}

#performance as barplots
barplot(reduced_inner_accuracies, names.arg = methods_to_loop) # no maxent inner accuracy is awailable
barplot(c(reduced_accuracy_maxent, reduced_confusionMatrix_accuracies), names.arg = c('Maxent', methods_to_loop))
barplot(c(reduced_sensitivity_maxent, reduced_sensitivities), names.arg = c('Maxent', methods_to_loop))
barplot(c(reduced_specificity_maxent, reduced_specificities), names.arg = c('Maxent', methods_to_loop))

#all the performance measures for initial dataset
reduced_performance <- cbind(c(reduced_accuracy_maxent, reduced_confusionMatrix_accuracies),
                             c(reduced_sensitivity_maxent, reduced_sensitivities),
                             c(reduced_specificity_maxent, reduced_specificities))
colnames(reduced_performance) <- c('Accuracy', 'Sensitivity', 'Specificity')
rownames(reduced_performance) <- (c('Maxent', methods_to_loop))

save(reduced_performance, file = '/home/mikhail/Documents/lsh17/worldclim/wc2.0_2.5m_bio/reduced_performance1.rda')

#making a list with predicted coordinates
reduced_caret_predicted_coord <- grep('reduced_predicted_coords', x = ls(), value = T)

#how many points were predicted?
sapply(reduced_caret_predicted_coord, function(x) {nrow(get(x[1]))})
pch = 11
lwd = 2

crimea <- map(ylim=c(44.3, 46), xlim=c(32.5,36.6), col='gray90', fill=TRUE)  

crimea_x_lims <- c(32.5, 36.6) 
crimea_y_lims <- c(44.4, 46)

par(mfrow = c(3,3),
    mar = rep(0.75, 4))
plot(x = crimea$x, y = crimea$y, xlim = c(32.5, 36.6), ylim = c(44.4, 46), xlab = 'Longitude', ylab = 'Latitude',type = 'l', asp = 1, main = 'Maxent')
points(x = P.tragium_res$Longitude, y = P.tragium_res$Latitude, col =1, pch = pch, lwd = lwd)
#points(x = maxnet_coords_over05$Longitude, y = maxnet_coords_over05$Latitude, col = adjustcolor("red", alpha=0.15),  pch = pch, lwd = lwd)
points(x = reduced_maxnet_coords_over_cutoff$Longitude, y = reduced_maxnet_coords_over_cutoff$Latitude, col = 2,  pch = pch, lwd = lwd)
#points(x = maxnet_coords_over2$Longitude, y = maxnet_coords_over2$Latitude, col = adjustcolor("red", alpha=1),  pch = pch, lwd = lwd)
## #dev.off()

for (i in reduced_caret_predicted_coord){
  plot(x = reduced_background_to_ML_scaled$Longitude, y = reduced_background_to_ML_scaled$Latitude, type = 'n', asp = 1, pch = 20, main = (gsub(pattern = 'predicted_coord_*', replacement = '', gsub(pattern = '*Fit1', replacement = '', i))), xlim = crimea_x_lims, ylim = crimea_y_lims)
  lines(x = crimea$x, y = crimea$y)
  points(x = P.tragium_res$Longitude, y = P.tragium_res$Latitude, col =1, pch = pch, lwd = lwd)
  #plot(x = crimea$x, y = crimea$y, xlim = c(32.5, 36.6), ylim = c(44.4, 46), xlab = 'Longitude', ylab = 'Latitude',type = 'l', asp = 1, main = 'Pimpinella tragium distribution')
  points(x = get(i)$Longitude, y = get(i)$Latitude, col = adjustcolor("red", alpha=1),  pch = pch, lwd = lwd)
  ## #dev.off()
}

#saving the file
reduced <- list()
reduced[1] <- list(reduced_maxnet_coords_over_cutoff)
for (i in seq_along(reduced_caret_predicted_coord)) {
  reduced[i + 1] <- list(get(reduced_caret_predicted_coord[i]))
}

names(reduced) <- c('Maxnet', methods_to_loop)
#save(initial, file = '/home/mikhail/Documents/lsh17/worldclim/wc2.0_2.5m_bio/initial_dataset_prediction.rda')


#performance barplots

par(mfrow = c(2,8),
    mar = c(2,2,2,.3))
barplot(initial_performance[1,], ylim = 0:1, col = terrain.colors(3), las = 2)
abline(h = 0.8, lty = 2)
for (i in 2:nrow(initial_performance)){
  barplot(initial_performance[i,], ylim = 0:1, yaxt = 'n', col = terrain.colors(3), las = 2)
  abline(h = 0.8, lty = 2)
}

barplot(reduced_performance[1,], ylim = 0:1, col = heat.colors(3), las = 2)
abline(h = 0.8, lty = 2)
for (i in 2:nrow(reduced_performance)){
  barplot(reduced_performance[i,], ylim = 0:1, yaxt = 'n', col = heat.colors(3), las = 2)
  abline(h = 0.8, lty = 2)
}

# both perforamnce sets

#cbind(initial_performance, reduced_performance)
perfomances_diff <- reduced_performance - initial_performance
save(perfomances_diff, file = paste0('/home/mikhail/Documents/Papers, presentations/worldclim_paper/perfomances_diff_upon_removing_', paste(no_to_remove, collapse = '_'), '_variables.rda'))

}



#loading performance from the big loop
setwd('/home/mikhail/Documents/Papers, presentations/worldclim_paper/')
for (i in dir()[grepl(pattern = 'perfomances_', dir())]){
  tmp <- get(load(i))
  assign(i, tmp)
}

performances_overnight <- lapply(grep('perfomances_diff', ls(), value = T), FUN = get)

names(performances_overnight) <-   dir()[grepl(pattern = 'perfomances_', dir())]

#sapply(performances_overnight, barplot)

accuracies_changes <- c()
sensitivties_changes <- c()
specificities_changes <- c()

for (i in rownames(performances_overnight[[1]])){
  tmp_a_method1 <- c()
  tmp_a_method2 <- c()
  tmp_a_method3 <- c()
  for (j in performances_overnight){
    tmp_a_method1 <- c(tmp_a_method1, j[i,'Accuracy'])
    tmp_a_method2 <- c(tmp_a_method2, j[i,'Sensitivity'])
    tmp_a_method3 <- c(tmp_a_method3, j[i,'Specificity'])
  }
  accuracies_changes <- cbind(accuracies_changes, tmp_a_method1)
  sensitivties_changes <- cbind(sensitivties_changes, tmp_a_method2)
  specificities_changes <- cbind(specificities_changes, tmp_a_method3)
}

xlab_raw <- dir()[grepl(pattern = 'perfomances_', dir())]

as.numeric(gsub("([0-9]+).*", "\\1", xlab_raw))
xlab <- sapply(regmatches(xlab_raw, gregexpr("[[:digit:]]+", xlab_raw)), FUN = function(x) {paste0(x, collapse = ', ')})
colnames(accuracies_changes) <- rownames(performances_overnight[[1]])


par(mfrow = c(3,1))
#AC
matplot(accuracies_changes, type = 'l', xaxt = "n")
#axis(1, at=seq_along(xlab), labels=xlab, srt = 45, adj = 1)
text(x = seq_along(xlab), y = -0.12, srt = 45, adj = 1, labels = xlab, xpd = T)
legend('topright', legend = rownames(performances_overnight[[1]]), lty = 1, col = seq_along(xlab))
#SE
matplot(sensitivties_changes, type = 'l', xaxt = "n")
#axis(1, at=seq_along(xlab), labels=xlab, srt = 45, adj = 1)
text(x = seq_along(xlab), y = -0.12, srt = 45, adj = 1, labels = xlab, xpd = T)
legend('topright', legend = rownames(performances_overnight[[1]]), lty = 1, col = seq_along(xlab))
#SP
matplot(specificities_changes, type = 'l', xaxt = "n")
#axis(1, at=seq_along(xlab), labels=xlab, srt = 45, adj = 1)
text(x = seq_along(xlab), y = -0.12, srt = 45, adj = 1, labels = xlab, xpd = T)
legend('topright', legend = rownames(performances_overnight[[1]]), lty = 1, col = seq_along(xlab))

plot(1, type="n", xlab="", ylab="", xlim=c(0, 10), ylim=c(0, 10))
  text(x = 5, y =seq(2,11, by = 0.5), cbind(seq_along(bioclim_vars), bioclim_vars))
