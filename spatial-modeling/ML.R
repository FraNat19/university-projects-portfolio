# HOMEWORK 4 - GROUP 12 (Agate, Bartocci, Natali)


# Part a - Year 2002 ------------------------------------------------------


# Using the shrimps data and the grid of the corresponding year:
# a) Estimate the total biomass on the estimation grid using 
# kriging in a maximum likelihood framework (choose variogram and trend, use the available covariates)
library(reshape2)
library(corrplot)
library(ade4)
library(gstat)
library(sp)
library(geoR)
library(dplyr)
library(akima)
library(sf)
library(ggplot2)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# Load the dataset and the file containing the estimation grids
load("C:/Users/frank/Downloads/shrimpsfull.RData")
load("C:/Users/frank/Downloads/AllGrids.RData")

# We now filter data for our years of interest, i.e. 2002 and 2008
shrimp_data_2002 <- shrimpsdata[shrimpsdata$ANNO == 2002, ]
shrimp_data_2008 <- shrimpsdata[shrimpsdata$ANNO == 2008, ]

colSums(is.na(shrimp_data_2002))
colSums(is.na(shrimp_data_2008))

# Display the structure of the dataset
str(shrimpsdata)
# Preview the first few rows
head(shrimpsdata)

#2002---------------------------------------------------------------

# We now start our analysis for the year 2002
summary(shrimp_data_2002)

#At first, we explore data with a multivariate analysis
pca1 <- dudi.pca(df = shrimp_data_2002[, -c(3:5)], scannf = FALSE, nf = 1)

# Perform PCA with at least two components
pca1 <- dudi.pca(df = shrimp_data_2002[, -c(3:5)], scannf = FALSE, nf = 2)

# Plot the PCA with the first two components, with the scatter() function plotting the first two axes
scatter(pca1, xax = 1, yax = 2)
scatter(pca1, xax = 1, yax = 2, clab.row = 0.7, clab.col = 0.7) #for better reading
#xax = 1, yax = 2 specifies that the first and second principal components (cs1 and cs2) should be plotted.


# Each principal component (PCs) represents a "dimension" of variability in the dataset, 
# with the first component (CS1) capturing the most variance, and the second component (CS2) capturing the second most.
# The values shown under CS1 and CS2 represent the loadings (or coefficients) of each variable 
# on the respective principal component, and they indicate how much each original variable 
# contributes to the component and, indirectly, to the data's main sources of variability.

# We can now visualize the loadings of the covariates with the principal components
print(pca1$c1)

# As we can notice, the loadings on CS1 are dominated by the salinity variables, 
# specifically salinity.minq3p, salinity.maxq3, and salinity.minq3. 
# High loadings on CS1 for these variables indicate that this component 
# likely represents spatial and environmental variability associated with salinity 
# and geographical position (X and Y coordinates).
# Since salinity and location are relatively stable environmental factors, 
# CS1 may reflect a spatial and stable environmental gradient in the dataset.

# CS2 on the other hand has high loadings on temp.maxq3p, temp.maxq3, and temp.minq2,
# as well as on bathymetry (bat) and distance from the coast (dist). 
# This suggests that CS2 represents variability more related to temperature and physical habitat features, 
# This second component may therefore represent a temperature gradient and other physical factors 
# that vary with depth and proximity to land, which could influence the shrimp total biomass distribution.


# Queste sono le variate con i loadings + alti nelle due componenti

# cs1
shrimp_data_2002$salinity.minq3 <- factor(shrimp_data_2002$salinity.minq3)
s.class(pca1$li, fac = shrimp_data_2002$salinity.minq3)

#cs2
shrimp_data_2002$temp.maxq3 <- factor(shrimp_data_2002$temp.maxq3)
s.class(pca1$li, fac = shrimp_data_2002$temp.maxq3)

# Visualizza s.class con una variabile categoriale, come ad esempio "Quarter" 
# (modifica secondo la variabile disponibile)
shrimp_data_2002$bat <- factor(shrimp_data_2002$bat)
s.class(pca1$li, fac = shrimp_data_2002$bat)

shrimp_data_2002$dist <- factor(shrimp_data_2002$dist)
s.class(pca1$li, fac = shrimp_data_2002$dist)

shrimp_data_2002$slope <- factor(shrimp_data_2002$slope)
s.class(pca1$li, fac = shrimp_data_2002$slope)


#2002
# Our dataset does not follow a normal distribution, and therefore, a transformation is necessary. 
#Shrimp biomass is likely to follow a log-normal distribution, a common pattern in ecological studies. 
#In such distributions, most values are concentrated in the lower range, with fewer extreme values observed at the higher end. 
#This skewed pattern is typical for population sizes in ecological contexts.
shrimp_data_2002$log_tot <- log(shrimp_data_2002$tot + 1)

shrimp_geodata_2002_log <- as.geodata(shrimp_data_2002,
                                      coords.col = c("X", "Y"),
                                      data.col = "log_tot",
                                      covar.col = c("bat","salinity.minq3","temp.maxq3","dist"))

# Now, based on this log transformation, we can plot our empirical variogram.
# Before doing that, we have also to choose our trend in order to have a stationary and isotropic variogram,
# which is the basis for the kriging;
# first we plot the variogram with first order trend

variogram_tot_2002 <- variog(shrimp_geodata_2002_log, trend = "1st", max.dist = 80)

plot(variogram_tot_2002,type="b")

# and we include our chosen covariates in the trend
trend_matrix_2002 <- cbind(shrimp_data_2002$bat, 
                           shrimp_data_2002$salinity.minq3, 
                           shrimp_data_2002$temp.maxq3,
                           shrimp_data_2002$dist)

variogram_tot_2002 <- variog(shrimp_geodata_2002_log, trend = "1st",trend.d = trend_matrix_2002, max.dist = 80)

# We now plot the variogram with second order trend
variogram_tot2_2002 <- variog(shrimp_geodata_2002_log, trend = "2nd", max.dist = 80)

plot(variogram_tot2_2002,type="b")

# and we include our chosen covariates in the trend
trend_matrix_2002 <- cbind(shrimp_data_2002$bat, 
                           shrimp_data_2002$salinity.minq3, 
                           shrimp_data_2002$temp.maxq3,
                           shrimp_data_2002$dist)

variogram_tot2_2002 <- variog(shrimp_geodata_2002_log, trend = "2nd",trend.d = trend_matrix_2002, max.dist = 80)

#If we now compare the two variograms, we can notice how their pverall behaviour is pretty similar and that 
# they have equal nugget effect;
# however, the second order trend has a slightly lower sill (which is still reached at distance 25 ca.)
# which might be indicating an overparametrization of the model.
# Therefore, we leave the first order trend as our preferred choice.


# For the choice of the better model for estimation, we can display the function eyefit
eyefit(variogram_tot_2002)

# or apply the function "variofit" for the models that look like a better fit for the variogrm
#fit_variogram <- variofit(variogram_tot_2002, cov.model = "matern", ini.cov.pars = c(5.16, 6.78), kappa=0.2)
#fit_variogram <- variofit(variogram_tot_2002, cov.model = "exponential", ini.cov.pars = c(5.37, 8.31))
#fit_variogram <- variofit(variogram_tot_2002, cov.model = "spherical", ini.cov.pars = c(4.94, 8.31))
#print(fit_variogram)

# The three models that give a better fit overall are the matern (with parameter k = 0.2), the exponential
# (i.e. the particular case of a matern with k = 0.5) and the spherical model
fit_matern_lik_2002 <- likfit(shrimp_geodata_2002_log, 
                               cov.model = "matern", 
                               ini.cov.pars = c(5.16, 6.78), 
                               nugget = 2, kappa = 0.2)
fit_exponential_lik_2002 <- likfit(shrimp_geodata_2002_log, 
                              cov.model = "exponential", 
                              ini.cov.pars = c(6.9, 4.6), 
                              nugget = 2)
fit_spherical_lik_2002 <- likfit(shrimp_geodata_2002_log, 
                                   cov.model = "spherical", 
                                   ini.cov.pars = c(6.8, 12.3), 
                                   nugget = 2)
# we  now plot our empirical variogram with all the three models for a better graphical comparison
plot(variogram_tot_2002, type = "b", xlab = 'Distance', ylab = 'Semivariance')

lines(fit_matern_lik_2002, col = "red", lwd = 2)
lines(fit_exponential_lik_2002, col = "green", lwd = 2)
lines(fit_spherical_lik_2002, col = "blue", lwd = 2)

# At first sight, the matern looks like the better fit overall; 
# we now compute the RMSE of each
vv.mat.2002<-xvalid(shrimp_geodata_2002_log,model=fit_matern_lik_2002)
vv.exp.2002<-xvalid(shrimp_geodata_2002_log,model=fit_exponential_lik_2002)
vv.sph.2002<-xvalid(shrimp_geodata_2002_log,model=fit_spherical_lik_2002)

# Calculate the Mean Squared Error (MSE) for each model
MSE_mat_2002 <- mean(vv.mat.2002$std.error^2)
MSE_exp_2002 <- mean(vv.exp.2002$std.error^2)
MSE_sph_2002 <- mean(vv.sph.2002$std.error^2)

# Calculate the Root Mean Squared Error (RMSE)
RMSE_mat_2002 <- sqrt(MSE_mat_2002)
RMSE_exp_2002 <- sqrt(MSE_exp_2002)
RMSE_sph_2002 <- sqrt(MSE_sph_2002)

# Print the RMSE values
cat("Matern Model: RMSE =", RMSE_mat_2002, "\n")
cat("Exponential Model: RMSE =", RMSE_exp_2002, "\n")
cat("Spherical Model: RMSE =", RMSE_sph_2002, "\n")

# as expected, the matern model is the one with the lowest RMSE and will therefore be our chosen model
# (with first order trend)


# We can now run kriging interpolation for the estimation of total biomass
# We include again the chosen covariates in the trend
trend.d.2002 <- trend.spatial(~ bat + salinity.minq3 + temp.maxq3 + dist, geodata = shrimp_geodata_2002_log)
trend.l.2002 <- trend.spatial(~ grid_2002$bat + grid_2002$salinity.minq3 + grid_2002$temp.maxq3 + grid_2002$dist)

# We control for the kriging
krige.2002 <- krige.control(
  cov.model = "matern",
  cov.pars = c(5.16, 6.78), 
  nugget = 2,
  trend.d = trend.d.2002,
  trend.l = trend.l.2002
)

# Run krige.conv with the updated locations
krig_2002 <- krige.conv(shrimp_geodata_2002_log, locations = as.matrix(grid_2002[, c("X", "Y")]), krige = krige.2002)

# And finally, we can visualize the kriging result for the matern
cc_mat_2002 <- data.frame(X = grid_2002$X, Y = grid_2002$Y, Z = krig_2002$predict)

ggplot(cc_mat_2002, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(
    title = "Kriging Interpolation of Total Shrimp Biomass (2002)",
    x = "X Coordinate",
    y = "Y Coordinate",
    fill = "Biomass Prediction") +
  theme_minimal()

#In 2002, the total shrimp biomass appears to be more concentrated in the central parts of the study area (Lazio e Toscana coast). 
#Brigest areas represent regions with the highest levels of biomass, while darkest areas indicate 
#regions with lower levels. High biomass values (ranging from 4.0 to 6.0) are observed in distinct 
#zones, possibly suggesting favorable environmental conditions that supported dense shrimp populations. 
#The coastal regions, especially in the middle and upper parts of the study area, have a gradient 
#where biomass decreases as you move offshore. It is possible to see that the areas of low biomass 
#intensity are in particular along the Liguria coast, this could be due to a combination of unsuitable features, 
#like the fact that the Ligurian coast is characterized by a steep continental shelf, leading to a rapid increase in depth offshore. 
#Shrimp species typically prefer shallower and more gradually sloping habitats where nutrient-rich waters are more abundant. 
#The steep bathymetric gradient in Liguria may reduce the extent of suitable habitats for shrimp. 
#The opposite situation can be noted in the coast of Tuscany region tends to have more stable and favorable temperature and salinity profiles compared to northern areas like Liguria. 
#These conditions are essential for shrimp breeding and larval development. Stable environmental conditions promote high recruitment and the maintenance of large shrimp populations


# Now, in order to obtain a better and more precise geographical representation,
# we plot the map of Italy to see how the shrimp biomass
# is distributed along the areas of the Tirrenean Sea
# Transform Italy to UTM (zone 32N, EPSG:32632) for consistency with your coordinates
italy <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf")
italy <- st_transform(italy, crs = 32632)

# Convert kriging results to a grid for ggplot
krig_result_df_02 <- data.frame(
  X = grid_2002$X * 1000, # Assuming the grid units need conversion to match Italy's CRS
  Y = grid_2002$Y * 1000,
  Z = krig_2002$predict
)



# Define the bounds for the area of interest based on your prediction data 
x_min <- min(krig_result_df_02$X) - 10000  # Adjust to give some margin 
x_max <- max(krig_result_df_02$X) + 10000 
y_min <- min(krig_result_df_02$Y) - 10000 
y_max <- max(krig_result_df_02$Y) + 10000 

ggplot() + # Plot the kriging predictions as filled contours 
  geom_raster( data =krig_result_df_02, aes(x = X, y = Y, fill = Z) ) + # Add a color scale for kriging predictions 
  scale_fill_viridis_c( option ="viridis", name = "Biomass Prediction" ) + # Overlay Italy's coastline 
  geom_sf( data = italy, fill = NA, color = "black",lwd = 0.7 ) + # Zoom into the area of interest 
  coord_sf( xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand =FALSE ) + # Add labels and theme 
  labs( title = "Interpolated Shrimp Biomass (2002)", x = "Longitude (km)", y ="Latitude (km)" ) + 
  theme_minimal()




# Assuming shrimp_geodata_2002_log contains the original coordinates and observed log-transformed biomass
observed_points_df02 <- data.frame(
  X = shrimp_geodata_2002_log$coords[, 1] * 1000,  # Convert to match UTM coordinates
  Y = shrimp_geodata_2002_log$coords[, 2] * 1000,  # Convert to match UTM coordinates
  Biomass = shrimp_geodata_2002_log$data           # Log-transformed biomass
)

# Plot with Italy map, kriging results, and observed points
ggplot() +
  # Kriging predictions as a raster layer
  geom_raster(data = krig_result_df_02, aes(x = X, y = Y, fill = Z)) +
  # Color scale for kriging predictions
  scale_fill_viridis_c(option = "viridis", name = "Biomass Prediction") +
  # Overlay Italy's coastline
  geom_sf(data = italy, fill = NA, color = "black", lwd = 0.7) +
  # Add observed biomass points, with size proportional to the biomass
  geom_point(data = observed_points_df02, aes(x = X, y = Y, size = Biomass), color = "red", alpha = 0.6) +
  # Size scale for observed points
  scale_size_continuous(name = "Observed Biomass", range = c(1, 5)) +  # Adjust size range as needed
  # Zoom into the area of interest
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE) +
  # Labels and theme
  labs(title = "Interpolated Shrimp Biomass (2002)", x = "Longitude (km)", y = "Latitude (km)") +
  theme_minimal()


#The figure provides a comprehensive visualization of the spatial distribution of shrimp biomass for the year 2002, 
#derived through kriging interpolation. The distribution of the red points indicates that higher observed biomass values 
#generally coincide with areas of high predicted biomass. This overlap suggests that the kriging model has reasonably captured the 
#underlying spatial variability in shrimp biomass. However, there are some discrepancies where observed biomass values do not align perfectly with the predictions. 
#These discrepancies may be due to limitations in the available data or model parameters

# PART B ------------------------------------------------------------------

# b) given a), build pointwise confidence intervals.

# c) map results including confidence intervals, and comment output.

#confidence intervals
low_2002 <- krig_2002$predict- 1.96*sqrt(krig_2002$krige.var)
up_2002 <- krig_2002$predict + 1.96*sqrt(krig_2002$krige.var)

hist(low_2002)
hist(up_2002)
hist(up_2002-low_2002)

cc_lower_2002 <- data.frame(X = grid_2002$X, Y = grid_2002$Y, Z = low_2002)
cc_upper_2002 <- data.frame(X = grid_2002$X, Y = grid_2002$Y, Z = up_2002)

ggplot(cc_lower_2002, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(title = "Lower 95% Confidence Bound 2002") +
  theme_minimal()
#This map illustrates the spatial variation in the lower 95% confidence bound for total biomass in 2002. The color gradient transitions from dark purple (representing lower values, between -5 and -1)
#to bright yellow (indicating higher values, between 1 and 1.5).
#Regions shaded in green to yellow mark areas with relatively higher conservative biomass estimates, suggesting these zones may serve as biomass hotspots. 
#These favorable areas likely correspond to optimal environmental conditions such as suitable bathymetry and salinity levels. In contrast, darker purple to blue regions denote locations with consistently lower biomass estimates,
#implying less favorable conditions, potentially due to factors like deeper water, suboptimal salinity, or nutrient scarcity.
#Overall, this spatial representation along the Gaeta-Genova coast highlights the variability in biomass distribution, emphasizing regions with robust conservative estimates and identifying areas where environmental conditions might limit biomass growth
ggplot(cc_upper_2002, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(title = "Upper 95% Confidence Bound 2002") +
  theme_minimal()
#The Upper 95% Confidence Bound map presents an optimistic projection of potential shrimp density (log scale) across spatial locations. The lighter shades, ranging from yellow to light green, denote regions with higher upper bounds (7–10 log scale), 
#indicating areas with the potential for significant total biomass under ideal environmental conditions. These zones are likely associated with optimal habitats, including shallower coastal shelves or areas near estuaries where depth and salinity are favorable. In contrast, darker shades, 
#from purple to blue, represent areas with lower upper bounds (4–6 log scale), suggesting limited potential for high biomass even under favorable conditions.
#Such areas may correspond to deeper waters, less suitable salinity gradients, or generally suboptimal habitats for shrimp.
#This map offers valuable insights into spatial variations in biomass potential along the Gaeta-Genova coastline, emphasizing areas with higher growth opportunities and identifying regions with limited biomass expectations

interval_width_2002 <- up_2002- low_2002
cc_interval <- data.frame(X = grid_2002$X, Y = grid_2002$Y, Z = interval_width_2002)
ggplot(cc_interval, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(title = "Width of 95% Confidence Interval (2002)") +
  theme_minimal()
#This map illustrates the biomass uncertainty along the coastal stretch from Gaeta to Genova in the northwestern Mediterranean during 2002. The width of the 95% confidence interval serves as a measure of the variability or uncertainty associated with the biomass estimates. Regions highlighted in yellow, with widths ranging from 10 to 11, indicate the highest uncertainty. These wide intervals suggest that biomass estimates in these areas are more variable or less reliable, possibly due to dynamic or poorly understood environmental conditions.
#In contrast, areas shaded green to blue, with widths between 6 and 9, represent moderate uncertainty. These regions likely experience more stable environmental factors, such as depth, salinity, or nutrient availability, leading to more consistent biomass predictions. The purple areas, with the narrowest confidence interval widths (7.5–8), signify zones of high certainty in biomass estimates. These areas are likely characterized by stable and predictable environmental conditions, such as consistent bathymetry, nutrient flows, or well-documented ecological factors.
#Overall, this map highlights the spatial variability in the reliability of biomass predictions, with certain regions showing greater uncertainty due to complex or less predictable ecological factors


# 2008 --------------------------------------------------------------------


#pca

# We now proceed in the same way for the year 2008
summary(shrimp_data_2008)

# Explore data with multivariate analysis
# Perform PCA with at least two components
pca2 <- dudi.pca(df = shrimp_data_2008[, -c(3:5)], scannf = FALSE, nf = 2)

# Plot the PCA with the first two components
scatter(pca2, xax = 1, yax = 2, clab.row = 0.7, clab.col = 0.7)

# Visualize the loadings
print(pca2$c1)

################################################################################## 
#####################################################SCELTA COVARIATE......

# cs1
shrimp_data_2008$salinity.maxq3 <- factor(shrimp_data_2008$salinity.maxq3) # loading = 0.2564
s.class(pca2$li, fac = shrimp_data_2008$salinity.maxq3)

#cs2
shrimp_data_2008$temp.maxq3 <- factor(shrimp_data_2008$temp.maxq3)
s.class(pca2$li, fac = shrimp_data_2008$temp.maxq3)

shrimp_data_2008$bat <- factor(shrimp_data_2008$bat)
s.class(pca2$li, fac = shrimp_data_2008$bat)

shrimp_data_2008$dist <- factor(shrimp_data_2008$dist)
s.class(pca2$li, fac = shrimp_data_2008$dist)

shrimp_data_2008$slope <- factor(shrimp_data_2008$slope)
s.class(pca2$li, fac = shrimp_data_2008$slope)


# we now proceed again with the log-transformation of total biomass
shrimp_data_2008$log_tot <- log(shrimp_data_2008$tot + 1)

shrimp_geodata_2008_log <- as.geodata(shrimp_data_2008,
                                      coords.col = c("X", "Y"),
                                      data.col = "log_tot",
                                      covar.col = c("bat","salinity.maxq3","temp.maxq3","slope"))

variogram_tot_2008 <- variog(shrimp_geodata_2008_log, trend = "1st", max.dist = 60)

plot(variogram_tot_2008,type="b")

trend_matrix_2008 <- cbind(shrimp_data_2008$bat, 
                           shrimp_data_2008$salinity.minq3, 
                           shrimp_data_2008$temp.maxq3,
                           shrimp_data_2008$slope)

variogram_tot_2008 <- variog(shrimp_geodata_2008_log, trend = "1st",trend.d = trend_matrix_2008, max.dist = 80)

# We now plot the variogram with second order trend
variogram_tot2_2008 <- variog(shrimp_geodata_2008_log, trend = "2nd", max.dist = 60)

plot(variogram_tot2_2008,type="b")

# and we include our chosen covariates in the trend
trend_matrix_2008 <- cbind(shrimp_data_2008$bat, 
                           shrimp_data_2008$salinity.minq3, 
                           shrimp_data_2008$temp.maxq3,
                           shrimp_data_2008$slope)

variogram_tot2_2008 <- variog(shrimp_geodata_2008_log, trend = "2nd",trend.d = trend_matrix_2008, max.dist = 60)

# OVERPARAM...... SCEGLIAMO FIRST ORDER TREND



eyefit(variogram_tot_2008)

#fit_variogram <- variofit(variogram_tot_2008, cov.model = "matern", ini.cov.pars = c(3.78, 6.24), kappa=1.5)
#fit_variogram <- variofit(variogram_tot_2008, cov.model = "exponential", ini.cov.pars = c(6.9, 6.3))
#fit_variogram <- variofit(variogram_tot_2008, cov.model = "spherical", ini.cov.pars = c(4.94, 8.31))
#print(fit_variogram)

fit_matern_lik_2008 <- likfit(shrimp_geodata_2008_log, 
                              cov.model = "matern", 
                              ini.cov.pars = c(4.9, 4.4), 
                              nugget = 3.78, kappa = 1.5)
fit_exponential_lik_2008 <- likfit(shrimp_geodata_2008_log, 
                                   cov.model = "exponential", 
                                   ini.cov.pars = c(6.9, 6.3), 
                                   nugget = 3.78)
fit_spherical_lik_2008 <- likfit(shrimp_geodata_2008_log, 
                                 cov.model = "spherical", 
                                 ini.cov.pars = c(6.4, 8.7), 
                                 nugget = 3.78)

plot(variogram_tot_2008, type = "b", xlab = 'Distance', ylab = 'Semivariance')

lines(fit_matern_lik_2008, col = "red", lwd = 2)
lines(fit_exponential_lik_2008, col = "green", lwd = 2)
lines(fit_spherical_lik_2008, col = "blue", lwd = 2)

vv.mat.2008<-xvalid(shrimp_geodata_2008_log,model=fit_matern_lik_2008)
vv.exp.2008<-xvalid(shrimp_geodata_2008_log,model=fit_exponential_lik_2008)
vv.sph.2008<-xvalid(shrimp_geodata_2008_log,model=fit_spherical_lik_2008)

# Calculate the Mean Squared Error (MSE) for each model
MSE_mat_2008 <- mean(vv.mat.2008$std.error^2)
MSE_exp_2008 <- mean(vv.exp.2008$std.error^2)
MSE_sph_2008 <- mean(vv.sph.2008$std.error^2)

# Calculate the Root Mean Squared Error (RMSE)
RMSE_mat_2008 <- sqrt(MSE_mat_2008)
RMSE_exp_2008 <- sqrt(MSE_exp_2008)
RMSE_sph_2008 <- sqrt(MSE_sph_2008)

# Print the RMSE values
cat("Matern Model: RMSE =", RMSE_mat_2008, "\n")
cat("Exponential Model: RMSE =", RMSE_exp_2008, "\n")
cat("Spherical Model: RMSE =", RMSE_sph_2008, "\n")

#better exponential

#krig

trend.d.2008 <- trend.spatial(~ bat + salinity.minq3 + temp.maxq3 + slope, geodata = shrimp_geodata_2008_log)
trend.l.2008 <- trend.spatial(~ grid_2008$bat + grid_2008$salinity.minq3 + grid_2008$temp.maxq3 + grid_2008$slope)

krige.2008 <- krige.control(
  cov.model = "exponential",
  cov.pars = c(6.9, 6.3), 
  nugget = 3.78,
  trend.d = trend.d.2008,
  trend.l = trend.l.2008
)

krig_2008 <- krige.conv(shrimp_geodata_2008_log, locations = as.matrix(grid_2008[, c("X", "Y")]), krige = krige.2008)


cc_exp_2008 <- data.frame(X = grid_2008$X, Y = grid_2008$Y, Z = krig_2008$predict)

ggplot(cc_exp_2008, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(
    title = "Kriging Interpolation of Total Shrimp Biomass (2008)",
    x = "X Coordinate",
    y = "Y Coordinate",
    fill = "Biomass Prediction") +
  theme_minimal()

#In 2008, the spatial distribution of biomass shows a broader spread of high biomass values. There is an apparent shift or expansion of areas
#with values between 5.0 and 6.5, especially towards the southern region. The gradient from high to low biomass is more diffused compared to 2002,
#suggesting changes in shrimp distribution patterns. This could indicate a response to changing environmental conditions. Shrimp may have moved
#to areas where environmental conditions remained within their optimal range. There is also a high concentration of biomass near the Ligurian
#coast with respect to the 2002 interpolation. The southern coast region benefits from relatively stable oceanographic conditions respect to the
#northen coast, including lower exposure to strong currents or deep-water upwellings. This stability provides a more suitable environment for
#shrimp, enabling them to concentrate in these areas without being dispersed by unfavorable water movements. Now, in order to obtain a
#better and more precise geographical representation, we plot the map of Italy to see how the shrimp biomass is distributed along the areas of the Tirrenean Sea.

# ADD THE MAP OF ITALY

# Transform Italy to UTM (zone 32N, EPSG:32632) for consistency with your coordinates
italy <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf")
italy <- st_transform(italy, crs = 32632)

# Convert kriging results to a grid for ggplot
krig_result_df_08 <- data.frame(
  X = grid_2008$X * 1000, # Assuming the grid units need conversion to match Italy's CRS
  Y = grid_2008$Y * 1000,
  Z = krig_2008$predict
)



# Define the bounds for the area of interest based on your prediction data 
x_min <- min(krig_result_df_08$X) - 10000  # Adjust to give some margin 
x_max <- max(krig_result_df_08$X) + 10000 
y_min <- min(krig_result_df_08$Y) - 10000 
y_max <- max(krig_result_df_08$Y) + 10000 

ggplot() + # Plot the kriging predictions as filled contours 
  geom_raster( data =krig_result_df_08, aes(x = X, y = Y, fill = Z) ) + # Add a color scale for kriging predictions 
  scale_fill_viridis_c( option ="viridis", name = "Biomass Prediction" ) + # Overlay Italy's coastline 
  geom_sf( data = italy, fill = NA, color = "black",lwd = 0.7 ) + # Zoom into the area of interest 
  coord_sf( xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand =FALSE ) + # Add labels and theme 
  labs( title = "Interpolated Shrimp Biomass (2002)", x = "Longitude (km)", y ="Latitude (km)" ) + 
  theme_minimal()

# Assuming shrimp_geodata_2002_log contains the original coordinates and observed log-transformed biomass
observed_points_df08 <- data.frame(
  X = shrimp_geodata_2008_log$coords[, 1] * 1000,  # Convert to match UTM coordinates
  Y = shrimp_geodata_2008_log$coords[, 2] * 1000,  # Convert to match UTM coordinates
  Biomass = shrimp_geodata_2002_log$data           # Log-transformed biomass
)

# Plot with Italy map, kriging results, and observed points
ggplot() +
  # Kriging predictions as a raster layer
  geom_raster(data = krig_result_df_08, aes(x = X, y = Y, fill = Z)) +
  # Color scale for kriging predictions
  scale_fill_viridis_c(option = "viridis", name = "Biomass Prediction") +
  # Overlay Italy's coastline
  geom_sf(data = italy, fill = NA, color = "black", lwd = 0.7) +
  # Add observed biomass points, with size proportional to the biomass
  geom_point(data = observed_points_df08, aes(x = X, y = Y, size = Biomass), color = "red", alpha = 0.6) +
  # Size scale for observed points
  scale_size_continuous(name = "Observed Biomass", range = c(1, 5)) +  # Adjust size range as needed
  # Zoom into the area of interest
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE) +
  # Labels and theme
  labs(title = "Interpolated Shrimp Biomass (2002)", x = "Longitude (km)", y = "Latitude (km)") +
  theme_minimal()

#The figure provides a comprehensive visualization of the spatial distribution of shrimp biomass for the year 2008, derived through
#kriging interpolation. There seems to be a good alignment between the observed biomass values (red points) and the predicted values
#(background gradient). The kriging framework has successfully interpolated shrimp biomass across the region, capturing spatial trends
#influenced by environmental factors. The visualization suggests that the model is reasonably predictive, but further validation (e.g., cross-validation) is necessary to confirm its accuracy. 


#confidence intervals

low_2008 <- krig_2008$predict- 1.96*sqrt(krig_2008$krige.var)
up_2008 <- krig_2008$predict + 1.96*sqrt(krig_2008$krige.var)

hist(low_2008)
hist(up_2008)
hist(up_2008-low_2008)

#The distribution of the lower bound (low) is skewed towards the lower end, with most values concentred between -8 and-4. This indicates that
#the lower bound of the biomass predictions is generally on the lower side, reflecting uncertainty and variability in areas with potentially
#low biomass. Then the upper confidence interval (up) values are more symmetric and centered around values between 9 and 11, with fewer
#extreme values compared to the lower bound. This shows that the upper limits of biomass predictions are more consistent and tend to reflect a higher level of potential shrimp biomass.
#The histogram of the difference (up - low) between the upper and lower bounds is unimodal and primarily concentrated between 12 and 18, indicating a consistent width of the confidence intervals across most of
#the grid. The difference histogram suggests that uncertainty in biomass predictions is fairly uniform across the study area, which can be a positive sign if the model is well-calibrated.

cc_lower_2008 <- data.frame(X = grid_2008$X, Y = grid_2008$Y, Z = low_2008)
cc_upper_2008 <- data.frame(X = grid_2008$X, Y = grid_2008$Y, Z = up_2008)

ggplot(cc_lower_2008, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(title = "Lower 95% Confidence Bound 2008") +
  theme_minimal()

#This map depicts the spatial variation in the lower 95% confidence bound for total biomass along the Gaeta-Genova coast in 2008. The lighter
#regions (green to yellow), with values ranging from approximately -1 to 0 (log scale), represent areas with higher biomass even under
#conservative estimates. These zones are likely associated with more productive or nutrient-rich habitats, such as shallow coastal zones or
#estuarine areas, where conditions like moderate salinity and sunlight penetration promote biomass growth. In contrast, darker regions (purple to blue), with values from approximately -9 to -4, indicate areas with 
#lower conservative biomass estimates. These regions may correspond to deeper waters or zones with less favorable salinity and nutrient
#conditions, which can limit the accumulation of biomass. Depth is a key factor influencing biomass distribution, as shallow coastal areas
#typically experience greater sunlight penetration, enabling photosynthesis and supporting a richer food web for marine life. The
#lighter zones on this map likely indicate biologically productive areas where biomass remains relatively high even when accounting for
#variability. These areas may benefit from moderate salinity levels influenced by freshwater inflows, creating environments that are
#conducive to sustaining higher biomass levels. Overall, the map highlights productive regions (in green to yellow) along the
#Gaeta-Genova coastline, likely reflecting the interplay of favorable bathymetry, salinity, and nutrient availability.

ggplot(cc_upper_2008, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(title = "Upper 95% Confidence Bound 2008") +
  theme_minimal()

#This map showcases the upper 95% confidence bound for shrimp density (log scale) along the Gaeta-Genova coast in 2008, offering an optimistic
#projection of potential biomass distribution. The regions shaded in green to yellow, with values ranging from 11 to 13.5 on the log scale,
#represent areas with the highest potential biomass. These zones are likely associated with optimal environmental conditions, such as
#favorable bathymetry (moderate depths), adequate salinity levels, and nutrient-rich waters, which together promote higher biomass
#productivity. Conversely, darker regions (shades of purple to blue), with values between 6.5 and 8.0, indicate areas where the potential
#biomass is lower, even under favorable conditions. These areas may be linked to less suitable environmental factors, such as greater depths,
#reduced sunlight penetration, or lower nutrient availability, which limit biomass accumulation. Bathymetry plays a central role in biomass
#distribution, as shallower areas are generally more productive due to greater sunlight availability that supports photosynthesis, fostering
#robust marine food webs. Additionally, moderate salinity, often influenced by freshwater inflows in estuarine or near-coastal areas,
#enhances nutrient availability, creating conditions conducive to higher biomass. This map reflects these influences, highlighting areas of
#potential biological richness (green to yellow zones) and offering insights into how depth and salinity gradients shape biomass distribution along the coastline.

interval_width_2008 <- up_2008- low_2008
cc_interval <- data.frame(X = grid_2008$X, Y = grid_2008$Y, Z = interval_width_2008)
ggplot(cc_interval, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(title = "Width of 95% Confidence Interval (2008)") +
  theme_minimal()

# This map illustrates the width of the 95% confidence interval for biomass estimates along the coastal stretch from Gaeta to Liguria in the
#northwestern Mediterranean in 2008. The width of the confidence interval serves as an indicator of the uncertainty associated with biomass
#predictions: High Uncertainty (Yellow, 17.5–18.0): Predominantly observed along the outer edges, these regions exhibit the widest #confidence intervals, reflecting significant variability in biomass
#estimates. Such uncertainty may arise from dynamic environmental factors, including fluctuating currents, variable nutrient availability,
#or human activities like fishing that impact local biomass levels. Moderate Uncertainty (Green, 13.0–16.5): Found in intermediate regions,
#these areas represent moderate confidence in biomass estimates. The variability may result from transitional oceanographic conditions,
#proximity to mixed ecological zones, or the influence of seasonal factors.Low Uncertainty (Purple to Blue, 10.5–12.5): Concentrated
#centrally along the coast, these regions have the narrowest confidence intervals, indicating stable and predictable conditions that enhance the
#precision of biomass estimates. These zones are likely characterized by consistent environmental factors, such as optimal depth, stable nutrient
#distribution, and salinity levels, which create favorable and predictable habitats for marine life.

