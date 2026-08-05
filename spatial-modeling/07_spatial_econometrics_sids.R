# Spatial Statistics and Statistical Tools for Environmental Data
# MSc in Statistical Methods and Applications, Sapienza University of Rome
# Group 12 — Agate, Bartocci, Natali
#
# Homework 7 — Areal spatial analysis (SIDS, North Carolina)
# Spatial weights, spatial autocorrelation and spatial regression on areal data with spdep/spatialreg.

# We first install and load the necessary libraries
# install.packages(c("spData", "sf", "spdep", "spatialreg", "tigris", "ggplot2"))
options(tigris_use_cache = TRUE)

library(spData)      # For spatial datasets
library(sf)          # For spatial objects and handling shapefiles
library(spdep)       # For spatial econometrics functions
library(tigris)      # For accessing shapefiles of US counties
library(ggplot2)     # For data visualization

# Load the Sudden Infant Death Syndrome (SIDS) data for North Carolina
data("nc.sids", package = "spData")
nc.sids <- as.data.frame(nc.sids)  # Convert to data frame for easier manipulation

# Load geographic data (North Carolina counties) from the "tigris" package
nc_geom <- counties(state = "NC", cb = TRUE, class = "sf") 

#After loading the dataset and the required geographic data,
# we can now plot the geographic data (county boundaries)
plot(st_geometry(nc_geom))

# Merge the SIDS data with the geographic data by county names
nc.sids$NAME <- rownames(nc.sids)  # Add county names to SIDS data
nc_merged <- merge(nc.sids, nc_geom, by = "NAME")  # Merge datasets by "NAME"
nc_merged <- st_as_sf(nc_merged)  # Convert to spatial object

# We can now visualize the spatial distribution of SIDS in 1974
ggplot(data = nc_merged) +
  geom_sf(aes(fill = SID74), color = "white") +  # Fill counties based on SIDS rate in 1974
  scale_fill_gradient(low = "green", high = "red", name = "SID74") +
  theme_minimal() +
  labs(
    title = "Distribution of Sudden Infant Deaths in 1974 (North Carolina)",
    subtitle = "SIDS by County",
    fill = "SIDS"
  )

# We plot the highest SIDS rates in 1974 and 1979
high_state_1974 <- nc_merged[which.max(nc_merged$SID74), c("NAME", "SID74")]
high_state_1979 <- nc_merged[which.max(nc_merged$SID79), c("NAME", "SID79")]
cat("Highest SIDS in 1974:", high_state_1974$NAME, "with", high_state_1974$SID74, "cases\n") 
#Mecklenburg with 44 cases
cat("Highest SIDS in 1979:", high_state_1979$NAME, "with", high_state_1979$SID79, "cases\n") 
#Cumberland with 57 cases

# Visualize the spatial distribution of SIDS in 1979
ggplot(data = nc_merged) +
  geom_sf(aes(fill = SID79), color = "white") +  # Fill counties based on SIDS rate in 1979
  scale_fill_gradient(low = "green", high = "red", name = "SID79") +
  theme_minimal() +
  labs(
    title = "Distribution of Sudden Infant Deaths in 1979 (North Carolina)",
    subtitle = "SIDS by County",
    fill = "SIDS"
  )

# Calculate total SIDS deaths in 1974 and 1979
total_sid74 <- sum(nc_merged$SID74, na.rm = TRUE)
total_sid79 <- sum(nc_merged$SID79, na.rm = TRUE)
cat("Total SID74:", total_sid74, "\n") #667
cat("Total SID79:", total_sid79, "\n") #836

# Explore correlation between SIDS and birth rates in the two years
correlation_1974 <- cor(nc_merged$SID74, nc_merged$NWBIR74, use = "complete.obs")
correlation_1979 <- cor(nc_merged$SID79, nc_merged$NWBIR79, use = "complete.obs")
cat("Correlation between SIDS and NWBIR (1974):", correlation_1974, "\n") #0.8977063
cat("Correlation between SIDS and NWBIR (1979):", correlation_1979, "\n") #0.8586045



# 1) Compare pseudo-likelihood estimates with those obtained using spautolm (library "spatialreg") --------

# At this point, we can now create a spatial neighborhood structure using contiguity-based neighbors
nc_nb <- poly2nb(nc_merged)  # Create neighbor list based on shared borders
summary(nc_nb)

# and convert the neighborhood list to a spatial weights matrix (row-standardized)
listw <- nb2listw(nc_nb, style = "W")

# We now load the libraries "lagsarlmtree" and "spatialreg" for spatial autoregressive models:
# to fit such autoregressive models we will utilize respectively the functions
# "lagsarlm" and "spautolm".
# The lagsarlm model uses pseudolikelihood estimation, which is an approximation to the full likelihood using 
# the product of conditional distributions, often used to simplify computations in spatial autoregressive models.
# The spautolm model instead uses full maximum likelihood estimation, which is theoretically more accurate
# but computationally intensive.
library(lagsarlmtree)
library(spatialreg) 


# Now we fit spatial autoregressive models (SAR) for 1974 and 1979 SIDS data
model_lagsarlm_74 <- lagsarlm(SID74 ~ NWBIR74, data = nc_merged, listw = listw)
model_spautolm_74 <- spautolm(SID74 ~ NWBIR74, data = nc_merged, listw = listw, method = "eigen")
model_lagsarlm_79 <- lagsarlm(SID79 ~ NWBIR79, data = nc_merged, listw = listw)
model_spautolm_79 <- spautolm(SID79 ~ NWBIR79, data = nc_merged, listw = listw, method = "eigen")

# and we print summaries of the models to compare them
summary(model_lagsarlm_74)
summary(model_spautolm_74)
summary(model_lagsarlm_79)
summary(model_spautolm_79)

# We now compare also the Akaike Information Criterion (AIC) and Bayesian Information Criterion (BIC) 
# for the spatial autoregressive models to evaluate their relative performance in terms of goodness-of-fit 
# and model complexity. 
data_74 <- data.frame(
  Model = c("lagsarlm", "spautolm"),
  AIC = c(AIC(model_lagsarlm_74), AIC(model_spautolm_74)),
  BIC = c(BIC(model_lagsarlm_74), BIC(model_spautolm_74))
)
data_79 <- data.frame(
  Model = c("lagsarlm", "spautolm"),
  AIC = c(AIC(model_lagsarlm_79), AIC(model_spautolm_79)),
  BIC = c(BIC(model_lagsarlm_79), BIC(model_spautolm_79))
)
print(data_74)
print(data_79)

# print(data_74)
# Model      AIC      BIC
# 1 lagsarlm 537.1068 547.5275
# 2 spautolm 537.1092 547.5299

# print(data_79)
# Model      AIC      BIC
# 1 lagsarlm 604.0076 614.4282
# 2 spautolm 604.6919 615.1126


# We notice that The AIC and BIC values for lagsarlm are slightly lower than those for spautolm in both years.
# This indicates that lagsarlm provides a marginally better fit to the data while penalizing model complexity.
# We also know that if the performance (AIC/BIC) of the pseudolikelihood method (lagsarlm) 
# is similar or better than that of the maximum likelihood method (spautolm), 
# it suggests that the approximation does not significantly compromise model quality.
# Hence, in conclusion for both datasets, lagsarlm offers a marginally better balance 
# between goodness-of-fit and complexity.



# Finally, we add additional analysis to check residuals and spatial autocorrelation

# Calculate residuals for the 1974 and 1979 models
residuals_lag_74 <- residuals(model_lagsarlm_74)
residuals_spa_74 <- residuals(model_spautolm_74)
residuals_lag_79 <- residuals(model_lagsarlm_79)
residuals_spa_79 <- residuals(model_spautolm_79)

# Plot residuals to visually inspect model fit (year 1974)
ggplot(nc_merged) +
  geom_sf(aes(fill = residuals_lag_74), color = "white") +
  scale_fill_gradient(low = "green", high = "darkred", name = "Residuals 1974") +
  theme_minimal() +
  labs(title = "Residuals of the SAR Model (1974)")

ggplot(nc_merged) +
  geom_sf(aes(fill = residuals_spa_74), color = "white") +
  scale_fill_gradient(low = "green", high = "darkred", name = "Residuals 1974 SPA") +
  theme_minimal() +
  labs(title = "Residuals of the SpatAutolm Model (1974)")


# To test whether the spatial autocorrelation of residuals is random or exhibits a pattern 
# (e.g., clustering or dispersion), we apply Moran's I test.
# High positive autocorrelation implies clustering of similar values, 
# while negative autocorrelation indicates dispersion or a checkerboard-like pattern.
moran_test_74_lag <- moran.test(residuals_lag_74, listw)
moran_test_74_spa <- moran.test(residuals_spa_74, listw)
moran_test_79_lag <- moran.test(residuals_lag_79, listw)
moran_test_79_spa <- moran.test(residuals_spa_79, listw)

# Print Moran's I test results for residuals
cat("Moran's I Test for 1974 Lag Model Residuals:\n")
print(moran_test_74_lag)
cat("Moran's I Test for 1974 SpatAutolm Model Residuals:\n")
print(moran_test_74_spa)
cat("Moran's I Test for 1979 Lag Model Residuals:\n")
print(moran_test_79_lag)
cat("Moran's I Test for 1979 SpatAutolm Model Residuals:\n")
print(moran_test_79_spa)





# 2) Choose the "best" neighborhood structure for a CAR model on these data and justify your choice --------

# Ensure complete cases for SID79 and NWBIR79
nc_merged <- nc_merged[complete.cases(nc.sids[, c("SID79", "NWBIR79")]), ]

# Now we create a spatial neighborhood structure not contiguity-based,
# but using distance-based neighbors (within 50km) instead
coords <- st_coordinates(st_centroid(nc_merged))
dist_nb <- dnearneigh(coords, 0, 50)  # 50 km distance threshold
listw_dist <- nb2listw(dist_nb)

# Fit spatial autoregressive models for 1974 and 1979 SIDS data using this distance-based neighborhood
model_lagsarlm_dist_74 <- lagsarlm(SID74 ~ NWBIR74, data = nc_merged, listw = listw_dist)
summary(model_lagsarlm_dist_74)
# AIC: 448.79 (AIC for lm: 535.2)

model_lagsarlm_dist_79 <- lagsarlm(SID79 ~ NWBIR79, data = nc_merged, listw = listw_dist)
summary(model_lagsarlm_dist_79)
# AIC: 567.55 (AIC for lm: 603.96)

# We now compare AIC and BIC for both contiguity-based (polygonal) and distance-based (Euclidean) 
# neighborhoods for 1974
data_74 <- data.frame(
  Neighborhood_Structure = c("Polygonal", "Euclidean"),
  AIC = c(AIC(model_lagsarlm_74), AIC(model_lagsarlm_dist_74)),
  BIC = c(BIC(model_lagsarlm_74), BIC(model_lagsarlm_dist_74))
)
names(data_74) = c("Neighborhood Structure", "AIC", "BIC")

# Compare AIC and BIC for both contiguity-based (polygonal) and distance-based (Euclidean) neighborhoods for 1979
data_79 <- data.frame(
  Neighborhood_Structure = c("Polygonal", "Euclidean"),
  AIC = c(AIC(model_lagsarlm_79), AIC(model_lagsarlm_dist_79)),
  BIC = c(BIC(model_lagsarlm_79), BIC(model_lagsarlm_dist_79))
)
names(data_79) = c("Neighborhood Structure", "AIC", "BIC")

# Print the AIC and BIC comparison tables for both years
print(data_74)
print(data_79)

# print(data_74)
# Neighborhood Structure      AIC      BIC
# 1              Polygonal 537.1068 547.5275
# 2              Euclidean 448.7907 459.2114

# print(data_79)
# Neighborhood Structure      AIC      BIC
# 1              Polygonal 604.0076 614.4282
# 2              Euclidean 567.5479 577.9686


# In both years, we can clearly notice how the distance-based neighborhood structure has 
# substantially lower AIC and BIC values, thus suggesting that a Euclidean-based 
# spatial weighting better captures the spatial relationships the dataset.





# Finally, we can compare the two neighborhood structures also in several other ways:

# A) Visualizing the neighborhood structures

# Contiguity (Polygonal neighborhood)
plot(st_geometry(nc_geom), main = "Polygonal Neighborhoods (Contiguity)")
plot(st_geometry(nc_merged), add = TRUE, col = "transparent", border = "red")
plot(st_geometry(nc_geom[poly2nb(nc_merged) != NULL,]), add = TRUE, col = "blue", lwd = 2) # Highlight neighbors

# Distance-based neighborhood (Euclidean) plot
# Plot the distance-based neighbors (50 km radius) on top of the map
plot(st_geometry(nc_geom), main = "Euclidean Neighborhoods (50 km Distance)")
plot(st_geometry(nc_merged), add = TRUE, col = "transparent", border = "red")
for (i in 1:length(dist_nb)) {
  plot(st_geometry(nc_merged[i,]), add = TRUE, col = ifelse(i %% 2 == 1, "green", "blue"))
}

# B) Comparison of residuals: Polygonal vs Euclidean Models
# Calculate residuals for both models (lag and distance)

residuals_74 <- data.frame(
  NAME = nc_merged$NAME,
  LagSAR = residuals(model_lagsarlm_74),
  SpAutolm = residuals(model_spautolm_74)
)

residuals_79 <- data.frame(
  NAME = nc_merged$NAME,
  LagSAR = residuals(model_lagsarlm_79),
  SpAutolm = residuals(model_spautolm_79)
)

# Add residuals to spatial data for visualization
nc_merged$res_lag_74 <- residuals_74$LagSAR
nc_merged$res_spa_74 <- residuals_74$SpAutolm
nc_merged$res_lag_79 <- residuals_79$LagSAR
nc_merged$res_spa_79 <- residuals_79$SpAutolm

# Plot residuals for 1974 and 1979 (SAR and SpAutolm models)

ggplot(nc_merged) +
  geom_sf(aes(fill = res_lag_74), color = "white") +
  scale_fill_gradient(low = "green", high = "darkred", name = "Residuals (LagSAR)") +
  theme_minimal() +
  labs(title = "Residuals of LagSAR Model (1974)")

ggplot(nc_merged) +
  geom_sf(aes(fill = res_spa_74), color = "white") +
  scale_fill_gradient(low = "green", high = "darkred", name = "Residuals (SpAutolm)") +
  theme_minimal() +
  labs(title = "Residuals of SpAutolm Model (1974)")

# Plotting histograms for the residuals of two models (Lag Model and Spatial Model) for the year 1974
ggplot(nc_merged) +
  geom_histogram(aes(x = res_lag_74, fill = "Lag Model"), 
                 alpha = 0.6, 
                 position = "dodge", 
                 bins = 30, 
                 color = "black", 
                 linewidth = 0.5) +  
  geom_histogram(aes(x = res_spa_74, fill = "Spatial Model"), 
                 alpha = 0.6, 
                 position = "dodge", 
                 bins = 30, 
                 color = "black", 
                 linewidth = 0.5) +  
  scale_fill_manual(values = c("Lag Model" = "#1f77b4", "Spatial Model" = "#ff7f0e")) +  
  labs(title = "Residual Comparison for 1974 Models",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.title = element_blank(), 
        legend.position = "top",
        legend.text = element_text(size = 12),  
        axis.title = element_text(size = 14),    
        axis.text = element_text(size = 12))  

#The residual analysis for the year 1974 shows that both models, Lag and Spatial, 
#exhibit a good concentration of residuals around zero, showing a significant overlap 
#between the residual distributions of the Lag model and the Spatial model, 
#indicating that both models explain a substantial portion of the data's variability and produce 
#a large number of predictions with small errors.


ggplot(nc_merged) +
  geom_sf(aes(fill = res_lag_79), color = "white") +
  scale_fill_gradient(low = "green", high = "darkred", name = "Residuals (LagSAR)") +
  theme_minimal() +
  labs(title = "Residuals of LagSAR Model (1979)")

ggplot(nc_merged) +
  geom_sf(aes(fill = res_spa_79), color = "white") +
  scale_fill_gradient(low = "green", high = "darkred", name = "Residuals (SpAutolm)") +
  theme_minimal() +
  labs(title = "Residuals of SpAutolm Model (1979)")

# Plotting histograms for the residuals of two models (Lag Model and Spatial Model) for the year 1979
ggplot(nc_merged) +
  geom_histogram(aes(x = res_lag_79, fill = "Lag Model"), 
                 alpha = 0.6, 
                 position = "dodge", 
                 bins = 30, 
                 color = "black", 
                 linewidth = 0.5) +  
  geom_histogram(aes(x = res_spa_79, fill = "Spatial Model"), 
                 alpha = 0.6, 
                 position = "dodge", 
                 bins = 30, 
                 color = "black", 
                 linewidth = 0.5) +  
  scale_fill_manual(values = c("Lag Model" = "#1f77b4", "Spatial Model" = "#ff7f0e")) +  
  labs(title = "Residual Comparison for 1979 Models",
       x = "Residuals",
       y = "Frequency") +
  theme_minimal() +
  theme(legend.title = element_blank(), 
        legend.position = "top",
        legend.text = element_text(size = 12),  
        axis.title = element_text(size = 14),    
        axis.text = element_text(size = 12))  

#While in 1974 both models showed a relatively similar performance, with a considerable overlap of residuals, 
#in 1979 the Spatial model appears to offer a significantly better fit, 
#with a more concentrated residual distribution and fewer extreme values compared to the Lag model.

# C) Moran's I test for residuals comparison

# As in Part 1, we now apply Moran's I test for spatial autocorrelation in residuals comparison

# For Polygonal Neighborhood (Lag Model)
moran_test_lag_74 <- moran.test(nc_merged$res_lag_74, listw)
moran_test_lag_79 <- moran.test(nc_merged$res_lag_79, listw)

cat("Moran's I for Residuals from SAR Model with Polygonal Neighborhood (1974):\n")
print(moran_test_lag_74)

cat("Moran's I for Residuals from SAR Model with Polygonal Neighborhood (1979):\n")
print(moran_test_lag_79)

# For Euclidean Neighborhood (Distance Model)
moran_test_dist_74 <- moran.test(nc_merged$res_spa_74, listw_dist)
moran_test_dist_79 <- moran.test(nc_merged$res_spa_79, listw_dist)

cat("Moran's I for Residuals from SAR Model with Euclidean Neighborhood (1974):\n")
print(moran_test_dist_74)

cat("Moran's I for Residuals from SAR Model with Euclidean Neighborhood (1979):\n")
print(moran_test_dist_79)

#In summary, the overall results show that the residuals from the SAR models, 
#both with polygonal and Euclidean proximity, do not exhibit significant spatial autocorrelation. 
#This suggests that the models may have been well-specified from a spatial perspective.
