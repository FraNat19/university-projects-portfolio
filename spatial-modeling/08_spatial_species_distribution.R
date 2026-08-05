# Spatial Statistics and Statistical Tools for Environmental Data
# MSc in Statistical Methods and Applications, Sapienza University of Rome
# Group 12 — Agate, Bartocci, Natali
#
# Homework 8 — Spatial species distribution modelling (mite data)
Presence-absence modelling of species with spatial predictors; evaluation with ROC/AUC.

# We first install and import the required libraries
library(vegan)
library(spdep)
library(sp)
library(spatialreg)
library(ggplot2)
library(sf)
library(pROC)
library(MLmetrics)
library(reshape2)
library(dplyr)
library(ROCR)

# Load datasets
data("mite")      # Species abundance data
data("mite.env")  # Environmental variables
data("mite.xy")   # Spatial coordinates

# 1. Explore the Environmental Variables -----------------------------------------
# Using the data in the mite datasets (from library vegan in R) create the presence-absence variable for species ONOV and LCIL

summary(mite.env)

# Summary of Environmental Variables in mite.env:

# 1. SubsDens (Substrate Density)
# Range: 21.17–80.59; Mean: 39.28; Median: 36.38; Q1: 30.01; Q3: 46.81.
# Higher density might influence mite populations by affecting habitat structure.

# 2. WatrCont (Water Content)
# Range: 134.1–827.0; Mean: 410.6; Median: 398.5; Quartiles as above.
# Moisture levels likely impact mite habitats, particularly in sensitive areas.

# 3. Substrate (Categorical)
# Types: Sphagn1 (25), Sphagn2 (11), Sphagn3 (1), Sphagn4 (2), Litter (2), Barepeat (2), Interface (27).
# Sphagn1 and Interface are the most frequent, reflecting prevalent substrate types.

# 4. Shrub (Categorical)
# Levels: None (19), Few (26), Many (25). Indicates shrub density at sites.

# 5. Topo (Topography)
# Levels: Blanket (44), Hummock (26). Represents terrain features (flat vs. uneven).

# The summary highlights key patterns in environmental variables, offering insight into potential factors affecting mite distribution.

# Now we see frequency of levels for categorical variables
lapply(mite.env[, sapply(mite.env, is.factor)], table)

# Visualize distributions
env_long <- melt(mite.env)  # Convert to long format for ggplot

ggplot(env_long, aes(x = value)) +
  geom_histogram(bins = 20, fill = "blue", color = "black") +
  facet_wrap(~variable, scales = "free") +
  theme_minimal() +
  labs(title = "Distributions of Environmental Variables", x = "Value", y = "Count")

#The histograms illustrate the distributions of the environmental variables. 
#The SubsDens variable is right-skewed, with most values concentrated between 20 and 50, 
#indicating a generally low density in most observations. Conversely, 
#the WatrCont variable shows a more symmetric, bell-shaped distribution, 
#with values predominantly around 400-500, suggesting more uniform variability. 
#These patterns highlight potential differences in variability and may require transformation or specific 
#modeling strategies to account for their characteristics.

# 1. Explore the available variables ---------------------------------------
# Now we create Presence-Absence Variables
mite_pa <- mite
mite_pa$ONOV <- ifelse(mite$ONOV > 0, 1, 0)
mite_pa$LCIL <- ifelse(mite$LCIL > 0, 1, 0)

ONOV <- mite_pa[, "ONOV"]  # Presence-absence for ONOV
LCIL <- mite_pa[, "LCIL"]  # Presence-absence for LCIL

# Check the first few rows of the modified dataset
head(mite_pa[, c("ONOV", "LCIL")])

# Summarize the presence-absence for the entire dataset
table(ONOV) # Counts for ONOV
#0 1 
#7 63
table(LCIL) # Counts for LCIL
#0 1 
# 15 55

# Alternatively, create a summary for both species together
colSums(mite_pa[, c("ONOV", "LCIL")])  # Total presence counts
colMeans(mite_pa[, c("ONOV", "LCIL")]) # Proportion of presence


# We now create a new dataframe for LCIL and ONOV with coordinates, presence status, water content, and substrate density
onov <- data.frame(x = mite.xy$x, y = mite.xy$y, water = mite.env$WatrCont,
                   presence = ONOV, density = mite.env$SubsDens)
lcil <- data.frame(x = mite.xy$x, y = mite.xy$y, water = mite.env$WatrCont,
                   presence = LCIL, density = mite.env$SubsDens)

# If we want to plot for instance the correlation between water content and mite presence we obtain
cor(lcil$water, lcil$presence)
# 0.5080447
cor(onov$water, onov$presence)
# -0.5486648


# 2. Build an Autologistic Model for each Species, choosing the best Neighbourhood Structure ----------------------------------------------

# The first crucial point is to build a suitable neighborhood structure: we try a first-order
# neighborhood at first
coordinates(lcil) <- ~x + y # Converting data for spatial analysis
lcilkn <- knearneigh(coordinates(lcil), k = 4)
# Finding nearest neighbors with k=4
plot(lcilkn$x, main = "First level neighborhood")
plot(knn2nb(lcilkn), lcilkn$x, add = TRUE)


# Performing neighborhood analysis for ONOV
coordinates(onov) <- ~x + y
onovkn <- knearneigh(coordinates(onov), k = 4)
plot(onovkn$x, main = "First level neighborhood - ONOV")
plot(knn2nb(onovkn), onovkn$x, add = TRUE)
# We obtain the same plot and, as expected, such neighborhood does not suit the model well

# Therefore, another feasible possibility is trying to coinstruct a distance-based neighborhood,
# based on Euclidean distance, through the "dnearneigh" function
lcil.dn <- dnearneigh(coordinates(lcil), 0, 0.8)
plot(lcilkn$x, main = "Distance-based neighborhood")
plot(lcil.dn, lcilkn$x, add = TRUE)

# Using distance-based neighborhood for ONOV
onov.dn <- dnearneigh(coordinates(onov), 0, 0.8)
plot(onovkn$x, main = "Distance-based neighborhood")
plot(onov.dn, onovkn$x, add = TRUE)


# Now we plot the summary of the obtained distance-based neighborhood structure for both species
#Summary for ONOV
summary(onov.dn)

#Neighbour list object:
 # Number of regions: 70 
#Number of nonzero links: 280 
#Percentage nonzero weights: 5.714286 
#Average number of links: 4 
#3 regions with no links:
 # 7, 10, 70
#5 disjoint connected subgraphs
#Link number distribution:
  
 # 0  1  2  3  4  5  6  7  8  9 
#3  2  9 15 20  5  8  3  4  1 
#2 least connected regions:
 # 1 20 with 1 link
#1 most connected region:
 # 39 with 9 links

summary(lcil.dn)

#Neighbour list object:
 # Number of regions: 70 
#Number of nonzero links: 280 
#Percentage nonzero weights: 5.714286 
#Average number of links: 4 
#3 regions with no links:
 # 7, 10, 70
#5 disjoint connected subgraphs
#Link number distribution:
  
 # 0  1  2  3  4  5  6  7  8  9 
#3  2  9 15 20  5  8  3  4  1 
#2 least connected regions:
 # 1 20 with 1 link
#1 most connected region:
 # 39 with 9 links


# We notice that three regions (i.e. 7, 10, 70) have no links:
# we therefore remove those three isolated points to avoid potential issues
#We first remove the isolated regions for ONOV
onov <- onov[-c(7, 10, 70),]
onov.dn <- dnearneigh(coordinates(onov), 0, 0.8)
# and for LCIL model
lcil <- lcil[-c(7, 10, 70),]
lcil.dn <- dnearneigh(coordinates(lcil), 0, 0.8)


# Now we can create an adjacency matrix based on the neighborhood structure
adj_lcil <- nb2mat(lcil.dn, style = "B", zero.policy = TRUE)
n_lcil <- nrow(lcil)
n1i_lcil <- rep(0, n_lcil)
numnbr_lcil <- rep(0, n_lcil)
for (i in 1:n_lcil) {
  n1i_lcil[i] <- sum(adj_lcil[,i] * lcil$presence[])
  numnbr_lcil[i] <- sum(adj_lcil[,i])
}
ni_lcil <- 4 * n1i_lcil / numnbr_lcil

# Create adjacency matrix for ONOV
adj_onov <- nb2mat(onov.dn, style = "B", zero.policy = TRUE)
n_onov <- nrow(onov)
n1i_onov <- rep(0, n_onov)
numnbr_onov <- rep(0, n_onov)
for (i in 1:n_onov) {
  n1i_onov[i] <- sum(adj_onov[,i] * onov$presence[])
  numnbr_onov[i] <- sum(adj_onov[,i])
}
ni_onov <- 4 * n1i_onov / numnbr_onov



# 2. Autologistic Models -----------------------------------------------------------------------------

# Since we are working with presence/absence data, we are unable to fit models using functions
# based on the Gaussian assumption (such as the "lagsarlm" function).
# Therefore, a feasible way to estimate the autologistic models is to utilize the glm function 
# (with a binomial family, which is typically used for binary response data).
# Now we will therefore build four logistic models for each species respectively:

#For ONOV
# Simple logistic model with only the neighborhood effects
onov1 <- glm(presence ~ ni_onov, data = onov, family = binomial)
pred1_onov <- predict(onov1, type = "response")
fit1_onov <- ifelse(pred1_onov > 0.5, 1, 0)

# A logistic model adding the water content effect
onov2 <- glm(presence ~ ni_onov + water, data = onov, family = binomial)
pred2_onov <- predict(onov2, type = "response")
fit2_onov <- ifelse(pred2_onov > 0.5, 1, 0)

# Logistic model with added substrate density effect
onov3 <- glm(presence ~ ni_onov + density, data = onov, family = binomial)
pred3_onov <- predict(onov3, type = "response")
fit3_onov <- ifelse(pred3_onov > 0.5, 1, 0)

# Logistic model with the effects of both covariates (water content and substrate density)
onov4 <- glm(presence ~ ni_onov + water + density, data = onov, family = binomial)
pred4_onov <- predict(onov4, type = "response")
fit4_onov <- ifelse(pred4_onov > 0.5, 1, 0)

# We proceed in the same way for LCIL
lcil1 <- glm(presence ~ ni_lcil, data = lcil, family = binomial)
pred1_lcil <- predict(lcil1, type = "response")
fit1_lcil <- ifelse(pred1_lcil > 0.5, 1, 0)


lcil2 <- glm(presence ~ ni_lcil + water, data = lcil, family = binomial)
pred2_lcil <- predict(lcil2, type = "response")
fit2_lcil <- ifelse(pred2_lcil > 0.5, 1, 0)


lcil3 <- glm(presence ~ ni_lcil + density, data = lcil, family = binomial)
pred3_lcil <- predict(lcil3, type = "response")
fit3_lcil <- ifelse(pred3_lcil > 0.5, 1, 0)


lcil4 <- glm(presence ~ ni_lcil + water + density,
             data = lcil, family = binomial)
pred4_lcil <- predict(lcil4, type = "response")
fit4_lcil <- ifelse(pred4_lcil > 0.5, 1, 0)


# We can also compare the models through the AIC
AIC(onov1, onov2, onov3, onov4, lcil1, lcil2, lcil3, lcil4)

#       df      AIC
#onov1  2 46.41517
#onov2  3 27.29905
#onov3  3 39.52509
#onov4  4 28.90828
#lcil1  2 58.72316
#lcil2  3 49.17438
#lcil3  3 60.48429
#lcil4  4 50.83430




# 3. Draw a predictive Map for all models and for the simple logistic model --------------------------------------------------------

# Firs choose good colors with viridis
colors_lcil <- viridis(length(unique(lcil$water)))
colors_onov <- viridis(length(unique(onov$water)))

# Function to create plots
create_plot <- function(data, colors, fit, main_title) {
  # Convert coordinates to degrees with sf
  data_sf <- st_as_sf(data, coords = c("x", "y"), crs = 4326)
  data_points <- data.frame(
    st_coordinates(data_sf),  # Extract coordinates in degrees
    density = data$density,
    fit = fit
  )
  
  # Create the plot with ggplot
  ggplot(data_points, aes(x = X, y = Y)) +
    geom_point(aes(color = density, size = density), alpha = 0.8) +
    scale_color_viridis_c(option = "viridis") +
    scale_size_continuous(range = c(3, 8), guide = "none") +
    theme_minimal() +
    labs(
      title = main_title,
      x = "Longitude (°)",
      y = "Latitude (°)"
    ) +
    coord_cartesian(xlim = c(0, 2.5), ylim = c(0, 10))  # Custom limits
}

# Plot for ONOV with the simple autologistic model
create_plot(onov, colors_onov, fit1_onov, "Fit of the simple autologistic model - ONOV")

# Plot for ONOV with the autologistic model with water content
create_plot(onov, colors_onov, fit2_onov, "Fit of the autologistic model with water content - ONOV")

# Plot for ONOV with the autologistic model with substrate density
create_plot(onov, colors_onov, fit3_onov, "Fit of the autologistic model with substrate density - ONOV")

# Plot for ONOV with the autologistic model incorporating both water content and substrate density
create_plot(onov, colors_onov, fit4_onov, "Fit of the autologistic model incorporating both variables - ONOV")

# Plot for LCIL with the simple autologistic model
create_plot(lcil, colors_lcil, fit1_lcil, "Fit of the simple autologistic model - LCIL")

# Plot for LCIL with the autologistic model with water content
create_plot(lcil, colors_lcil, fit2_lcil, "Fit of the autologistic model with water content - LCIL")

# Plot for LCIL with the autologistic model with substrate density
create_plot(lcil, colors_lcil, fit3_lcil, "Fit of the autologistic model with substrate density - LCIL")

# Plot for LCIL with the autologistic model incorporating both water content and substrate density
create_plot(lcil, colors_lcil, fit4_lcil, "Fit of the autologistic model incorporating both variables - LCIL")


# 4. Verify the Predictive Ability of all models ----------------------------------------------

# First, we evaluate the models using confusion matrices

conf_matrix_onov1 <- table(onov$presence, fit1_onov)
print("Confusion Matrix for ONOV model 1:")
print(conf_matrix_onov1)
#   1
#0  7
#1 60

# The confusion matrix for ONOV Model 1 compares the predicted presence (1) and absence (0) with the actual presence (1) and absence (0) from the dataset
# shows 7 false negatives and 60 true positives, indicating good model performance in predicting species presence, though a few false negatives remain.

conf_matrix_onov2 <- table(onov$presence, fit2_onov)
print("Confusion Matrix for ONOV model 2:")
print(conf_matrix_onov2)
#    0  1
# 0  4  3
# 1  2 58

# The confusion matrix for ONOV Model 2 shows 4 false negatives, 3 false positives, 2 true negatives, 
# and 58 true positives, indicating good performance with a high number of correct presence predictions

conf_matrix_onov3 <- table(onov$presence, fit3_onov)
print("Confusion Matrix for ONOV model 3:")
print(conf_matrix_onov3)
#    0  1
# 0  1  6
# 1  1 59

# The confusion matrix for ONOV Model 3, which includes substrate density as a covariate, 
# shows 1 false negative, 6 false positives, 1 true negative, and 59 true positives. 
# The model performs similarly to the second, with slightly more true positives and fewer false negatives

conf_matrix_onov4 <- table(onov$presence, fit4_onov)
print("Confusion Matrix for ONOV model 4:")
print(conf_matrix_onov4)

#    0  1
# 0  4  3
# 1  2 58

# For ONOV Model 4, which includes both water content and substrate density as covariates, 
# the confusion matrix shows results similar to Model 2: 4 false negatives, 3 false positives, 
# 2 true negatives, and 58 true positives

# Again, we can comprare the models through the Akaike Information Criterion
AIC(onov1, onov2, onov3, onov4)
#          df      AIC
#   onov1  2 46.41517
#   onov2  3 27.29905
#   onov3  3 39.52509
#   onov4  4 28.90828

# This shows that overall, the best fitting model for ONIV is the one 
# with the added water content effect as the sole covariate, with the model with both covariates 
# having a close (though higher) AIC value


# For the LCIL models we obtain
conf_matrix_lcil1 <- table(lcil$presence, fit1_lcil)
print("Confusion Matrix for LCIL model 1:")
print(conf_matrix_lcil1)
#      0  1
#   0  2 11
#   1  3 51

# The confusion matrix for LCIL Model 1 shows 2 false negatives, 11 false positives, 3 true negatives, 
# and 51 true positives. The model has a higher number of false positives compared to true negatives

conf_matrix_lcil2 <- table(lcil$presence, fit2_lcil)
print("Confusion Matrix for LCIL model 2:")
print(conf_matrix_lcil2)
#     0  1
#  0  7  6
#  1  2 52

# The confusion matrix for LCIL Model 2, with water content as a covariate, shows 7 false negatives, 
# 6 false positives, 2 true negatives, and 52 true positives. 
# This model is more balanced than Model 1, with fewer false positives and negatives

conf_matrix_lcil3 <- table(lcil$presence, fit3_lcil)
print("Confusion Matrix for LCIL model 3:")
print(conf_matrix_lcil3)
#     0  1
#  0  7  6
#  1  2 52

# The confusion matrix for LCIL Model 3, incorporating substrate density as a covariate, 
# shows results identical to Model 2. This suggests that substrate density does not significantly 
# improve the model's performance

conf_matrix_lcil4 <- table(lcil$presence, fit4_lcil)
print("Confusion Matrix for LCIL model 4:")
print(conf_matrix_lcil4)
#    0  1
# 0  7  6
# 1  2 52

# The confusion matrix for LCIL Model 4, using both water content and substrate density as covariates, 
# yields results identical to Models 2 and 3, indicating no additional improvement in predictive performance

# comparison of the models for LCIL
AIC(lcil1, lcil2, lcil3, lcil4)

#         df      AIC
#  lcil1  2 58.72316
#  lcil2  3 49.17438
#  lcil3  3 60.48429
#  lcil4  4 50.83430

# Analogously to the ONOV case, the best fitting model is the one couting for water content effect,
# with the fourth model having a slightly higher AIC value



# We can also confirm out obtained results through other tools, such as the accuracy of the models
# First, we plot the accuracy for ONOV models (True Positives + True Negatives / Total Observations)
accuracy_onov1 <- sum(diag(conf_matrix_onov1)) / sum(conf_matrix_onov1)
cat("Accuracy for ONOV model 1:", accuracy_onov1, "\n") 
# 0.1044776

accuracy_onov2 <- sum(diag(conf_matrix_onov2)) / sum(conf_matrix_onov2)
cat("Accuracy for ONOV model 2:", accuracy_onov2, "\n") 
# 0.9253731

accuracy_onov3 <- sum(diag(conf_matrix_onov3)) / sum(conf_matrix_onov3)
cat("Accuracy for ONOV model 3:", accuracy_onov3, "\n") 
# 0.8955224


accuracy_onov4 <- sum(diag(conf_matrix_onov4)) / sum(conf_matrix_onov4)
cat("Accuracy for ONOV model 4:", accuracy_onov4, "\n") 
# 0.9241

# As expected, the second autologistic model has the highest accuracy (with model 4 having a very close value),
# whereas the fit for the simple autologistic model gives the worst results overall.


# Accuracy for LCIL models
accuracy_lcil1 <- sum(diag(conf_matrix_lcil1)) / sum(conf_matrix_lcil1)
cat("Accuracy for LCIL model 1:", accuracy_lcil1, "\n") 
# 0.7910448


accuracy_lcil2 <- sum(diag(conf_matrix_lcil2)) / sum(conf_matrix_lcil2)
cat("Accuracy for LCIL model 2:", accuracy_lcil2, "\n") 
# 0.880597

accuracy_lcil3 <- sum(diag(conf_matrix_lcil3)) / sum(conf_matrix_lcil3)
cat("Accuracy for LCIL model 3:", accuracy_lcil3, "\n") 
# 0.8358209

accuracy_lcil4 <- sum(diag(conf_matrix_lcil4)) / sum(conf_matrix_lcil4)
cat("Accuracy for LCIL model 4:", accuracy_lcil4, "\n") 
# 0.880597

# Again, the model with water content is the best fitting one (with model 4 having the exact same value) 
# with a satisfactory value close to 90% accuracy












# ONOV Models
roc_onov1 <- roc(onov$presence, fit1_onov)
roc_onov2 <- roc(onov$presence, fit2_onov)
roc_onov3 <- roc(onov$presence, fit3_onov)
roc_onov4 <- roc(onov$presence, fit4_onov)

# LCIL Models
roc_lcil1 <- roc(lcil$presence, fit1_lcil)
roc_lcil2 <- roc(lcil$presence, fit2_lcil)
roc_lcil3 <- roc(lcil$presence, fit3_lcil)
roc_lcil4 <- roc(lcil$presence, fit4_lcil)

# Plot ROC curves for ONOV
plot(roc_onov1, col = "blue", main = "ROC Curves for ONOV Models", 
     xlim = c(0, 1), ylim = c(0, 1), lwd = 2)
lines(roc_onov2, col = "red", lwd = 2)
lines(roc_onov3, col = "green", lwd = 2)
lines(roc_onov4, col = "purple", lwd = 2)
legend("bottomright", legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
       col = c("blue", "red", "green", "purple"), lwd = 2)

# Plot ROC curves for LCIL
plot(roc_lcil1, col = "blue", main = "ROC Curves for LCIL Models", 
     xlim = c(0, 1), ylim = c(0, 1), lwd = 2)
lines(roc_lcil2, col = "red", lwd = 2)
lines(roc_lcil3, col = "green", lwd = 2)
lines(roc_lcil4, col = "purple", lwd = 2)
legend("bottomright", legend = c("Model 1", "Model 2", "Model 3", "Model 4"),
       col = c("blue", "red", "green", "purple"), lwd = 2)

# Compute and display AUC for each model
cat("AUC for ONOV Model 1:", auc(roc_onov1), "\n")
cat("AUC for ONOV Model 2:", auc(roc_onov2), "\n")
cat("AUC for ONOV Model 3:", auc(roc_onov3), "\n")
cat("AUC for ONOV Model 4:", auc(roc_onov4), "\n")

cat("AUC for LCIL Model 1:", auc(roc_lcil1), "\n")
cat("AUC for LCIL Model 2:", auc(roc_lcil2), "\n")
cat("AUC for LCIL Model 3:", auc(roc_lcil3), "\n")
cat("AUC for LCIL Model 4:", auc(roc_lcil4), "\n")
