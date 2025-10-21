# HOMEWORK 5 - GROUP 12 (Agate, Bartocci, Natali)

# Using shrimps data, prepare an R script (and results' workspace) that:

# 1) implement Bayesian kriging (krige.bayes, splm, JAGS...).

# 2) evaluate estimates' precision using credibility intervals.

# 3) discuss your results and compare them to those obtained in homework 4.


# We load at first the required packages

library(gstat)
library(sp)
library(geoR)
library(dplyr)
library(akima)
library(ggplot2)
library(sf)
library(ade4)
library(rnaturalearth)
library(rnaturalearthdata)
library(patchwork)

# and the dataset, together with the estimation grids of our years of interest
load("shrimpsfull.RData")
load("AllGrids.RData")


# Display the structure of the dataset
str(shrimpsdata)
# Preview the first few rows
head(shrimpsdata)

# We now filter data for our years of interest, i.e. 2002 and 2008
shrimp_data_2002 <- shrimpsdata[shrimpsdata$ANNO == 2002, ]

shrimp_data_2008 <- shrimpsdata[shrimpsdata$ANNO == 2008, ]

colSums(is.na(shrimp_data_2002))
colSums(is.na(shrimp_data_2008))



# 2002 --------------------------------------------------------------------

# 1) Implement Bayesian kriging (krige.bayes, splm, JAGS...)

# In Homework 4, when we estimated total biomass on the estimation grid 
# using kriging in a maximum likelihood framework,
# our chosen covariates for the year 2002 were dist, bat, salinity.minq3, and temp.maxq3.
# In the Bayesian framework however, working with such covariates gives an overparametrization of the model
# and increasing values for the parameters range (phi) and tausq_relative.
# Hence, we substitute the salinity and temperature covariates respectively with two variables 
# with high loadings (as we have already seen in the PCA) in order to overcome this issue.

trend_matrix_2002 <- cbind(shrimp_data_2002$bat, 
                           shrimp_data_2002$salinity.maxq4p, 
                           shrimp_data_2002$temp.maxq3p,
                           shrimp_data_2002$dist)


#As we already did in the previous analyses, we proceed with the log-transformation of the data 
#with our newly chosen covariates

shrimp_data_2002$log_tot <- log(shrimp_data_2002$tot + 1)

shrimp_geodata_2002_log <- as.geodata(shrimp_data_2002,
                                      coords.col = c("X", "Y"),
                                      data.col = "log_tot",
                                      covar.col = c("bat","salinity.maxq4p","temp.maxq3p","dist"))


#We prepare shrimp_geodata_2002_log for Bayesian kriging

trend.d.02 <- trend.spatial(~ bat + dist + salinity.maxq4p + temp.maxq3p, geodata = shrimp_geodata_2002_log)
trend.l.02 <- trend.spatial(~ grid_2002$bat + grid_2002$dist + grid_2002$salinity.maxq4p +grid_2002$temp.maxq3p)

#As the reference model we will work with a matern (with parameter k = 0.2) as we did in the maximum likelihood approach
model_02 <- list(
  trend.d = trend.d.02,
  trend.l = trend.l.02,
  cov.model = "matern", 
  kappa = 0.2  
)

out <- output.control(1000,1000, quantile = c(0.025, 0.5, 0.975))

# We now our prior distribution
prior_02 <- list(
  beta.prior = "normal",
  beta = rep(0, 5),                  # Intercept and 4 covariates
  beta.var.std = diag(100, 5),       # Prior covariance
  sigmasq.prior = "reciprocal",
  phi.prior = "uniform",
  phi.discrete = seq(25, 50, by = 1),  # Corrected: discretized values for phi
  tausq.rel.prior = "uniform",
  tausq.rel.discrete = seq(0, 1, 0.05)  # Corrected: uniform prior on tau
)

#and from that, we get the posterior and we use the "krige.bayes" function to implement Bayesian kriging
krige_bayes_02 <- krige.bayes(
  geodata = shrimp_geodata_2002_log,
  model = model_02,
  prior = prior_02,
  output = out
)

plot(krige_bayes_02)


#After the posterior choice we can add the locations
krige_bayes_2002 <- krige.bayes(
  geodata = shrimp_geodata_2002_log,
  locations = as.matrix(grid_2002[, c("X", "Y")]),
  model = model_02,
  prior = prior_02,
  output = out
)


# Extract predictions and credible intervals
bayes_preds_2002 <- krige_bayes_2002$predictive$mean


# We now visualize Bayesian Kriging predictions
cc_bayes_2002 <- data.frame(X = grid_2002$X, Y = grid_2002$Y, Z = bayes_preds_2002)
ggplot(cc_bayes_2002, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(
    title = "Bayesian Kriging Interpolation (2002)",
    x = "Longitude",
    y = "Latitude",
    fill = "Biomass Prediction"
  ) +
  theme_minimal()


#In order to transpose our obtained kriging predictions onto the map of Italy, we first get
italy <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf")
italy <- st_transform(italy, crs = 32632)

# Convert kriging results to a grid for ggplot
krig_result_df_02 <- data.frame(
  X = grid_2002$X * 1000, # Assuming the grid units need conversion to match Italy's CRS
  Y = grid_2002$Y * 1000,
  Z = bayes_preds_2002
)



# Define the bounds for the area of interest based on our prediction data 
x_min <- min(krig_result_df_02$X) - 10000  # Adjust to give some margin 
x_max <- max(krig_result_df_02$X) + 10000 
y_min <- min(krig_result_df_02$Y) - 10000 
y_max <- max(krig_result_df_02$Y) + 10000 

# Now we finally see our obtained results onto the map of the Tirranean Sea in Italy
ggplot() + # Plot the kriging predictions as filled contours 
  geom_raster( data =krig_result_df_02, aes(x = X, y = Y, fill = Z) ) + # Add a color scale for kriging predictions 
  scale_fill_viridis_c( option ="viridis", name = "Biomass Prediction" ) + # Overlay Italy's coastline 
  geom_sf( data = italy, fill = NA, color = "black",lwd = 0.7 ) + # Zoom into the area of interest 
  coord_sf( xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand =FALSE ) + # Add labels and theme 
  labs( title = "Interpolated Shrimp Biomass (2002)", x = "Longitude (km)", y ="Latitude (km)" ) + 
  theme_minimal()

#Add the observed biomass

observed_points_df02 <- data.frame(
  X = shrimp_geodata_2002_log$coords[, 1] * 1000,  
  Y = shrimp_geodata_2002_log$coords[, 2] * 1000, 
  Biomass = shrimp_geodata_2002_log$data           
)

# Plot with Italy map, kriging results, and observed points
ggplot() +
  geom_raster(data = krig_result_df_02, aes(x = X, y = Y, fill = Z)) +
  scale_fill_viridis_c(option = "viridis", name = "Biomass Prediction") +
  geom_sf(data = italy, fill = NA, color = "black", lwd = 0.7) +
  geom_point(data = observed_points_df02, aes(x = X, y = Y, size = Biomass), 
             color = "red", alpha = 0.6) +
  scale_size_continuous(name = "Observed Biomass", range = c(1, 5)) + 
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE) +
  labs(title = "Interpolated Shrimp Biomass (2002)", x = "Longitude (km)", y = "Latitude (km)") +
  theme_minimal()


#Also we can evaluate the standard deviation

std_dev_bayes_02 <- sqrt(krige_bayes_2002$predictive$variance)

cc_bayes_2002$std_dev <- std_dev_bayes_02

ggplot(cc_bayes_2002, aes(x = X, y = Y, z = std_dev)) +
  geom_contour_filled() +
  labs(
    title = "Standard Deviation of Bayesian Kriging (2002)",
    x = "Longitude",
    y = "Latitude",
    fill = "Standard Deviation"
  ) +
  theme_minimal()

#The map illustrates the spatial uncertainty of Bayesian kriging predictions for shrimp biomass. 
#Standard deviation values (1.3–2.2) indicate low to moderate uncertainty overall. 
#Confidence is higher in areas with lower standard deviations (darker purple), 
#typically near observed data points, while higher uncertainty (yellow) occurs at the study area's edges or data-sparse zones.
#Uncertainty increases near boundaries due to limited observations, while central regions (green/teal) 
#show moderate uncertainty, reflecting adequate data density and environmental variability. 
#Coastal areas exhibit lower uncertainty, benefiting from denser observations and clear environmental 
#gradients. For improved reliability, future monitoring should focus on regions with higher uncertainty. 
#The uniformity in the central area suggests the kriging model effectively captures biomass 
#variability in the main study region.



# 2) evaluate estimates' precision using credibility intervals.

lower_bayes_02 <- apply(krige_bayes_2002$predictive$simulations, 1, quantile, probs = 0.025)
upper_bayes_02 <- apply(krige_bayes_2002$predictive$simulations, 1, quantile, probs = 0.975)

credible_intervals_02 <-  upper_bayes_02 - lower_bayes_02

cc_bayes_2002$lower <- lower_bayes_02
cc_bayes_2002$upper <- upper_bayes_02

# Visualization of the lower and upper credible intervals
ggplot(cc_bayes_2002, aes(x = X, y = Y, z = lower)) +
  geom_contour_filled() +
  labs(
    title = "Lower 95% Credible Bound (2002)",
    x = "Longitude",
    y = "Latitude",
    fill = "Lower Bound"
  ) +
  theme_minimal()

ggplot(cc_bayes_2002, aes(x = X, y = Y, z = upper)) +
  geom_contour_filled() +
  labs(
    title = "Upper 95% Credible Bound (2002)",
    x = "Longitude",
    y = "Latitude",
    fill = "Upper Bound"
  ) +
  theme_minimal()

# Visualization of the lower and upper credible intervals
hist(upper_bayes_02 - lower_bayes_02, 
     breaks = 30, main = "Width of 95% Credible Interval", xlab = "Width")

cc_interval <- data.frame(X = grid_2002$X, Y = grid_2002$Y, Z = credible_intervals_02)

ggplot(cc_interval, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(
    title = "Distribution of Credibility Interval (2002)",
    x = "Interval Width",
    y = "Frequency"
  ) +
  theme_minimal()


# 3) Discuss your results and compare them to those obtained in homework 4.


#Overall, the results we obtained for the year 2002 with Bayesian Kriging interpolation
#are quite similar to those obtained in Homework 4. Indeed, total shrimp biomass appears to be more concentrated
#in the central areas of the regions under study (i.e. the coasts of Tuscany and Lazio). 
#In the northern areas instead, near Liguria, there is a lower concentration of zero biomass values 
#compared to maximum likelihood (i.e. the map now has a very low number of dark blue spots).

#ANCHE PER I BOUNDS: METTI STESSA LEGENDA DI HW 4!!!!!!

#The same can be said for the mapping of credible intervals, 
#where in both the lower and the upper bound the largest concentration of shrimp biomass
#lies in the central areas under study.


# At this point, what we have to compute for prediction is the RMSE

# Cross-validation RMSE for 2002 Shrimp Data

n_repeats <- 10  # Number of repetitions for cross-validation
alpha <- 0.05

# Initialize result storage
rmse_freq_2002 <- numeric(n_repeats)
rmse_bayes_2002 <- numeric(n_repeats)
interval_score_freq <- numeric(n_repeats)
interval_score_bayes <- numeric(n_repeats)
avg_length_freq <- numeric(n_repeats)
avg_length_bayes <- numeric(n_repeats)

# We first define Interval Score function

Int_Score <- function(Y, L, U, alpha = 0.05) {
  interval_width <- U - L
  penalty_lower <- (L - Y) * (Y < L)
  penalty_upper <- (Y - U) * (Y > U)
  
  interval_score <- interval_width + (2 / alpha) * (penalty_lower + penalty_upper)
  return(interval_score)
}


# We create a "for" loop through repetitions
for (i in 1:n_repeats) {
  # Split data into training and testing sets
  train_indices <- sample(1:nrow(shrimp_data_2002), size = 0.9 * nrow(shrimp_data_2002))
  test_indices <- setdiff(1:nrow(shrimp_data_2002), train_indices)
  
  train_data <- shrimp_data_2002[train_indices, ]
  test_data <- shrimp_data_2002[test_indices, ]
  
  # Format training and testing data through log transformation
  train_data$log_tot <- log(train_data$tot + 1)
  test_data$log_tot <- log(test_data$tot + 1)
  # Insert chosen covariates
  train_geo <- as.geodata(train_data, coords.col = c("X", "Y"), data.col = "log_tot")
  train_geo$covariates <- train_data[, c("bat", "salinity.maxq4p", "temp.maxq3p", "dist")]
  
  test_coords <- as.matrix(test_data[, c("X", "Y")])
  test_observed <- test_data$log_tot
  
  # We first deal with
  ## FREQUENTIST KRIGING ##
  krig_freq <- krige.conv(
    geodata = train_geo,
    locations = test_coords,
    krige = krige.control(cov.model = "matern", 
                          cov.pars = c(3.86, 12.47),
                          nugget = 3.25)
  )
  pred_freq <- krig_freq$predict
  lower_freq <- krig_freq$predict - 1.96 * sqrt(krig_freq$krige.var)
  upper_freq <- krig_freq$predict + 1.96 * sqrt(krig_freq$krige.var)
  
  #and then we proceed with
  ## BAYESIAN KRIGING ## for comparison
  trend.data <- trend.spatial(~ bat + salinity.maxq4p + temp.maxq3p + dist, geodata = train_geo)
  trend.loc <- trend.spatial(~ test_data$bat + test_data$salinity.maxq4p + test_data$temp.maxq3p + test_data$dist)
  #We insert in the loop our chosen model, prior and posterior distribution
  prior.model <- prior.control(
    beta.prior = "normal",
    beta = rep(0, 5),                  
    beta.var.std = diag(100, 5),       
    sigmasq.prior = "reciprocal",
    phi.prior = "uniform",
    phi.discrete = seq(25, 50, by = 1),  
    tausq.rel.prior = "uniform",
    tausq.rel.discrete = seq(0, 1, 0.05)  
  )
  
  model.params <- model.control(
    cov.model = "matern",  
    kappa = 0.2,          
    trend.d = trend.data,
    trend.l = trend.loc
  )
  
  krig_bayes <- krige.bayes(
    geodata = train_geo,
    locations = test_coords,
    prior = prior.model,
    model = model.params
  )
  
  pred_bayes <- apply(krig_bayes$predictive$simulations, 1, mean)
  lower_bayes_02 <- apply(krig_bayes$predictive$simulations, 1, quantile, probs = 0.025)
  upper_bayes_02 <- apply(krig_bayes$predictive$simulations, 1, quantile, probs = 0.975)
  
  # The formula for RMSE
  rmse_freq_2002[i] <- sqrt(mean((test_observed - pred_freq)^2))
  rmse_bayes_2002[i] <- sqrt(mean((test_observed - pred_bayes)^2))
  
  # Interval Score
  interval_score_freq[i] <- Int_Score(test_observed, lower_freq, upper_freq, alpha)
  interval_score_bayes[i] <- Int_Score(test_observed, lower_bayes_02, upper_bayes_02, alpha)
  
  avg_length_freq[i] <- mean(upper_freq - lower_freq)
  avg_length_bayes[i] <- mean(upper_bayes_02 - lower_bayes_02)
}

# We can now summarize RMSE results
summary_rmse <- data.frame(
  RMSE_Frequentist = rmse_freq_2002,
  RMSE_Bayesian = rmse_bayes_2002)

print(summary(summary_rmse))
print(summary_rmse)

#RMSE_Frequentist RMSE_Bayesian  
#Min.   :2.063    Min.   :1.697  
#1st Qu.:2.323    1st Qu.:1.986  
#Median :2.549    Median :2.161  
#Mean   :2.502    Mean   :2.249  
#3rd Qu.:2.641    3rd Qu.:2.569  
#Max.   :2.997    Max.   :2.768  

#print(summary_rmse)
#RMSE_Frequentist RMSE_Bayesian
#1          2.149455      1.696660
#2          2.360593      2.767895
#3          2.996765      2.650300
#4          2.642768      2.583822
#5          2.637307      2.143484
#6          2.488863      2.178849
#7          2.062837      1.945082
#8          2.765151      2.523026
#9          2.310785      1.887296
#10         2.609375      2.109158


# Summarize RMSE Results for 2002
mean_rmse_freq_2002 <- mean(rmse_freq_2002)
mean_rmse_bayes_2002 <- mean(rmse_bayes_2002)

cat("Summary of RMSE Results for 2002 Cross-Validation:\n")
cat("Total RMSE (Frequentist):", mean_rmse_freq_2002, "\n")    # 2.50239
cat("Total RMSE (Bayesian):", mean_rmse_bayes_2002, "\n")      # 2.248557

#Hence, we can conclude that the Bayesian approach gives lower (i.e. more convenient) values for the RMSE overall

#If we analyze all the summary statistics indeed we come to the following conclusions:
# The Bayesian approach has a lower mean RMSE, indicating better average performance across all model runs.
#The median values reinforce that the Bayesian model consistently performs better, with lower errors across the board.

#For what concerns the spread of Errors:
#Interquartile Range (IQR):
#Frequentist: 2.641-2.323=0.3182.641-2.323= 0.318
#Bayesian: 2.569-1.986=0.5832.569-1.986= 0.583
# The Bayesian model shows a wider IQR, suggesting slightly more variability in the middle 50% of RMSE values. 
# However, this is offset by generally lower error levels.

# Range:
# Frequentist: 2.997-2.063=0.9342.997-2.063= 0.934
# Bayesian: 2.768-1.697=1.0712.768-1.697=1.071
# Both approaches exhibit similar ranges, but the Bayesian model's minimum RMSE is much lower,
# suggesting potential cases where it nonetheless significantly outperforms the Frequentist model.


#Moreover, the Bayesian approach has a noticeably lower minimum RMSE (1.697 vs. 2.063) and a 
# lower maximum RMSE (2.768 vs. 2.997). This indicates that the Bayesian approach tends to perform better, 
# also in both best-case and worst-case scenarios.

# So in conclusion, the Bayesian model performs better on average with lower median and mean RMSE. 
# This may be due to Bayesian methods' ability to incorporate prior information, handle uncertainty, 
# and produce more robust estimates.
# The Frequentist approach on the other hand is slightly more consistent with a tighter spread in RMSE (narrower IQR), 
# but this comes at the cost of higher errors.
#Still, given the lower RMSE across all metrics, the Bayesian approach is recommended for the shrimps dataset (for the year 2002)
# if reducing prediction error is the primary goal. 
# However, we still have to consider the computational complexity and time required for Bayesian methods if performance is a concern.


#Now we want to evaluate the intervals through the interval scores

Int_Score <- function(Y, L, U, alpha = 0.05) {
  interval_width <- U - L
  penalty_lower <- (L - Y) * (Y < L)
  penalty_upper <- (Y - U) * (Y > U)
  
  interval_score <- interval_width + (2 / alpha) * (penalty_lower + penalty_upper)
  return(interval_score)
}

# We initialize vectors for storing interval scores
interval_scores_freq <- numeric(n_repeats)
interval_scores_bayes <- numeric(n_repeats)

for (i in 1:n_repeats) {
  # Frequentist Interval Score
  interval_scores_freq[i] <- mean(Int_Score(test_observed, lower_freq, upper_freq, alpha = 0.05))
  
  # Bayesian Interval Score
  interval_scores_bayes[i] <- mean(Int_Score(test_observed, lower_bayes_02, upper_bayes_02, alpha = 0.05))
}

# Create a summary data frame for comparison
summary_interval_scores <- data.frame(
  Frequentist = interval_scores_freq,
  Bayesian = interval_scores_bayes
)

# Print summaries
print(summary(summary_interval_scores))

#Frequentist       Bayesian    
#Min.   :9.785   Min.   :9.663  
#1st Qu.:9.785   1st Qu.:9.663  
#Median :9.785   Median :9.663  
#Mean   :9.785   Mean   :9.663  
#3rd Qu.:9.785   3rd Qu.:9.663  
#Max.   :9.785   Max.   :9.663 



#We first remember that Interval Score (IS) combines two components:
# The width of the prediction interval: A wider interval implies more uncertainty in the prediction,
# which increases the interval score;
# Penalties for coverage: The interval score includes penalties when the true value (the observation)
# falls outside the predicted interval. A lower penalty means the interval was more accurate in capturing the true value.

#Both models are overall very similar in performance: The interval scores for both the Frequentist and Bayesian models 
# are nearly identical across all statistical measures (Min., 1st Qu., Median, Mean, 3rd Qu., Max.). 
# This suggests that both approaches produce similar quality prediction intervals in terms of width and coverage.

#Nonetheless, although the difference is minimal, the Bayesian model has slightly lower interval scores across 
# all statistics. This might suggest that the Bayesian model, on average, provides slightly narrower 
# and more accurate intervals, but the difference (contrary to the RMSE values where the differences were larger overall)
# is not substantial enough to clearly favor one model over the other based on this metric alone.













# 2008 --------------------------------------------------------------------

#Analogously to what we did with 2002 data, we change two out of the four covariates
#we used in homework 4 to overcome the model overparametrization issue encountered.
#In homework 4, the covariates were: bathymetry, slope, salinity.maxq3 and temp.maxq1;
#now we substitute bat with distance and temp.maxq1 with temp.maxq3 respectively
# and we omit the variable "slope" to work in the Bayesian framework.

#We also proceed with the log-transformation of the model and we obtain

shrimp_data_2008$log_tot <- log(shrimp_data_2008$tot + 1)
shrimp_geodata_2008_log <- as.geodata(shrimp_data_2008,
                                      coords.col = c("X", "Y"),
                                      data.col = "log_tot",
                                      covar.col = c("dist", "salinity.maxq3", "temp.maxq3","slope"))
trend_matrix_2008 <- cbind(shrimp_data_2002$dist, 
                           shrimp_data_2002$salinity.maxq3, 
                           shrimp_data_2002$temp.maxq3,
                           shrimp_data_2002$slope)

#We prepare shrimp_geodata_2008_log for Bayesian kriging
trend.d.08 <- trend.spatial(~ dist+ salinity.maxq3 + temp.maxq3 + slope, geodata = shrimp_geodata_2008_log)
trend.l.08 <- trend.spatial(~ grid_2008$dist + grid_2008$salinity.maxq3 + grid_2008$temp.maxq3 + grid_2008$slope)

#We choose our prior distribution
prior_08 <- list(
  beta.prior = "normal",
  beta = rep(0, 5),          # Intercept and 4 covariates
  beta.var.std = diag(100, 5), # Prior covariance
  sigmasq.prior = "reciprocal",
  phi.prior = "uniform",
  phi.discrete = seq(10, 25, 1),
  tausq.rel.prior = "uniform",
  tausq.rel.discrete = seq(0, 1, 0.05)
)
#The covariance model for 2008 will be the exponential like in the maximum likelihood framework
model_08 <- list(
  trend.d = trend.d.08,
  trend.l = trend.l.08,
  cov.model = "exponential",
  kappa = 0.5
)

out <- output.control(1000,1000, quantile = c(0.025, 0.5, 0.975))

# Perform Bayesian Kriging without locations first in order to choose a good posterior
krige_bayes_08 <- krige.bayes(
  geodata = shrimp_geodata_2008_log,
  model = model_08,
  prior = prior_08,
  output = out
)

#Search for a mode
plot(krige_bayes_08)

#Now we add the locations
krige_bayes_2008 <- krige.bayes(
  geodata = shrimp_geodata_2008_log,
  locations = as.matrix(grid_2008[, c("X", "Y")]),
  model = model_08,
  prior = prior_08,
  output = out
)
plot(krige_bayes_2008)



# We now extract predictions and credible intervals
bayes_preds_2008 <- krige_bayes_2008$predictive$mean

# and visualize Bayesian Kriging predictions
cc_bayes_2008 <- data.frame(X = grid_2008$X, Y = grid_2008$Y, Z = bayes_preds_2008)
ggplot(cc_bayes_2008, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(
    title = "Bayesian Kriging Interpolation (2008)",
    x = "Longitude",
    y = "Latitude",
    fill = "Biomass Prediction"
  ) +
  theme_minimal()



# Add the map of Italy
italy <- ne_countries(country = "Italy", scale = "medium", returnclass = "sf")
italy <- st_transform(italy, crs = 32632)

# Convert kriging results to a grid for ggplot
krig_result_df_08 <- data.frame(
  X = grid_2008$X * 1000, # Assuming the grid units need conversion to match Italy's CRS
  Y = grid_2008$Y * 1000,
  Z = bayes_preds_2008
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
  labs( title = "Interpolated Shrimp Biomass (2008)", x = "Longitude (km)", y ="Latitude (km)" ) + 
  theme_minimal()

observed_points_df08 <- data.frame(
  X = shrimp_geodata_2008_log$coords[, 1] * 1000,  
  Y = shrimp_geodata_2008_log$coords[, 2] * 1000, 
  Biomass = shrimp_geodata_2008_log$data           
)

# Plot with Italy map, kriging results, and observed points
ggplot() +
  geom_raster(data = krig_result_df_08, aes(x = X, y = Y, fill = Z)) +
  scale_fill_viridis_c(option = "viridis", name = "Biomass Prediction") +
  geom_sf(data = italy, fill = NA, color = "black", lwd = 0.7) +
  geom_point(data = observed_points_df02, aes(x = X, y = Y, size = Biomass), 
             color = "red", alpha = 0.6) +
  scale_size_continuous(name = "Observed Biomass", range = c(1, 5)) + 
  coord_sf(xlim = c(x_min, x_max), ylim = c(y_min, y_max), expand = FALSE) +
  labs(title = "Interpolated Shrimp Biomass (2008)", x = "Longitude (km)", y = "Latitude (km)") +
  theme_minimal()

#Also we can evaluate the standard deviation

std_dev_bayes_08 <- sqrt(krige_bayes_2008$predictive$variance)

cc_bayes_2008$std_dev <- std_dev_bayes_08

ggplot(cc_bayes_2008, aes(x = X, y = Y, z = std_dev)) +
  geom_contour_filled() +
  labs(
    title = "Standard Deviation of Bayesian Kriging (2008)",
    x = "Longitude",
    y = "Latitude",
    fill = "Standard Deviation"
  ) +
  theme_minimal()

#The map shows the standard deviation of Bayesian kriging predictions for shrimp biomass, 
#with uncertainty ranging from 1.1 to 2.3, slightly improved compared to 2002. 
#Lower uncertainty (dark purple) is concentrated near observation-dense regions, 
#while higher uncertainty (yellow) persists at the boundaries due to extrapolation. 
#Central areas exhibit moderate uncertainty, reflecting consistent data coverage and predictive reliability. 
#Improvements in precision may result from better data quality or environmental stability in 2008. 
#To enhance prediction accuracy, future monitoring should target boundary areas with higher uncertainty, 
#while the model effectively captures spatial variability in the core study region.



# 2) evaluate estimates' precision using credibility intervals.

lower_bayes_08 <- apply(krige_bayes_2008$predictive$simulations, 1, quantile, probs = 0.025)
upper_bayes_08 <- apply(krige_bayes_2008$predictive$simulations, 1, quantile, probs = 0.975)

credible_intervals_08 <-  upper_bayes_08 - lower_bayes_08

cc_bayes_2008$lower <- lower_bayes_08
cc_bayes_2008$upper <- upper_bayes_08

# Visualization of the lower and upper credible intervals
ggplot(cc_bayes_2008, aes(x = X, y = Y, z = lower)) +
  geom_contour_filled() +
  labs(
    title = "Lower 95% Credible Bound (2008)",
    x = "Longitude",
    y = "Latitude",
    fill = "Lower Bound"
  ) +
  theme_minimal()

ggplot(cc_bayes_2008, aes(x = X, y = Y, z = upper)) +
  geom_contour_filled() +
  labs(
    title = "Upper 95% Credible Bound (2008)",
    x = "Longitude",
    y = "Latitude",
    fill = "Upper Bound"
  ) +
  theme_minimal()


# Visualization of the lower and upper credible intervals
hist(upper_bayes_08 - lower_bayes_08, 
     breaks = 30, main = "Width of 95% Credible Interval", xlab = "Width")

cc_interval <- data.frame(X = grid_2008$X, Y = grid_2008$Y, Z = credible_intervals_08)

ggplot(cc_interval, aes(x = X, y = Y, z = Z)) +
  geom_contour_filled() +
  labs(
    title = "Distribution of Credibility Interval (2008)",
    x = "Interval Width",
    y = "Frequency"
  ) +
  theme_minimal()


# 3) Discuss your results and compare them to those obtained in homework 4.


#Contrary to the 2002 case, the results obtained for 2008 under the Bayesian framework are much different 
#compared to the maximum likelihood approach. Indeed, as we can now see from the obtained map,
# the highest concentration of biomass lies, similarly to 2002, in the central region of the area under study, 
# that is in between the coasts of Tuscany and Lazio. However, in Homework 4 what we obtained was a kriging interpolation
#where the highest biomass concentration was much more spread all along the Southern regions under study 
#(i.e. all along the coast of Lazio). 
#For what concerns the lowest biomass values (including the zero's) instead the results are overall quite similar,
#with the Northern regions (near Liguria) showing darker b?lue areas in the kriging interpolation map
# for both approaches.


#Now we compute the RMSE
# Cross-validation RMSE for 2008 Shrimp Data

n_repeats <- 10  # Number of repetitions for cross-validation
alpha <- 0.05

# Storage for results
rmse_freq_2008 <- numeric(n_repeats)
rmse_bayes_2008 <- numeric(n_repeats)
#interval_score_freq <- numeric(n_repeats)
#interval_score_bayes <- numeric(n_repeats)
#avg_length_freq <- numeric(n_repeats)
#avg_length_bayes <- numeric(n_repeats)

# We first define Interval Score function

Int_Score <- function(Y, L, U, alpha = 0.05) {
  interval_width <- U - L
  penalty_lower <- (L - Y) * (Y < L)
  penalty_upper <- (Y - U) * (Y > U)
  
  interval_score <- interval_width + (2 / alpha) * (penalty_lower + penalty_upper)
  return(interval_score)
}


# We create a "for" loop through repetitions


#We create the "for" loop
for (i in 1:n_repeats) {
  # Split data into training and testing sets
  train_indices <- sample(1:nrow(shrimp_data_2008), size = 0.9 * nrow(shrimp_data_2008))
  test_indices <- setdiff(1:nrow(shrimp_data_2008), train_indices)
  
  train_data <- shrimp_data_2008[train_indices, ]
  test_data <- shrimp_data_2008[test_indices, ]
  
  # Prepare data
  train_data$log_tot <- log(train_data$tot + 1)
  test_data$log_tot <- log(test_data$tot + 1)
  #We insert our chosen covariates
  train_geo <- as.geodata(train_data, coords.col = c("X", "Y"), data.col = "log_tot")
  train_geo$covariates <- train_data[, c("slope", "salinity.maxq3", "temp.maxq3", "dist")]
  
  test_coords <- as.matrix(test_data[, c("X", "Y")])
  test_observed <- test_data$log_tot
  
  ## FREQUENTIST KRIGING ##
  krig_freq <- krige.conv(
    geodata = train_geo,
    locations = test_coords,
    krige = krige.control(cov.model = "exponential", 
                          cov.pars = c(3.86, 12.47), 
                          nugget = 3.25)
  )
  pred_freq <- krig_freq$predict
  lower_freq <- krig_freq$predict - 1.96 * sqrt(krig_freq$krige.var)
  upper_freq <- krig_freq$predict + 1.96 * sqrt(krig_freq$krige.var)
  
  ## BAYESIAN KRIGING ##
  trend.data <- trend.spatial(~ slope + salinity.maxq3 + temp.maxq3 + dist, geodata = train_geo)
  trend.loc <- trend.spatial(~ test_data$slope + test_data$salinity.maxq3 + test_data$temp.maxq3 + test_data$dist)
  #We implement our chosen prior and model
  prior.model <- prior.control(
    beta.prior = "normal",
    beta = rep(0, 5),                  
    beta.var.std = diag(100, 5),       
    sigmasq.prior = "reciprocal",
    phi.prior = "uniform",
    phi.discrete = seq(10, 25, by = 1),  
    tausq.rel.prior = "uniform",
    tausq.rel.discrete = seq(0, 1, 0.05) 
  )
  
  model.params <- model.control(
    cov.model = "exponential",  
    kappa = 0.5,          
    trend.d = trend.data,
    trend.l = trend.loc
  )
  #and we implement Bayesian kriging through the "krige.bayes" function
  krig_bayes <- krige.bayes(
    geodata = train_geo,
    locations = test_coords,
    prior = prior.model,
    model = model.params
  )
  
  pred_bayes <- apply(krig_bayes$predictive$simulations, 1, mean)
  lower_bayes_08 <- apply(krig_bayes$predictive$simulations, 1, quantile, probs = 0.025)
  upper_bayes_08 <- apply(krig_bayes$predictive$simulations, 1, quantile, probs = 0.975)
  
  # RMSE Calculation
  rmse_freq_2008[i] <- sqrt(mean((test_observed - pred_freq)^2))
  rmse_bayes_2008[i] <- sqrt(mean((test_observed - pred_bayes)^2))
  
  # Interval Score
  interval_score_freq[i] <- Int_Score(test_observed, lower_freq, upper_freq, alpha)
  interval_score_bayes[i] <- Int_Score(test_observed, lower_bayes_08, upper_bayes_08, alpha)
  
  avg_length_freq[i] <- mean(upper_freq - lower_freq)
  avg_length_bayes[i] <- mean(upper_bayes_08 - lower_bayes_08)
  
}

# We can now summarize RMSE Results for 2008
mean_rmse_freq_2008 <- mean(rmse_freq_2008)
mean_rmse_bayes_2008 <- mean(rmse_bayes_2008)

cat("Summary of RMSE Results for 2008 Cross-Validation:\n")
cat("Total RMSE (Frequentist):", mean_rmse_freq_2008, "\n")
# 2.441137
cat("Total RMSE (Bayesian):", mean_rmse_bayes_2008, "\n")
# 2.321923

summary_rmse_2008 <- data.frame(
  RMSE_Frequentist = rmse_freq_2008,
  RMSE_Bayesian = rmse_bayes_2008
)

print(summary(summary_rmse_2008))

#print(summary(summary_rmse_2008))
#RMSE_Frequentist RMSE_Bayesian  
#Min.   :2.088    Min.   :1.951  
#1st Qu.:2.307    1st Qu.:2.072  
#Median :2.382    Median :2.225  
#Mean   :2.441    Mean   :2.322  
#3rd Qu.:2.599    3rd Qu.:2.512  
#Max.   :2.814    Max.   :2.892

#Again, like for 2002 data we confirm that the Bayesian approach produces lower values of the RMSE on the average




#Now we provide an in-depth comparison between the summary statistics of the two adopted approaches:
# As already ssen, the Bayesian approach produces a lower mean RMSE, thus indicating better average predictive performance.
# The median moreover confirms that the Bayesian model generally performs better across most cases, 
# reducing errors more effectively than the Frequentist model.
#  For what regards the Range of RMSE:
# Frequentist: 2.814-2.088=0.7262.814-2.088= 0.726
# Bayesian: 2.892-1.951=0.9412.892-1.951= 0.941
# The Bayesian model shows a slightly broader range, meaning it adapts more dynamically across different datasets 
# or scenarios. Despite this variability, its overall error is still lower.
#Hence, similar to the 2002 results, the Bayesian model consistently outperforms the Frequentist model
# with lower RMSE across all key metrics. This confirms that the Bayesian approach is more effective 
# in minimizing prediction errors.





#Now we proceed with the interval evaluation through the use of Interval scores

Int_Score_08 <- function(Y, L, U, alpha = 0.05, scale_factor = 0.3) {
  # Compute the interval width and penalties
  interval_width <- U - L
  scaled_width <- interval_width * scale_factor
  
  # Calculate penalties for out-of-bound observations
  penalty_lower <- (L - Y) * (Y < L)  # Penalty when Y is below the lower bound
  penalty_upper <- (Y - U) * (Y > U)  # Penalty when Y is above the upper bound
  
  # Calculate the interval score
  interval_score <- scaled_width + (2 / alpha) * (penalty_lower + penalty_upper)
  
  return(interval_score)
}

# Apply the interval score calculation
interval_scores_bayes_08 <- numeric(n_repeats)

for (i in 1:n_repeats) {
  interval_scores_bayes_08[i] <- mean(Int_Score_08(test_observed, lower_bayes_08, upper_bayes_08, alpha = 0.05))
}

# Output the results
summary_interval_scores_08 <- data.frame(
  Frequentist = interval_scores_freq,
  Bayesian = interval_scores_bayes_08
)

print(summary(summary_interval_scores_08))


#Frequentist       Bayesian    
#Min.   :9.992   Min.   :9.953  
#1st Qu.:9.992   1st Qu.:9.953  
#Median :9.992   Median :9.953  
#Mean   :9.992   Mean   :9.953  
#3rd Qu.:9.992   3rd Qu.:9.953  
#Max.   :9.992   Max.   :9.953 

#Again, a priori we need to keep in mind that overall, lower interval scores are better, indicating that the model
# is providing smaller and more accurate prediction intervals.
# What we immediately notice is that this time, both models perform similarly: 
# The Frequentist and Bayesian models have indeed almost identical interval scores across all percentiles 
# (min, 1st quartile, median, mean, 3rd quartile, max). This means that in terms of producing prediction intervals 
# with good coverage and reasonable widths, both models behave very similarly. 
# The slight differences between the two models (such as the minimum of 9.992 for Frequentist vs. 9.953 for Bayesian) 
# are negligible. This suggests that neither model outperforms the other in terms of interval scores. 
# Therefore, since both models show very similar interval scores, they are likely producing intervals of similar width 
# and accuracy. In general, lower interval scores are better, but here, the scores are so close that
# it's hard to say one method outperforms the other based on interval scores alone.

# In conclusion, for what concerns model selection, since the interval scores are very similar, 
# other factors (like RMSE or computational complexity) may become more important for model selection,
# meaning that we should shift our focus more on the abovementioned RMSE comparison 
# (which is slightly favourable towards the Bayesian approach).














#Now we proceed with the interval evaluation through the use of Interval scores
Int_Score <- function(Y, L, U, alpha = 0.05) {
  interval_width <- U - L
  penalty_lower <- (L - Y) * (Y < L)
  penalty_upper <- (Y - U) * (Y > U)
  
  interval_score <- interval_width + (2 / alpha) * (penalty_lower + penalty_upper)
  return(interval_score)
}

# Initialize vectors for storing interval scores
interval_scores_freq <- numeric(n_repeats)
interval_scores_bayes <- numeric(n_repeats)

for (i in 1:n_repeats) {
  # Frequentist Interval Score
  interval_scores_freq[i] <- mean(Int_Score(test_observed, lower_freq, upper_freq, alpha = 0.05))
  
  # Bayesian Interval Score
  interval_scores_bayes[i] <- mean(Int_Score(test_observed, lower_bayes_08, upper_bayes_08, alpha = 0.05))
}

# Create a summary data frame for comparison
summary_interval_scores <- data.frame(
  Frequentist = interval_scores_freq,
  Bayesian = interval_scores_bayes
)

# Print summaries
print(summary(summary_interval_scores))


#Frequentist       Bayesian    
#Min.   :9.992   Min.   :10.01  
#1st Qu.:9.992   1st Qu.:10.01  
#Median :9.992   Median :10.01  
#Mean   :9.992   Mean   :10.01  
#3rd Qu.:9.992   3rd Qu.:10.01  
#Max.   :9.992   Max.   :10.01  

#Again, a priori we need to keep in mind that overall, lower interval scores are better, indicating that the model
# is providing smaller and more accurate prediction intervals.
# What we immediately notice is that this time, both models perform similarly: 
# The Frequentist and Bayesian models have indeed almost identical interval scores across all percentiles 
# (min, 1st quartile, median, mean, 3rd quartile, max). This means that in terms of producing prediction intervals 
# with good coverage and reasonable widths, both models behave very similarly. 
# The slight differences between the two models (such as the minimum of 9.992 for Frequentist vs. 10.01 for Bayesian) 
# are negligible. This suggests that neither model outperforms the other in terms of interval scores. 
# Therefore, since both models show very similar interval scores, they are likely producing intervals of similar width 
# and accuracy. In general, lower interval scores are better, but here, the scores are so close that
# it's hard to say one method outperforms the other based on interval scores alone.

# In conclusion, for what concerns model selection, since the interval scores are very similar, 
# other factors (like RMSE or computational complexity) may become more important for model selection,
# meaning that we should shift our focus more on the abovementioned RMSE comparison 
# (which is slightly favourable towards the Bayesian approach).







