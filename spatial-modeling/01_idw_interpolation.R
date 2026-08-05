# Spatial Statistics and Statistical Tools for Environmental Data
# MSc in Statistical Methods and Applications, Sapienza University of Rome
# Group 12 — Agate, Bartocci, Natali
#
# Homework 1 — Inverse Distance Weighting (IDW) interpolation
Choosing the smoothing parameter p by spatial cross-validation (RMSE), using gstat/sp.

# PART A ------------------------------------------------------------------

# a) Build an R function to choose the parameter p in the inverse distance weighting estimator.
# Hint: use the idw implemented in gstat.



# When we have not made any assumptions regarding our spatial data and the appropriate model,
# a useful method to explore data is to plot surfaces based on the IDW interpolator:
# the Inverse Distance estimator (IDW) is a weighted mean with weights proportional to the inverse of the distance
# between unobserved and observed locations.
# The rationale behind this interpolator is Tobler's Law,
# i.e. closer locations are more correlated than distant ones.
# The weights between observed and unobserved locations (s_0, s_1) are given by
# w(s_0, s_1) = 1/[dist(s_0, s_1)^p], where "p" is the smoothing parameter of the function:
# the larger p, the larger the spatial clusters hence the smoother the surface.

# I first load the required packages
require(gstat) #the "idw" function is stored in this package
require(sp)


# For the IDW interpolator, a proper method we can use is cross-validation for spatial data
# in order to choose the optimal smoothing parameter p.
# After that, we will evaluate our prediction through an appropriate score, 
# that is, the Root Mean Square Error (RMSE)

choose_best_p <- function(observations, p_values = seq(0.5, 4, 0.5), nmax = 12) {
  
  # An important assumption to make is that we are working with spatial data (i.e. a spatial data frame)
  if (!inherits(observations, "SpatialPointsDataFrame")) {
    stop("Input must be a SpatialPointsDataFrame")
  }
  
  # Extract coordinates and the variable of interest
  coords <- coordinates(observations)
  z <- observations$z # Assuming the column 'z' holds the variable to interpolate, 
  #that is, z(s) is our spatial data of interest
  
  # Store the RMSE for each 'p' value
  rmse_values <- numeric(length(p_values))
  
  # Loop over p values
  for (i in seq_along(p_values)) {
    p <- p_values[i]
    
    # To work with cross-validation, we apply the "Leave one out" approach:
    # I remove one of the points from the set, I then interpolate this value with the others
    # before recording the estimate x_hat and the observed value x for every point.
    # I then compute the (root) MSE using these two values.
    # This method works when we are dealing with dependent data.
    predictions <- numeric(length(z))
    for (j in 1:length(z)) {
      # Training data excluding a generic point j
      train_data <- observations[-j,]
      
      # Perform IDW with current 'p' value
      idw_model <- idw(z ~ 1, train_data, newdata = observations[j,], idp = p, nmax = nmax)
      
      # We store the obtained prediction
      predictions[j] <- idw_model$var1.pred
    }
    
    # We finally compute the RMSE for current value of 'p'
    rmse_values[i] <- sqrt(mean((z - predictions)^2))
  }
  
  # We find the 'p' that minimizes RMSE. 
  # Let us also underline that the value of 'p' that we will choose is the one which minimizes also the cross validation coefficient CV 
  best_p <- p_values[which.min(rmse_values)]
  
  # Finally we return the optimal 'p' and the corresponding RMSE
  return(list(best_p = best_p, rmse = rmse_values[which.min(rmse_values)]))
}


# PART B ------------------------------------------------------------------

#b) Work on the Wolfcamp dataset. Interpolate Wolfcamp data on a 20x20 grid (built from the coordinates) 
# using IDW and multilevel bi-splines (function mba.surf from MBA) and compare the results using a proper score.



# We first load the required packages
require(geoR)
require(MBA)

# As a first step, we load the wolfcamp dataset from the geoR package
data("wolfcamp")

# A summary of the dataset
summary(wolfcamp)

# We can print the first few rows of the data set
head(wolfcamp)

# Display the structure of the data
str(wolfcamp)

# We notice that the data set contains 85 observations
n <- nrow(wolfcamp$coords)
n #85

# We can also plot the data to visualize the spatial points also with an histogram
plot(wolfcamp)

# Convert the Wolfcamp data into a SpatialPointsDataFrame (spatial format)
coords <- cbind(wolfcamp$coords[,1], wolfcamp$coords[,2])  # Extract coordinates (x, y)
wolfcamp_sp <- SpatialPointsDataFrame(coords = coords, data = data.frame(z = wolfcamp$data))  # Create spatial object

## geoR spatial interpolation:
# For the wolfcamp dataset, we can see a spatial trend (function of coordinates)
# which means that the only covariates we need are the Coordinates.
# From such coordinates we now build a 20 x 20 grid to interpolate the Wolfcamp data.
# We first create two vectors for the two axes using a sequence of 20 coordinates each, 
# that lie between the minimum observed and the maximum observed.
x <- seq(min(wolfcamp$coords[,1]), max(wolfcamp$coords[,1]), length = 20)
y <- seq(min(wolfcamp$coords[,2]), max(wolfcamp$coords[,2]), length = 20)

# Create a grid by expanding the x and y coordinates into a full grid
grid <- expand.grid(x, y)
plot(grid)  # Visualize the grid (20x20 cells)

# Convert the grid to SpatialPoints object
coordinates(grid) <- ~Var1+Var2  # 'Var1' and 'Var2' represent the x and y grid points
gridded(grid) <- TRUE  # Mark the grid as a gridded spatial object

# Apply the IDW interpolation using the best 'p' value obtained from cross-validation
best_p <- choose_best_p(wolfcamp_sp)$best_p  # Get the best 'p' value by running the function
idw_result <- idw(z ~ 1, wolfcamp_sp, newdata = grid, idp = best_p)  # Perform IDW with the best 'p'

# An interpolator category is the one of "not fully modeled", where uncertainty can be evaluated 
# only by cross validation or resampling methods such as bootstrap.
# Approaches belonging to this category are for instance multilevel bi-splines
# and several types of weighted averages, like the IDW, respectively.

# Prepare the data for MBA interpolation (needs data in (x, y, z) format)
xyz <- cbind(wolfcamp$coords[,1], wolfcamp$coords[,2], wolfcamp$data)  # Combine coordinates and data

# Apply the MBA (multilevel B-splines) interpolation with a 20x20 grid
mba_result <- mba.surf(xyz, no.X = 20, no.Y = 20, extend = TRUE)  # 'extend = TRUE' to avoid NAs on the edges

# Convert the MBA result into a SpatialPixelsDataFrame for comparison
mba_grid <- expand.grid(mba_result$xyz.est$x, mba_result$xyz.est$y)  # Grid of MBA interpolated points
mba_sp <- SpatialPixelsDataFrame(mba_grid, data = data.frame(z = as.vector(mba_result$xyz.est$z)))  # Create spatial object

# Extract the observed values from the Wolfcamp data (for RMSE calculation)
observed <- wolfcamp$data

# Extract the predicted values from the IDW result
idw_predictions <- over(wolfcamp_sp, idw_result)$var1.pred  # IDW predicted values for each point

# Extract the predicted values from the MBA result
mba_predictions <- over(wolfcamp_sp, mba_sp)$z  # MBA predicted values for each point

# Compute the RMSE for the IDW method
rmse_idw <- sqrt(mean((observed - idw_predictions)^2))

# Compute the RMSE for the MBA method
rmse_mba <- sqrt(mean((observed - mba_predictions)^2))

# Print out the RMSE values for both methods (IDW and MBA)
cat("RMSE for IDW: ", rmse_idw, "\n") #83.25
cat("RMSE for MBA: ", rmse_mba, "\n") #20.72

# Plot the IDW and MBA interpolation surfaces side by side for comparison
par(mfrow = c(1, 2))  # Divide the plot area into 1 row, 2 columns
spplot(idw_result, main = "IDW Interpolation")  # Plot the IDW result
spplot(mba_sp, main = "MBA Interpolation")  # Plot the MBA result


# The lower RMSE value for the MBA (20.72) suggests that this method better models the Wolfcamp dataset, 
# capturing spatial variability more accurately and providing more precise estimates.
# The higher RMSE value for the IDW (83.25) suggests that this method has a lower predictive ability than the MBA 
# for this specific case.




