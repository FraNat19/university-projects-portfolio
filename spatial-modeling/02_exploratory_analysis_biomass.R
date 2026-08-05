# Spatial Statistics and Statistical Tools for Environmental Data
# MSc in Statistical Methods and Applications, Sapienza University of Rome
# Group 12 — Agate, Bartocci, Natali
#
# Homework 2 — Exploratory analysis of shrimp biomass and environmental descriptors
EDA, correlations and multivariate analysis on the shrimp data (years 2002 and 2008).

# HOMEWORK 2 - GROUP 12 (Agate, Bartocci, Natali)

# PART A ------------------------------------------------------------------

# 1) Run an exploratory data analysis and comment on the results. 
# In particular, analyze the relationship between the total biomass (variable tot) 
# and the environmental descriptors.

# First we load the necessary libraries
library(ggplot2)
library(dplyr)
library(reshape2)
library(corrplot)
library(ade4)
library(sp)

# We also load the dataset and the file containing the estimation grids
load("shrimpsfull.RData")
load("AllGrids.RData")

# We now filter data for our years of interest, i.e. 2002 and 2008
shrimp_data_2002 <- shrimpsdata[shrimpsdata$ANNO == 2002, ]
shrimp_data_2008 <- shrimpsdata[shrimpsdata$ANNO == 2008, ]

colSums(is.na(shrimp_data_2002))
colSums(is.na(shrimp_data_2008))

# Display the structure of the dataset
str(shrimpsdata)
# Preview the first few rows
head(shrimpsdata)


# We now start our analysis by displaying descriptive Statistics for all the variables in both years
summary(shrimp_data_2002)
summary(shrimp_data_2008)

#In the summary statistics, Salinity averages are relatively stable between years (all values are around 38 mg/l).
#However, we can notice that there is a slight decrease in minimum salinity in 2008.
#The minimum and maximum temperatures also show some stability, with slight variations between the two years.
#We might assume that the increase in average temperature could have implications for shrimp populations and their total biomass.
#The mean depth (bathymetry) became less deep in 2008, which could also affect the total biomass,
#while both the average distance to the coast and the slope of the seabed remained substantially unchanged.


#For what regards the total biomass instead, in 2002 it has great variability,
#with the mean significantly higher than the median, 
#suggesting the presence of some extremely high values influencing the mean.
#Again in 2008, the mean is much higher than the median, indicating a similar distribution to that of 2002.
#Moreover, there is a significant increase in the average total biomass from 2002 to 2008, thus indicating
#that 2008 might have some outliers with even more extreme values on the right tail of the distribution.
#We hence plot an histogram for both years to verify our assumption


hist(shrimp_data_2002$tot, main="Distribution of Total Biomass in 2002", xlab="Total Biomass") 
#higher frequency of values <200 biomass gives an explanation for the lower mean
hist(shrimp_data_2008$tot, main="Distribution of Total Biomass in 2008", xlab="Total Biomass")

#As we can see in the histograms (with both distributions being right skewed),
#the values of <200 biomass have the highest frequency in both years, 
#but 2008 presents as expected more variability due to extreme outliers,
# which consequently increase the mean.




#We now regroup all the temperature and salinity quarters into single objects
temp_vars <- c("temp.minq3p", "temp.minq4p", "temp.minq1", "temp.minq2", "temp.minq3", 
               "temp.maxq3p", "temp.maxq4p", "temp.maxq1", "temp.maxq2", "temp.maxq3")

salinity_vars <- c("salinity.minq3p", "salinity.minq4p", "salinity.minq1", "salinity.minq2", "salinity.minq3", 
                   "salinity.maxq3p", "salinity.maxq4p", "salinity.maxq1", "salinity.maxq2", "salinity.maxq3")


#and we now choose the following variables of interest to analyze the correlation between the total biomass (tot)
#and the chosen covariates
vars_of_interest <- c("tot", "bat", "dist", "slope", 
                      "salinity.minq3p", "salinity.minq4p", "salinity.minq1", "salinity.minq2", "salinity.minq3", 
                      "salinity.maxq3p", "salinity.maxq4p", "salinity.maxq1", "salinity.maxq2", "salinity.maxq3",
                      "temp.minq3p", "temp.minq4p", "temp.minq1", "temp.minq2", "temp.minq3",
                      "temp.maxq4p", "temp.maxq1", "temp.maxq2", "temp.maxq3")

#Remove possible missing values for the selected variables
shrimp_data_2002_clean <- shrimp_data_2002 %>% select(all_of(vars_of_interest)) %>% na.omit()
shrimp_data_2008_clean <- shrimp_data_2008 %>% select(all_of(vars_of_interest)) %>% na.omit()


#we can now obtain and plot the correlation matrix between the total biomass and the environmental indicators;
#first we plot the matrix for 2002 data
cor_matrix_2002 <- cor(shrimp_data_2002_clean)
cor_matrix_2002
corrplot(cor_matrix_2002, method = "color", type = "lower", tl.cex = 0.7, main = "Correlation Matrix - 2002")

#Among environmental factors, temperature (in particular temp.minq2, temp.minq3, temp.maxq3) 
#never shows negative correlations in all quarters overall, and in particular it shows mild positive correlations 
#with biomass (around 0.21-0.23), indicating some moderate influence.
#Salinity also has positive correlations with tot, but they are much weaker compared to temperature, 
#mostly below 0.1, thus indicating that it might not play a strong role in driving biomass changes.
#Overall, stronger positive correlations can be found between salinity and temperature variables,
#since biologically speaking higher temperature usually implies higher salinity.
#The distance from the coast and the slope moreover show an almost negligible negative correlation, 
#while the depth depicts a weak correlation with positive sign


#We can now plot a linear regression for each of the last three variables

ggplot(shrimp_data_2002, aes(x = shrimp_data_2002$bat, y = shrimp_data_2002$tot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "green") +
  labs(title = "Total Biomass vs Bathymetry (2002)", x = "Bathymetry (m)", y = "Total Biomass")

ggplot(shrimp_data_2002, aes(x = shrimp_data_2002$dist, y = shrimp_data_2002$tot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  labs(title = "Total Biomass vs Distance from Coast (2002)", x = "Distance from Coast (m)", y = "Total Biomass")

ggplot(shrimp_data_2002, aes(x = shrimp_data_2002$slope, y = shrimp_data_2002$tot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "orange") +
  labs(title = "Total Biomass vs Slope (2002)", x = "Slope (degrees)", y = "Total Biomass")

#The regression lines respectively confirm that bathymetry has a correlation with positive sign 
#(i.e. total biomass decreases with increasing depth), something which contradicts the usual behaviour of such type of shrimps,
#whereas dist has an almost flat relationship and the slope a slightly negative correlation (that is, 
# total biomass decreases on average with steeper seabeds)


# Combine temperature and salinity variables in one dataframe
shrimp_data_long_2002 <- melt(shrimp_data_2002, id.vars = "tot", 
                         measure.vars = c(temp_vars, salinity_vars),
                         variable.name = "variable",
                         value.name = "value")

# Creating facets to compare each temperature and salinity variable with the total biomass
ggplot(shrimp_data_long_2002, aes(x = value, y = tot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "green") +
  facet_wrap(~variable, scales = "free_x") +
  labs(title = "Total Biomass vs Temperature and Salinity Variables (2002)", 
       x = "Environmental Variable Value", y = "Total Biomass")

#In almost all facets, the relationship between total biomass and environmental variables appears to be weak. 
#The regression lines (in green) tend to be rather flat, indicating that there is no strong or obvious 
#relationship between these variables and biomass.
#For temperature, some slight increases are observed in facets such as temp.minq3 and temp.minq2,
#but the correlation seems to be very weak or almost absent in many other facets.
#This suggests that temperature does not have a clearly identifiable impact on total biomass in this dataset, 
#at least not in a linear or direct way. The same applies to salinity, where only a few facets, 
#such as salinity.maxq3p and salinity.minq2, show slight variations, but the relationship is still rather weak.







# We now plot the correlation matrix for the year 2008
cor_matrix_2008 <- cor(shrimp_data_2008_clean)
cor_matrix_2008
corrplot(cor_matrix_2008, method = "color", type = "lower", tl.cex = 0.7, main = "Correlation Matrix - 2008")

#The correlation between total biomass (tot) and bathymetry (bat) though having a negative sign 
#appears to be positive, contrary to the 2002 data,
#suggesting that shrimp biomass tends to increase in deeper areas. 
#This makes sense if deeper areas provide more suitable habitats for shrimps, 
#possibly due to temperature, pressure, or predator/prey dynamics.
#Indeed, according to the biology of the "gambero rosa" (i.e.parapeneus longirostris), it actually confirms
#that when the shrimps' biomass increases, it looks like they move towards deeper habitats (mostly to around 450 m)

#The relationship between biomass and distance from the coast (dist) also shows a weak positive correlation
#(light blue color) and also in this case we have a different behaviour compared to 2002. 
#This indicates a slight tendency for biomass to increase as shrimps move further 
#from the coast, though the relationship is not very strong.
#Very similar values are shown also in the correlation between biomass and the seabed's slope,
#which once again depicts an opposite behaviour to 2002.

#On the other hand, salinity and temperature show moderately stronger correlations with total biomass,
#but their signs show some variability: contrary to the majority of quarters, one quarter for temperature (temp.maxq2)
#shows negative correlation, whereas for salinity they are divided equally between negative (reddish color)
#and positive correlations (light blue color) thus leaving an ambiguous effect overall of such covariate on total biomass.


#We now plot the regression lines to confirm our analysis in a graphical way

ggplot(shrimp_data_2008, aes(x = shrimp_data_2008$bat, y = shrimp_data_2008$tot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "green") +
  labs(title = "Total Biomass vs Bathymetry (2008)", x = "Bathymetry (m)", y = "Total Biomass")

ggplot(shrimp_data_2008, aes(x = shrimp_data_2008$dist, y = shrimp_data_2008$tot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "purple") +
  labs(title = "Total Biomass vs Distance from Coast (2008)", x = "Distance from Coast (m)", y = "Total Biomass")

ggplot(shrimp_data_2008, aes(x = shrimp_data_2008$slope, y = shrimp_data_2008$tot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "orange") +
  labs(title = "Total Biomass vs Slope (2008)", x = "Slope (degrees)", y = "Total Biomass")


# Combine temperature and salinity variables in one dataframe
shrimp_data_long_2008 <- melt(shrimp_data_2008, id.vars = "tot", 
                         measure.vars = c(temp_vars, salinity_vars),
                         variable.name = "variable",
                         value.name = "value")

# Creating facets to compare each temperature and salinity variable with the total biomass
ggplot(shrimp_data_long_2008, aes(x = value, y = tot)) +
  geom_point() +
  geom_smooth(method = "lm", se = FALSE, color = "green") +
  facet_wrap(~variable, scales = "free_x") +
  labs(title = "Total Biomass vs Temperature and Salinity Variables (2008)", 
       x = "Environmental Variable Value", y = "Total Biomass")

#The analysis shows that there are no very strong correlations between total biomass, temperature and salinity. 
#Some temperature variables, such as ‘temp.minq1’ and ‘temp.maxq4p’, indicate a slight increase in biomass 
#with increasing temperature, while salinity seems to have less of an impact. 
#The scattering of data suggests that other factors influence biomass, necessitating further analysis 
#with more complex models or by including additional variables.




# PART B ------------------------------------------------------------------


# 2) Using cross-validation, choose the optimal value for p in the IDW estimator.

#We first ensure the package "gstat" is loaded for IDW
library(gstat)

# We build an improved function to choose the best p for IDW
choose_p_idw_cv <- function(observed_data, value_col, p_range = seq(0.5, 4, 0.25), 
                            k = 5, seed = 1234) {
  
  set.seed(seed)
  n <- nrow(observed_data) # Number of observations
  
  # Create random fold assignments for k-fold cross-validation
  fold_indices <- sample(rep(1:k, length.out = n)) 
  

  errors <- numeric(length(p_range))
  
  #We now create the loop
  for (i in seq_along(p_range)) {
    p <- p_range[i]
    fold_rmse <- numeric(k)  
    
    for (fold in 1:k) {
      # Split data into training and testing sets
      test_indices <- which(fold_indices == fold)
      train_indices <- setdiff(1:n, test_indices)
      
      training_data <- observed_data[train_indices, ]
      test_data <- observed_data[test_indices, ]
      
      # Perform IDW interpolation on the test set using the training set
      idw_result <- idw(formula = as.formula(paste(value_col, "~ 1")),
                        locations = training_data,
                        newdata = test_data,
                        idp = p)
      
      
      # Calculate residuals and RMSE for this fold
      observed_values <- test_data[[value_col]]
      predicted_values <- idw_result$var1.pred
      residuals <- observed_values - predicted_values
      
      fold_rmse[fold] <- sqrt(mean(residuals^2))  # Root mean square error for this fold
    }
    
    
    errors[i] <- mean(fold_rmse)
  }
  
  # Find the best p value (the one with the lowest RMSE)
  optimal_p <- p_range[which.min(errors)]
  
  return(list(optimal_p = optimal_p, rmse = min(errors)))
}

p_values <- seq(0.5, 4, by = 0.02)


#We remind that the choice of the smoothing parameter p has to be done year by year
#since IDW, as a spatial interpolator, does not let us add time before using one year at a time

#Hence, we start with 2002

# Prepare observed_data as a spatial object for 2002:
coords <- cbind(shrimp_data_2002$X, shrimp_data_2002$Y)  # Extract coordinates (x, y)
shrimp_data_2002_sp <- SpatialPointsDataFrame(coords = coords, 
                                              data = data.frame(z = shrimp_data_2002$tot), 
                                              proj4string = CRS("+proj=utm +zone=33 +datum=WGS84")) # Create spatial object with CRS

# Prepare grid using the provided grid_2002
grid_2002_sp <- grid_2002
coordinates(grid_2002_sp) <- ~X + Y
gridded(grid_2002_sp) <- TRUE

# Make sure the grid has the same CRS as shrimp_data_2002_sp:
proj4string(grid_2002_sp) <- proj4string(shrimp_data_2002_sp)


# Apply the IDW interpolation using the best 'p' value obtained from cross-validation
best_p_2002 <- choose_p_idw_cv(shrimp_data_2002_sp, value_col = "z", p_range = p_values, k = 5)$optimal_p 
# Get the best 'p' value by running the function
best_p_2002




#2008 analysis

# Prepare observed_data as a spatial object for 2008:
coords <- cbind(shrimp_data_2008$X, shrimp_data_2008$Y)  # Extract coordinates (x, y)
shrimp_data_2008_sp <- SpatialPointsDataFrame(coords = coords, 
                                              data = data.frame(z = shrimp_data_2008$tot), 
                                              proj4string = CRS("+proj=utm +zone=33 +datum=WGS84")) # Create spatial object with CRS

# Prepare grid using the provided grid_2008
grid_2008_sp <- grid_2008
coordinates(grid_2008_sp) <- ~X + Y
gridded(grid_2008_sp) <- TRUE

# Assign the same CRS to the grid as shrimp_data_2008_sp to ensure consistency.
proj4string(grid_2008_sp) <- proj4string(shrimp_data_2008_sp)


# Apply the IDW interpolation using the best 'p' value obtained from cross-validation
best_p_2008 <- choose_p_idw_cv(shrimp_data_2008_sp, value_col = "z", p_range = p_values, k = 5)$optimal_p  # Get the best 'p' value by running the function
best_p_2008 



# PART C ------------------------------------------------------------------


# 3) Using the IDW estimator and the optimal p, map the total biomass values on the grid 
# stored in the object grid.


#Again, we proceed year by year and we start from 2002

# Assuming `idw_result` already contains the interpolated values
idw_result_2002 <- idw(z ~ 1, locations = shrimp_data_2002_sp, newdata = grid_2002_sp, idp = best_p_2002)


# Extract the predicted values and add to grid_2002_sp
grid_2002_sp$predicted_tot <- idw_result_2002$var1.pred

# Convert the grid to a data frame for easier plotting with ggplot2.
grid_2002_df <- as.data.frame(grid_2002_sp)
grid_2002_df <- grid_2002_df[, c("X", "Y", "predicted_tot")]

# We now plot the interpolated map for 2002
ggplot(grid_2002_df, aes(x = X, y = Y, fill = predicted_tot)) +
  geom_tile() +
  scale_fill_viridis_c() +  # Optional, for a nice color scale
  labs(title = "IDW Interpolated Total Biomass (2002)",
       x = "X Coordinate", y = "Y Coordinate", fill = "Biomass") +
  theme_minimal()



#2008

idw_result_2008 <- idw(z ~ 1, locations = shrimp_data_2008_sp, newdata = grid_2008_sp, idp = best_p_2008)

# Extract the predicted values and add to grid_2008_sp
grid_2008_sp$predicted_tot <- idw_result_2008$var1.pred

# Convert the grid_2008_sp to a data frame for easier plotting.
grid_2008_df <- as.data.frame(grid_2008_sp)
grid_2008_df <- grid_2008_df[, c("X", "Y", "predicted_tot")]

# Plot the interpolated map for 2008.
ggplot(grid_2008_df, aes(x = X, y = Y, fill = predicted_tot)) +
  geom_tile() +
  scale_fill_viridis_c() +  # Optional, for a nice color scale
  labs(title = "IDW Interpolated Total Biomass (2008)",
       x = "X Coordinate", y = "Y Coordinate", fill = "Biomass") +
  theme_minimal()


#In both 2002 and 2008, clusters of higher biomass are aligned along the diagonal representing the Italian coastline (from Genoa to Gaeta). 
#However, in 2008, the high biomass region shifts slightly southeast, suggesting a change in shrimp distribution. 
#Biomass values are higher in 2008, with peaks over 300 compared to 150 in 2002, indicating a possible increase in total biomass.
#Environmental factors likely played a role in these changes: maximum salinity remained stable in years, 
#though a slight decrease in the minimum value could have influenced shrimp distribution.
#Temperature increased slightly through years, which may explain the rise in biomass in certain areas, 
#as warmer waters can affect shrimp spawning and recruitment.
#Meanwhile bathymetry (depth) became shallower, potentially creating more suitable habitats closer to the coast.
#Overall, the spatial distribution of high biomass is more concentrated in 2008, 
#likely reflecting the combined effects of temperature and depth changes on shrimp populations.







