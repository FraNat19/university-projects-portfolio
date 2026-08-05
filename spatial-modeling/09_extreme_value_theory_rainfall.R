# Spatial Statistics and Statistical Tools for Environmental Data
# MSc in Statistical Methods and Applications, Sapienza University of Rome
# Group 12 — Agate, Bartocci, Natali
#
# Homework 9 — Extreme value theory on rainfall
Block-maxima approach and Generalized Extreme Value (GEV) modelling of yearly rainfall maxima.

# Import necessary libraries
library(dplyr)
library(extRemes)
library(evir)
library(ggplot2)
library(gridExtra)

# Load the dataset
load("datiVE.RData")

# GEV Model: 1. Define an appropriate time window to use in maxima building (maxima should be iid) ------------------------------------------------------------- 

# Combine rainfall data from the four different stations into one unified dataset;
# Following the block maxima approach, we extract yearly maxima to ensure 
# the data fits GEV assumptions of independence.
rainfall_data <- bind_rows(
  T4219 %>% mutate(Station_Code = "S4219", Rainfall = S4219) %>% select(date, Station_Code, Rainfall),
  T4220 %>% mutate(Station_Code = "S4220", Rainfall = S4220) %>% select(date, Station_Code, Rainfall),
  T4237 %>% mutate(Station_Code = "S4237", Rainfall = S4237) %>% select(date, Station_Code, Rainfall),
  T4238 %>% mutate(Station_Code = "S4238", Rainfall = S4238) %>% select(date, Station_Code, Rainfall)
) %>% mutate(Year = as.character(format(date, "%Y")))

# Handle missing data
rainfall_data <- rainfall_data %>% filter(!is.na(Rainfall))

# Calculate the yearly maximum rainfall values for each station
yearly_max_rainfall <- rainfall_data %>%
  group_by(Station_Code, Year) %>%
  summarise(Max_Rainfall = max(Rainfall), .groups = "drop")



# Ensuring independence and identically distributed (iid) of maxima is critical for valid 
# GEV model application.
# To verify independence:
# 1. Use autocorrelation plots to check for temporal correlations in yearly maxima.
# Plot the autocorrelation function for yearly maxima
acf(yearly_max_rainfall$Max_Rainfall, main = "Autocorrelation of Yearly Maxima")


# 2. Apply statistical tests (like the Ljung-Box test) to formally test independence.
# Perform Ljung-Box test to check for independence
ljung_box_test <- Box.test(yearly_max_rainfall$Max_Rainfall, lag = 10, type = "Ljung-Box")
print(ljung_box_test)

# Box-Ljung test
# 
# data:  yearly_max_rainfall$Max_Rainfall
# X-squared = 16.546, df = 10, p-value = 0.08503

# Therefore, based on the Ljung-Box test, we do not observe
# significant temporal dependence in the yearly maxima of rainfall data, meaning that the yearly maxima 
# can be considered independent for the purposes of the Generalized Extreme Value (GEV) model.



# 2. Estimate the GEV model for each rainfall station and comment on possible differences------------------------------------------------------

# Prepare storage for GEV model results and parameter data
gev_results_list <- list()
gev_params_df <- data.frame(Station = character(),
                            Location_Parameter = numeric(), Scale_Parameter = numeric(),
                            Shape_Parameter = numeric())

# Fit GEV models for each station and extract parameters
for (station in unique(yearly_max_rainfall$Station_Code)) {
  # Extract the rainfall data for the current station
  station_data <- filter(yearly_max_rainfall, Station_Code == station)$Max_Rainfall
  
  # Fit a GEV distribution to the station's data
  gev_model <- fevd(station_data, type = "GEV")
  gev_results_list[[station]] <- gev_model
  
  # Extract model parameters (location, scale, and shape)
  model_params <- gev_model$results$par
  gev_params_df <- rbind(gev_params_df,
                         data.frame(Station = station, Location_Parameter = model_params["location"],
                                    Scale_Parameter = model_params["scale"], Shape_Parameter = model_params["shape"]))
}


# Display the parameters of the GEV model for each station
row.names(gev_params_df) <- 1:nrow(gev_params_df)
print(gev_params_df)

#   Station      Location_Parameter  Scale_Parameter  Shape_Parameter
#1   S4219          104.37234             31.33693       0.19141874
#2   S4220           77.28865             27.90420       0.08023894
#3   S4237           66.89426             20.65030       0.71867937
#4   S4238          106.79034             39.62161       0.17561354

# We remind that the GEV distribution is a unique formalization for three distributions:
# the Gumbel, the Frèchet, and the Weibull.
# The shape parameter (csi) is overall the most relevant one, since it allows us to understand
# to which member of the family of distributions the single distribution belongs.
# All shape parameters are greater than -0.5, hence ML estimators are regular and they follow the common
# properties of MLE's (like consistency for instance).
# In particular, distributions of extreme values of each station could be associated to the Frèchet distribution 
# since csi>0 in every station, but
# the second station has the lowest shape parameter (0.08), hence, between all 4 distributions, the second is the
# nearest to the Gumbel distribution (where the value of the parameter tends to zero).
# In any case, since we have no negative values, none of the distributions can be a Weibull.
# We also remind that the Frèchet distribution is overall appropriate when the 
# distribution of the maxima has a heavy tail, while the Gumbel
# is more appropriate when the distribution of M_n has an exponential tail.
# The first and fourth distributions have their locations centered around 105 (mu parameter), 
# while the second and the third are centered around 70.
# Lastly, the fourth distribution has the highest scale parameter sigma, meaning that it has the highest variability, while
# the third station has the lowest scale.


# 3. Compute return levels of order 100, 200, 500 and 1000  ------------------------------------------------

# Define return periods and calculate return levels based on GEV
return_periods <- c(100, 200, 500, 1000)

# Calculate return levels for each station using the GEV model
gev_return_levels <- data.frame()  # Create an empty data frame to store the results

for (station in names(gev_results_list)) {
  # Compute return levels for the current station
  levels <- return.level(gev_results_list[[station]], return_periods)
  
  # Transform the results into a data frame and add the station name
  station_data <- data.frame(
    Station = station,
    Return_Period = return_periods,
    Return_Level = as.vector(levels)
  )
  
  # Append the data to the overall data frame
  gev_return_levels <- rbind(gev_return_levels, station_data)
}

# Print the results
print(gev_return_levels)

# Note: Higher return levels are sensitive to the shape parameter.
# A positive shape parameter indicates heavy-tailed distributions.

# [1] "Return Levels for period units in years" for Station S4219
# 100-year level  200-year level  500-year level 1000-year level 
#   335.5680        391.8167        478.4653        554.8308 

# [1] "Return Levels for period units in years" for Station S4220
# 100-year level  200-year level  500-year level 1000-year level 
#   232.5472        261.4240        302.0730        334.8431 

# [1] "Return Levels for period units in years" for Station S4237
# 100-year level  200-year level  500-year level 1000-year level 
#   821.9201       1330.3081       2537.1990       4152.2453 

# [1] "Return Levels for period units in years" for Station S4238
# 100-year level  200-year level  500-year level 1000-year level 
#   387.2511        453.0129        553.0240        640.0583 




# Now, to fit a GEV model for every station,
# we build a loop through each unique station in the rainfall data
for (station in unique(rainfall_data$Station_Code)) {
  # Extract the data of annual rainfall maxima for the current station
  station_data <- filter(yearly_max_rainfall, Station_Code == station)$Max_Rainfall
  
  # Fit the GEV model
  gev_model <- fevd(station_data, type = "GEV")
  
  # Create diagnostic plots with custom titles for each station
  par(mfrow = c(2, 2)) # Arrange plots in a 2x2 grid
  
  # Probability plot
  plot(gev_model, type = "prob", main = paste("Probability Plot - Station", station), cex.main = 1)
  
  # Quantile plot
  plot(gev_model, type = "qq", main = paste("Quantile Plot - Station", station), cex.main = 1)
  
  # Return level plot
  plot(gev_model, type = "rl", main = paste("Return Level Plot - Station", station), cex.main = 1)
  
  # Density plot
  plot(gev_model, type = "density", main = paste("Density Plot - Station", station), cex.main = 1)
}

# The purpose of the diagnostic is to see how well the model adapts to the data.
# The density plot shows the observed histogram and the GEV curve,
# in the probability plot we plot model estimated probabilities against empirical ones,
# and points should be on a line when the model is appropriate.
# Finally, in the return level plot we plot estimated return levels of order N
# against a range of N-values with their asymptotic confidence intervals.





#4. On the same data, estimate GPD model  -----------------------------------------------------------------

# Determine the 95th percentile threshold for each station:
# Following the Pickands-Balkema-de Haan theorem, thresholds must balance enough exceedances 
# while maintaining asymptotic behavior.
# Calculate the 95th percentile threshold and the number of exceedances for each station
station_thresholds <- rainfall_data %>%
  group_by(Station_Code) %>%
  summarise(
    Threshold_95 = quantile(Rainfall, 0.95, na.rm = TRUE),
    Exceedances = sum(Rainfall > quantile(Rainfall, 0.95, na.rm = TRUE))
  )

# Print the table with thresholds and exceedance counts for each station
print(station_thresholds)

# Merge the threshold data with the main dataset and filter for exceedances above the threshold
exceedance_data <- rainfall_data %>%
  left_join(station_thresholds, by = "Station_Code") %>%
  filter(Rainfall > Threshold_95) %>% # Retain only values exceeding the threshold
  mutate(Excess_Above_Threshold = Rainfall - Threshold_95) 

# Compute the average number of exceedances per year
year_range <- as.numeric(range(rainfall_data$Year))
total_years <- year_range[2] - year_range[1] + 1
avg_exceedances <- nrow(exceedance_data) / total_years
cat("Average number of exceedances per year:", avg_exceedances, "\n")  
# 72.95 
# 73 is overall a reasonable number for exceedances.


# Loop through each unique station in the exceedance data
for (station in unique(exceedance_data$Station_Code)) {
  # Extract the data of exceedances for the current station
  station_exceedances <- filter(exceedance_data, Station_Code == station)$Excess_Above_Threshold
  
  # Plot the ACF 
  acf(station_exceedances, main = paste("Autocorrelation of exceedances - Station", station), cex.main = 1)
}


# Peaks for other lags are low and do not exceed the threshold, indicating that there are 
# no significant correlations in the exceedances, which is a good sign of independence.

# Ljung-Box test for autocorrelation in the exceedances for each station
for (station in unique(exceedance_data$Station_Code)) {
  station_exceedances <- filter(exceedance_data, Station_Code == station)$Excess_Above_Threshold
  ljung_box_test <- Box.test(station_exceedances, lag = 10, type = "Ljung-Box")
  print(paste("Ljung-Box test for station", station))
  print(ljung_box_test)
}

# 
# [1] "Ljung-Box test for station S4219"
# 
# Box-Ljung test
# 
# data:  station_exceedances
# X-squared = 13.597, df = 10, p-value = 0.1922
# 
# [1] "Ljung-Box test for station S4220"
# 
# Box-Ljung test
# 
# data:  station_exceedances
# X-squared = 12.777, df = 10, p-value = 0.2364
# 
# [1] "Ljung-Box test for station S4237"
# 
# Box-Ljung test
# 
# data:  station_exceedances
# X-squared = 13.091, df = 10, p-value = 0.2186
# 
# [1] "Ljung-Box test for station S4238"
# 
# Box-Ljung test
# 
# data:  station_exceedances
# X-squared = 23.63, df = 10, p-value = 0.008647


# The p-values are high, indicating that there is no evidence of significant autocorrelation in these stations.
# In other words, we can assume the exceedances are independent.
# The fact that autocorrelation exceeds the threshold only at lags 0 and 1 suggests weak temporal dependence. 
# Often, in environmental or meteorological data, some short-term dependence is observed
# (e.g., correlation between consecutive exceedances).
# If this dependence is limited to immediate lags (0 and 1), 
# it might not overly compromise the independence of the exceedances in the long run.

# The autocorrelation plot shows significant peaks (high values) for higher lags, indicating potential dependence in the exceedances.


# Prepare storage for GPD model results and parameter data
gpd_results_list <- list()
gpd_params_df <- data.frame(Station_Code = character(),
                            Scale_Parameter = numeric(), Shape_Parameter = numeric(),
                            Threshold_Value = numeric())

# Now, we fit GPD models for each station and extract parameters
gpd_return_levels <- data.frame()

# Helper function to calculate return levels
get_return_levels <- function(fit_model, return_periods) {
  return.level(fit_model, return_periods)
}

for (station in unique(rainfall_data$Station_Code)) {
  station_data <- filter(rainfall_data, Station_Code == station)$Rainfall
  threshold_value <- as.numeric(quantile(station_data, 0.95, na.rm = TRUE))
  
  # Fit a GPD distribution to the exceedances above the threshold
  gpd_model <- fevd(x = station_data, threshold = threshold_value, type = "GP", 
                    period.basis = "year", time.units = "72.95/year")
  gpd_results_list[[station]] <- gpd_model
  
  # Extract model parameters (scale and shape) and the threshold
  model_params <- gpd_model$results$par
  gpd_params_df <- rbind(gpd_params_df, data.frame(Station_Code = station,
                                                   Scale_Parameter = model_params["scale"],
                                                   Shape_Parameter = model_params["shape"],
                                                   Threshold_Value = threshold_value))
  
  # Calculate return levels for the GPD model
  station_levels <- get_return_levels(gpd_model, return_periods)
  gpd_return_levels <- rbind(gpd_return_levels, 
                             data.frame(Station = station, Return_Period = return_periods, Return_Level = station_levels))
  row.names(gpd_return_levels) <- NULL
}


# In the GPD approach, the main issue is the choice of the threshold u.
# Indeed, it has to be large enough to justify the use of the GDP theory,
# and at the same time, small enough to allow a reasonably large amount of exceedances.

# The Mean Excess Plot is used to identify the optimal threshold 
# for the Generalized Pareto Distribution (GPD)
# Hence, we build a function to calculate the mean excess for a given threshold u using the corrected formula
mean_excess_plot <- function(station_data, u_values) {
  mean_excess_values <- sapply(u_values, function(u) {
    # Filter the exceedances above the threshold u
    exceedances <- station_data[station_data > u] - u
    
    # If there are any exceedances, compute the mean of those exceedances
    if (length(exceedances) > 0) {
      mean_excess <- sum(exceedances, na.rm = TRUE) / length(exceedances)
    } else {
      mean_excess <- NA  # Return NA if there are no exceedances above u
    }
    return(mean_excess)
  })
  return(mean_excess_values)
}

# Define a range of threshold values u (from the 90th to the 99th percentile)
u_values <- seq(from = quantile(rainfall_data$Rainfall, 0.90), 
                to = max(rainfall_data$Rainfall), 
                length.out = 50)

# Prepare an empty data frame to store the results
mean_excess_results <- data.frame()

# Apply the function for each station
for (station in unique(rainfall_data$Station_Code)) {
  # Filter the data for the station
  station_data <- filter(rainfall_data, Station_Code == station)$Rainfall
  
  # Calculate the mean excess for the station
  mean_excess_values <- mean_excess_plot(station_data, u_values)
  
  # Store the results for this station
  station_results <- data.frame(Station_Code = station, u = u_values, Mean_Excess = mean_excess_values)
  mean_excess_results <- rbind(mean_excess_results, station_results)  # Append results to the main data frame
}

# Plot the results with different lines for each station
ggplot(mean_excess_results, aes(x = u, y = Mean_Excess, color = Station_Code, linetype = Station_Code)) +
  geom_line(size = 1) +  # Draw lines with size 1
  labs(x = "Threshold Value (u)", y = "Mean Excess", title = "Mean Excess Plot per Station") +  # Add labels and title
  theme_minimal() +  # Use minimal theme
  theme(legend.title = element_blank())  # Remove legend title

# The Mean Excess Plot is used to identify the optimal threshold for the Generalized Pareto Distribution (GPD). 
# It shows how the mean excess changes as the threshold (u) increases. 
# On the y-axis we have the average of the exceedances, whereas 
# on the horizontal axis we find the thresholds.
# A linear relationship between mean excess and threshold suggests an appropriate threshold
# for modeling extreme values. 
# Significant deviations from linearity may indicate the optimal threshold
# for analyzing extreme rainfall events.


# Station S4219 shows a steady linear trend up to a threshold of around u=150, 
# after which it starts to fluctuate. This suggests that thresholds up to u=150
# could be appropriate for fitting a GPD.
# Station S4220: Exhibits a relatively linear trend but with a slight upward slope
# and variability starting at u=75. It becomes noisy as the threshold value increases.
# Station S4237: Displays a similar trend to Station S4220, with a stable 
# and slightly upward linear trend until u=125, after which variability increases.
# Station S4238: This station shows instead a different behavior with a steeper increase 
# in mean excess, especially after u=100, and becomes more irregular at higher thresholds.
# Overall, station S4238 exhibits the most variability at higher thresholds, 
# suggesting fewer data points for modeling extreme events.
# Station S4219 seems to have the most stable and linear trend over a range of thresholds, 
# making it a potentially strong candidate for GPD modeling.
# We can conclude that the stations generally exhibit a linear trend in the mean excess plot, 
# suggesting that the GPD might be an appropriate model for exceedances.


# Now we check if the shape parameter (csi) is less than 1 for GPD, which suggests appropriate behavior
gpd_shape_param <- gpd_model$results$par["shape"]
cat("Shape parameter:", gpd_shape_param, "\n")
# 0.1913807


# Display the parameters of the GPD model for each station
row.names(gpd_params_df) <- 1:nrow(gpd_params_df)
print(gpd_params_df)

#  Station_Code  Scale_Parameter      Shape_Parameter Threshold_Value
#1   S4219          18.03906              0.3193445         16.000
#2   S4220          18.45509              0.1574803         13.785
#3   S4237          13.80078              0.2253858         22.600
#4   S4238          14.42833              0.4352364         19.000



# 5. Compute return levels of order 100, 200, 500 and 1000-----------------------------------------------

# Calculate and display return levels for each station using the GPD model
for (station in names(gpd_results_list)) {
  cat("Return levels for Station", station, "using GPD model:\n")
  print(get_return_levels(gpd_results_list[[station]], return_periods))
}

# Note: GPD return levels can differ significantly for higher return periods
# due to sensitivity to the shape parameter.

# [1] "Return Levels for period units in years" for Station S4219
# 100-year level  200-year level  500-year level 1000-year level 
#   331.0795        423.1403        580.7383        734.6554 

# [1] "Return Levels for period units in years" for Station S4220
# 100-year level  200-year level  500-year level 1000-year level 
#   193.4222        227.6574        279.0476        323.1587

# [1] "Return Levels for period units in years" for Station S4237
# 100-year level  200-year level  500-year level 1000-year level 
#    192.7668        231.8939        293.9501        350.1863 

# [1] "Return Levels for period units in years" for Station S4238
# 100-year level  200-year level  500-year level 1000-year level 
#    416.7949        568.5448        854.0889       1159.8242 


# 6. Discuss differences between the two approaches -------------------------------------------------

# The main difference between the two models is that the GEV analyses 
# the distribution of maximum values through years, considering the whole sample of observed data,
# while GPD considers daily values that go above a fixed threshold.
# The main drawback of the GEV distribution approach is indeed the fact that
# we need at least maxima to be independent, i.e. we need to choose time windows that are quite large.
# GEV follows the Fisher-Tippett-Gnedenko theorem, for which, under particular assumptions, the unknown CDF
# F belongs either to the Gumbel, the Frèchet or the Weibull family. 
# On the other hand, GPD follows the Pickands-Balkema-de Haan theorem, for which the unknown CDF, 
# conditional to the threshold u, is approximated to a Generalized Pareto distribution.



# We now compare GEV and GPD return levels for Stations S4219 and S4238:
# Station S4219
gev_params_df[1, ]
#     Station Location_Parameter Scale_Parameter Shape_Parameter
# 1   S4219           104.3723        31.33693       0.1914187

gpd_params_df[ 1, ]
#     Station_Code Scale_Parameter Shape_Parameter Threshold_Value
# 1        S4219        18.03906       0.3193445              16


cat("Comparison of GEV and GPD return levels for Station S4219:\n")
cat("GEV return levels:\n")
print(get_return_levels(gev_results_list$S4219, return_periods))
# fevd(x = station_data, type = "GEV")
# get(paste("return.level.fevd.", newcl, sep = ""))(x = x, return.period = return.period)
# 
# GEV model fitted to  station_data  
# Data are assumed to be  stationary 
# [1] "Return Levels for period units in years"
# 100-year level  200-year level  500-year level 1000-year level 
# 335.5680        391.8167        478.4653        554.8308 
cat("GPD return levels:\n")
print(get_return_levels(gpd_results_list$S4219, return_periods))
# fevd(x = station_data, threshold = threshold_value, type = "GP", 
#      time.units = "72.95/year", period.basis = "year")
# get(paste("return.level.fevd.", newcl, sep = ""))(x = x, return.period = return.period)
# 
# GP model fitted to  station_data  
# Data are assumed to be  stationary 
# [1] "Return Levels for period units in years"
# 100-year level  200-year level  500-year level 1000-year level 
# 331.0795        423.1403        580.7383        734.6554 


# Station S4238
gev_params_df[4,]
#    Station Location_Parameter Scale_Parameter Shape_Parameter
#4   S4238           106.7903        39.62161       0.1756135

gpd_params_df[4,]
#     Station_Code Scale_Parameter Shape_Parameter Threshold_Value
# 4        S4238        14.42833       0.4352364              19

cat("Comparison of GEV and GPD return levels for Station S4238:\n")
cat("GEV return levels:\n")
print(get_return_levels(gev_results_list$S4238, return_periods))
# fevd(x = station_data, type = "GEV")
# get(paste("return.level.fevd.", newcl, sep = ""))(x = x, return.period = return.period)
# 
# GEV model fitted to  station_data  
# Data are assumed to be  stationary 
# [1] "Return Levels for period units in years"
# 100-year level  200-year level  500-year level 1000-year level 
# 387.2511        453.0129        553.0240        640.0583 
cat("GPD return levels:\n")
print(get_return_levels(gpd_results_list$S4238, return_periods))
# fevd(x = station_data, threshold = threshold_value, type = "GP", 
#      time.units = "72.95/year", period.basis = "year")
# get(paste("return.level.fevd.", newcl, sep = ""))(x = x, return.period = return.period)
# 
# GP model fitted to  station_data  
# Data are assumed to be  stationary 
# [1] "Return Levels for period units in years"
# 100-year level  200-year level  500-year level 1000-year level 
# 416.7949        568.5448        854.0889       1159.8242 


# In the case of the first and the last station, the return level of order 100 of both models is almost the same, but
# for higher orders the return level of GPD increases more rapidly than GEV; this could be caused by the higher value
# of the estimated shape parameter.



# Station S4220
gev_params_df[2,]
#   Station Location_Parameter Scale_Parameter Shape_Parameter
# 2   S4220           77.28865         27.9042      0.08023894

gpd_params_df[3,]
# Station_Code Scale_Parameter Shape_Parameter Threshold_Value
# 3   S4220        18.45509       0.1574803          13.785

cat("Comparison of GEV and GPD return levels for Station S4220:\n")
cat("GEV return levels:\n")
print(get_return_levels(gev_results_list$S4220, return_periods))
# fevd(x = station_data, type = "GEV")
# get(paste("return.level.fevd.", newcl, sep = ""))(x = x, return.period = return.period)
# 
# GEV model fitted to  station_data  
# Data are assumed to be  stationary 
# [1] "Return Levels for period units in years"
# 100-year level  200-year level  500-year level 1000-year level 
# 232.5472        261.4240        302.0730        334.8431
cat("GPD return levels:\n")
print(get_return_levels(gpd_results_list$S4220, return_periods))
# fevd(x = station_data, threshold = threshold_value, type = "GP", 
#      time.units = "72.95/year", period.basis = "year")
# get(paste("return.level.fevd.", newcl, sep = ""))(x = x, return.period = return.period)
# 
# GP model fitted to  station_data  
# Data are assumed to be  stationary 
# [1] "Return Levels for period units in years"
# 100-year level  200-year level  500-year level 1000-year level 
# 193.4222        227.6574        279.0476        323.1587 


# The second station is the best case between all four, because return levels of the two models behave similarly
# during return periods



# Station S4237
gev_params_df[3,]
# Station Location Scale     Shape
# 3 S4237 66.89426 20.6503 0.7186794

gpd_params_df[4,]
# Station Scale    Shape Threshold
# 3 S4237 13.80078 0.2253858 22.6

cat("Comparison of GEV and GPD return levels for Station S4237:\n")
cat("GEV return levels:\n")
print(get_return_levels(gev_results_list$S4237, return_periods))
# fevd(x = station_data, type = "GEV")
# get(paste("return.level.fevd.", newcl, sep = ""))(x = x, return.period = return.period)
# 
# GEV model fitted to  station_data  
# Data are assumed to be  stationary 
# [1] "Return Levels for period units in years"
# 100-year level  200-year level  500-year level 1000-year level 
# 821.9201       1330.3081       2537.1990       4152.2453 
cat("GPD return levels:\n")
print(get_return_levels(gpd_results_list$S4237, return_periods))
# fevd(x = station_data, threshold = threshold_value, type = "GP", 
#      time.units = "72.95/year", period.basis = "year")
# get(paste("return.level.fevd.", newcl, sep = ""))(x = x, return.period = return.period)
# 
# GP model fitted to  station_data  
# Data are assumed to be  stationary 
# [1] "Return Levels for period units in years"
# 100-year level  200-year level  500-year level 1000-year level 
# 192.7668        231.8939        293.9501        350.1863 


# The third station is the only one for which return levels are the most different, 
# meaning that the models can’t describe correctly the distribution of extreme values.



# Finally, we can combine GEV and GPD return levels for comparison
comparison_data <- rbind(
  data.frame(Station = gev_return_levels$Station, Model = "GEV",
             Return_Period = gev_return_levels$Return_Period, Return_Level = gev_return_levels$Return_Level),
  data.frame(Station = gpd_return_levels$Station, Model = "GPD",
             Return_Period = gpd_return_levels$Return_Period, Return_Level = gpd_return_levels$Return_Level)
)


# and plot return levels for each station

for (station in unique(comparison_data$Station)) {
  station_data <- filter(comparison_data, Station == station)
  print(  # Ensure the plot is printed in a loop
    ggplot(station_data, aes(x = Return_Period, y = Return_Level, color = Model)) +
      geom_line(linewidth = 1.2) +  # Updated to use `linewidth`
      geom_point(size = 2) +
      scale_x_log10() +  # Logarithmic scale for return periods
      labs(title = paste("Comparison of GEV and GPD Return Levels for Station", station),
           x = "Return Period (years)",
           y = "Return Level (mm)",
           color = "Model") +
      theme_minimal()
  )
}


# These plots confirm what we previously obtained, that is, the first and fourth station 
# show a consistent increasing trend for GPD return levels (which are consistently higher than GEV),
# whereas stations S4220 and S4237 show higher GEV return levels overall.







plots <- list()

# Loop per creare i grafici e salvarli nella lista
for (station in unique(comparison_data$Station)) {
  station_data <- filter(comparison_data, Station == station)
  plot <- ggplot(station_data, aes(x = Return_Period, y = Return_Level, color = Model)) +
    geom_line(linewidth = 1.2) +
    geom_point(size = 2) +
    scale_x_log10() +
    labs(
      title = paste("Comparison of GEV and GPD Return Levels for Station", station),
      x = "Return Period (years)",
      y = "Return Level (mm)",
      color = "Model"
    ) +
    theme_minimal()
  
  # Aggiungi il grafico alla lista
  plots[[station]] <- plot
}

# Mostra i grafici in una griglia 2x2
grid.arrange(grobs = plots, ncol = 2, nrow = 2)

