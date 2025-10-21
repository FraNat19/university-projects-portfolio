# HOMEWORK 3 - Group 12 (Agate, Bartocci, Natali)

# Using wolfcamp data, prepare a script developing the following tasks:


# PART 1 ------------------------------------------------------------------

# 1) Choose the variogram model (justify your choice).

#We first load the dataset and the required package
library(geoR)
data("wolfcamp")

class(wolfcamp) #geodata
summary(wolfcamp)  #we get a summary of the coordinates, 
# of the min. and max. distances among the points, and then a data summary

plot(wolfcamp)
#The top left plot, the plot of values of coordinates, shows a spatial trend, so we know for sure that 
#the mean is not constant


#We want to understand the spatial variation of this data
#, so we now obtain a plot of the empirical variogram, with dotsw and lines:
plot(variog(wolfcamp), type="b") 
#We can notice that there is no nugget effect, but we don't see an asymptote either, 
#hence we are dealing with a non-stationary phenomenon.


#We can now start modelling the mean, since we know it is not constant, by choosing a trend
plot(variog(wolfcamp, trend="1st"),type="b")
#we obtain something that increases with a distance from the origin, then becomes flat, and then jumps down.
#This happens because when distances are large, we have less pairs of differences, 
#and after a certain distance
#we might have a less robust estimation of this function.

#Moreover, points very far away from each other might have a negligible correlation,
#so we can now add the option "max distance" to the function.
#In particular, we now estimate the empirical variogram up to distance 250, in  order to cut out the erratic part
vv <- variog(wolfcamp, trend="1st",max.dist = 250)
plot(vv,type="b")
#Now this looks like a stationary variogram.
#Indeed, with the "trend 1st" option we have been fitting a linear trend in the coordinates,
#which removes the non-stationary behaviour we previously had (i.e. we have now obtained a stationary process).
#Furthermore, we now have a clear nugget effect, since with reparametrization increase
#if, for instance, we use as trend the second order polinomyal,
#we increase the nugget effect even more, until we might get a constant variogram;
#this is especially true when our second term in  the model is a random term (i.e. not a fixed effect).
#This means that we have already overparametrized the model, that is, it has likely 
#"absorbed" all the information in the data, something which would not enable us
#to run our spatial interpolation.


#To estimate the covariates we can now use eyefit
eyefit(vv) 
#graphical estimation (try to fit a theoretical variogram)

#Before checking the shape of several different variogram models, 
#we can already exclude some options through some a priori considerations:
#a) Since we know that we do not have a periodical behaviour in space, the wave model is not our best choice;
#b) since the variogram does not start with a straight line, we can claim
# that also the Gaussian and the cubic models do not fit at all the shape of our empirical variogram;
#c) We have now obtained a stationary empirical variogram, therefore we can exclude 
#non-stationary models, such as the linear or the power


#After all these considerations, we now first give a look at the exponential and spherical models, 
#also with graphical representations,
#by estimating with the likelihood function
vv.exp<-likfit(wolfcamp,trend = "1st",cov.model="exponential", ini.cov.pars = c(2000,100))
lines(vv.exp,col=2) #to see how well it fits with our empirical estimation
summary(vv.exp)

vv.sph<-likfit(wolfcamp,trend = "1st",cov.model="spherical", ini.cov.pars = c(2000,100))
lines(vv.sph,col=3)
summary(vv.sph)

#The exponential is reaching the sill asymptotically,
#wheres for the spherical model the sill is reached for a finite (i.e. one) value of the range.



#Hence, we can claim that we may get even better results with a more careful choice of the covariance model.
#Through eyefit, we see whether the Mat?rn could have an even better fit; these are the parameter obtained by using eyefit

#   cov.model sigmasq    phi   tausq kappa kappa2   practicalRange
# 1    matern  2884.5 110.44 1562.44   0.3   <NA> 265.833503445485


vv.mat <- likfit(wolfcamp,trend="1st", cov.model = "matern", ini.cov.pars = c(2884.5, 110.44), kappa=0.3)
lines(vv.mat,col=4)

#In the end, after looking at the various possible models,
#the best choice overall seems indeed to be the Mat?rn model,
#so we can now choose it as our variogram model

#We finally give a summary of the likelihood estimates
summary(vv.mat)


# PART 2 ------------------------------------------------------------------

# 2) Run a spatial interpolation using kriging in the likelihood framework.

# Kriging with the "geoR" package under MLE can be considered a 2-step procedure:
# Step 1 involves the choice of the variogram model, which we have already made in Part 1 above.

# Step 2 instead means to run kriging with fixed covariance 


#Now we need to run interpolation
#We need a grid in which run the estimation, so we create our interpolation grid:
x <- seq(min(wolfcamp$coords[,1]), max(wolfcamp$coords[,1]), length = 20)
# 20 between the min. observed and the max. observed
y <- seq(min(wolfcamp$coords[,2]), max(wolfcamp$coords[,2]), length = 20)

grid <- expand.grid(x,y)

plot(grid)
# We obtain a 20 x 20 grid (400 points)


# Now we check how to implement things in the kriging function
# We first control for the kriging:
# We have to build a krige.control list telling to the krig.conv
# all the elements that are required.
kk.c <- krige.control(trend.d = "1st", trend.l = "1st", obj.model = vv.mat)

# For the interpolation procedure, we have to exactly define the model on the new locations,
# to specify the trend on both the model and the new location.
# We need to set the locations where interpolation happens, i.e. they are just cooordinates.
# Finally, we need to connect to the control
krig1 <- krige.conv(wolfcamp, locations = grid, krige = kk.c)

# now we can picture this
image(krig1)
# We can also add
contour(krig1,add=T,nlev=20)
# and see what we have in
names(krig1)

# We now define the kriging prediction variance
krige_var_matrix <- matrix(krig1$krige.var, nrow = length(x), ncol = length(y))

# and then we plot it
image(x, y, krige_var_matrix, main = "Kriging Prediction Variance")
contour(x, y, krige_var_matrix, add = TRUE, nlevels = 10)


#We can also build confidence intervals at the 95% level
low <- krig1$predict - 1.96*sqrt(krig1$krige.var)
up <- krig1$predict + 1.96*sqrt(krig1$krige.var)

#We can at first build one map for the lower bound and one map for the other, and then we can compare both

hist(low)
hist(up)

#We can also build a map for the length of these distances (i.e. for the differences)

hist(up-low)

# We now reshape lower and upper bounds into matrices
lower_bound_matrix <- matrix(low, nrow = length(x), ncol = length(y))
upper_bound_matrix <- matrix(up, nrow = length(x), ncol = length(y))

# Plot the lower bound
image(x, y, lower_bound_matrix, main = "95% Confidence Interval Lower Bound")
contour(x, y, lower_bound_matrix, add = TRUE, nlevels = 20)

# Plot the upper bound
image(x, y, upper_bound_matrix, main = "95% Confidence Interval Upper Bound")
contour(x, y, upper_bound_matrix, add = TRUE, nlevels = 20)

# We print the values for the first 10 locations in both bounds
print(low[1:10]) 
print(up[1:10])


# We also reshape the predicted values into a matrix for plotting
predict_matrix <- matrix(krig1$predict, nrow = length(x), ncol = length(y))

# Plot the predicted values
image(x, y, predict_matrix, main = "Kriging Predictions with 95% Confidence Interval")
contour(x, y, predict_matrix, add = TRUE, nlevels = 20)

# We add lower confidence bound as a contour
contour(x, y, lower_bound_matrix, add = TRUE, nlevels = 20, col = "blue", lwd = 2, lty = 2)

# and then add the upper confidence bound as a contour
contour(x, y, upper_bound_matrix, add = TRUE, nlevels = 20, col = "red", lwd = 2, lty = 2)

# Finally, we add a legend to distinguish the prediction and the confidence interval bounds
legend("topright", legend = c("Predictions", "Lower Bound", "Upper Bound"), 
       col = c("black", "blue", "red"), lty = c(1, 2, 2), lwd = c(1, 2, 2))





