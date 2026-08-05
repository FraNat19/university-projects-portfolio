# Spatial Statistics: Kriging Models for Mediterranean Shrimp Biomass

![R](https://img.shields.io/badge/R-geoR%20%7C%20gstat%20%7C%20spdep-blue)
![Course](https://img.shields.io/badge/Sapienza-Spatial%20Statistics%202024%2F25-8b1a1a)

Semester project for **Spatial Statistics and Statistical Tools for Environmental Data**, MSc in Statistical Methods and Applications, Sapienza University of Rome.

The central work predicts deep-water rose shrimp (*Parapenaeus longirostris*) biomass across the Italian Tyrrhenian coast from 120 trawl survey stations, comparing maximum-likelihood kriging against Bayesian kriging. Four further assignments apply the same geostatistical toolkit to other datasets.

`Spatial_Presentation.pdf` is the final deliverable and carries the maps, priors, and comparison tables discussed below.

## What is in this folder

The eight scripts are weekly assignments. Four use the MEDITS shrimp data, four use other datasets that the course provided for specific techniques.

| Script | Dataset | What it does |
|---|---|---|
| `01_idw_interpolation.R` | Wolfcamp aquifer | Inverse distance weighting; picks the exponent *p* by spatial cross-validation |
| `02_exploratory_analysis_biomass.R` | MEDITS shrimp | EDA, correlation structure, PCA to select covariates |
| `03_variogram_modelling.R` | Wolfcamp aquifer | Trend removal, stationarity checks, empirical and fitted variograms |
| `04_kriging_biomass_estimation.R` | MEDITS shrimp | ML kriging with covariates; prediction maps and 95% confidence bounds |
| `05_bayesian_kriging.R` | MEDITS shrimp | Bayesian kriging via `krige.bayes`; credible intervals; comparison with script 04 |
| `07_spatial_econometrics_sids.R` | SIDS, North Carolina | Areal data: spatial weights, Moran's I, spatial regression |
| `08_spatial_species_distribution.R` | Mite data | Presence-absence modelling with spatial predictors, evaluated by ROC/AUC |
| `09_extreme_value_theory_rainfall.R` | Yearly rainfall | Block maxima and GEV fitting |

`shrimp2002.csv` and `shrimp2008.csv` hold the 120 stations for each survey year. `ML.RData` and `Bayesian.rdata` hold the fitted kriging objects, and the four `*_Kriging_*.pdf` files are the prediction maps they produce.

## The shrimp problem

MEDITS (Mediterranean International Trawl Survey) runs standardised trawls to monitor demersal species. Each survey year gives 120 sampled stations along the coast from Liguria to Lazio, with 29 variables per station: bathymetry, distance from coast, slope, and quarterly minima and maxima for temperature and salinity.

The response is total biomass in kg/km². Two things make it awkward:

- **Exactly half the stations record zero biomass**, 60 of 120 in both 2002 and 2008. The distribution is far from Gaussian.
- **Spawners and recruits sum to total biomass**, so neither can serve as a covariate.

The response is modelled as `log(tot + 1)`, which keeps the zeros at zero and leaves the non-zero values roughly symmetric. Covariates are selected from the loadings of the first two principal components: salinity, bathymetry, distance from coast, slope, and quarterly temperature, with the specific quarters differing between 2002 and 2008.

## Model

The spatial process is a mixed-effect model:

$$Z(s) = X(s)\beta + W(s) + \varepsilon(s)$$

with $W(s) \sim N_n(0, \sigma^2 H_{11}(\phi))$ the spatially correlated effect and $\varepsilon(s) \sim N_n(0, \tau^2)$ the nugget.

The ML route fits the variogram with `likfit` and predicts with `krige.conv`. The Bayesian route puts a flat normal prior on $\beta$, a reciprocal prior on $\sigma^2$, and a uniform prior on $\phi$ over a discretised range, then predicts with `krige.bayes`. Both are validated by repeated train-test splits.

## Results

Averages over 10 cross-validation repetitions:

| | ML kriging | Bayesian kriging |
|---|---|---|
| RMSE 2002 | 2.2617 | 2.0759 |
| RMSE 2008 | 2.3812 | 2.1481 |
| Interval score 2002 | 9.731 | 8.384 |
| Interval score 2008 | 9.708 | 9.95 |

Bayesian kriging gives lower RMSE in both years, and a clearly better interval score in 2002. In 2008 the interval scores are effectively tied, with the ML approach marginally ahead. The gap is small enough that the two methods are close to interchangeable on this dataset, and the Bayesian fit costs considerably more computation. What the priors bought here was a modest gain in point accuracy, not a decisive one.

Spatially, both approaches agree on where the shrimp are. Biomass concentrates along the Tuscany and Lazio coasts at 200-400 m depth, and stays low off Liguria, where the bathymetric gradient is steep. Between 2002 and 2008 the hotspots shift south-east, and 2008 shows a wider spread of high values. The two methods disagree most in 2008: ML spreads the peak toward the Lazio coast while the Bayesian fit keeps it concentrated in the central area.

## Running it

Written for R with `geoR`, `gstat`, `sp`, `sf`, `spdep`, `spatialreg`, `vegan`, `extRemes`, and `ggplot2`.

Scripts 02, 04, and 05 open `shrimpsfull.RData` (the full multi-year survey panel) and `AllGrids.RData` (the prediction grids), and script 09 opens `datiVE.RData`. Those three files are course data and are not redistributed here, so the shrimp scripts will not run end to end from a clone. The two CSVs cover the 2002 and 2008 stations, and `ML.RData` and `Bayesian.rdata` contain the fitted objects, which is enough to inspect the models and reproduce the maps without refitting. Scripts 01, 03, 07, and 08 run as they are: their datasets ship with `geoR`, `spData`, and `vegan`.

## References

- MEDITS Programme: [sibm.it](https://www.sibm.it/SITO%20MEDITS/Medits%20programme.htm)
- Diggle, P.J. & Ribeiro Jr, P.J. (2007). *Model-based Geostatistics*. Springer.
- Ribeiro Jr, P.J. & Diggle, P.J. (2001). geoR: A package for geostatistical analysis. *R-NEWS*, 1(2), 15-18.
