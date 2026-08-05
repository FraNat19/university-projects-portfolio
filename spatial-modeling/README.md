<div align="center">

# 🦐 Spatial Prediction of Mediterranean Shrimp Biomass

**Bayesian and Maximum-Likelihood Kriging on MEDITS Trawl Survey Data**

<img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R">
<img src="https://img.shields.io/badge/geoR-kriging-1a7f37?style=for-the-badge" alt="geoR">
<img src="https://img.shields.io/badge/Sapienza-Spatial%20Statistics-8b1a1a?style=for-the-badge" alt="Sapienza">
<img src="https://img.shields.io/badge/A.Y.-2024%2F2025-555?style=for-the-badge" alt="Year">

<sub>MSc in Statistical Methods and Applications · Sapienza University of Rome</sub>

</div>

---

## 📋 Overview

Where do deep-water rose shrimp (*Parapenaeus longirostris*) actually live along the Italian Tyrrhenian coast, and how confident can we be about the places nobody sampled?

This project answers that with geostatistics. Starting from 120 trawl stations per survey year, it builds a spatial mixed-effect model of shrimp biomass, predicts it across the entire coastal grid from Liguria to Lazio, and quantifies the uncertainty of every prediction two different ways: **maximum-likelihood kriging** with confidence intervals, and **Bayesian kriging** with credible intervals.

The two approaches are then put head to head on RMSE and interval score across two survey years, 2002 and 2008.

| | |
|---|---|
| **Survey** | MEDITS (Mediterranean International Trawl Survey) |
| **Study area** | Tyrrhenian Sea, Liguria to Lazio |
| **Years** | 2002 and 2008 |
| **Stations** | 120 per year, 29 variables each |
| **Response** | Total biomass `tot`, kg/km² |
| **Methods** | IDW · variogram modelling · ML kriging · Bayesian kriging |
| **Deliverable** | [`Spatial_Presentation.pdf`](Spatial_Presentation.pdf) |

---

## 🗂️ Repository map

Eight weekly assignments. Four build the shrimp pipeline end to end; four apply the same toolkit to other datasets the course used to teach specific techniques.

| Script | Dataset | Technique |
|---|---|---|
| [`01_idw_interpolation.R`](01_idw_interpolation.R) | Wolfcamp aquifer | Inverse distance weighting, exponent *p* chosen by spatial cross-validation |
| [`02_exploratory_analysis_biomass.R`](02_exploratory_analysis_biomass.R) | 🦐 MEDITS shrimp | EDA, correlation structure, PCA for covariate selection |
| [`03_variogram_modelling.R`](03_variogram_modelling.R) | Wolfcamp aquifer | Trend removal, stationarity, empirical and fitted variograms |
| [`04_kriging_biomass_estimation.R`](04_kriging_biomass_estimation.R) | 🦐 MEDITS shrimp | ML kriging with covariates, prediction maps, 95% confidence bounds |
| [`05_bayesian_kriging.R`](05_bayesian_kriging.R) | 🦐 MEDITS shrimp | `krige.bayes`, credible intervals, comparison against script 04 |
| [`07_spatial_econometrics_sids.R`](07_spatial_econometrics_sids.R) | SIDS, North Carolina | Areal data: spatial weights, Moran's I, spatial regression |
| [`08_spatial_species_distribution.R`](08_spatial_species_distribution.R) | Mite data | Presence-absence modelling, ROC/AUC evaluation |
| [`09_extreme_value_theory_rainfall.R`](09_extreme_value_theory_rainfall.R) | Yearly rainfall | Block maxima, Generalized Extreme Value fitting |

**Data and outputs:** `shrimp2002.csv` / `shrimp2008.csv` (the 120 stations per year) · `ML.RData` / `Bayesian.rdata` (fitted kriging objects) · four `*_Kriging_*.pdf` prediction maps.

---

## 📊 The data, and what makes it hard

MEDITS runs standardised bottom trawls across the Mediterranean to monitor demersal and benthic fish, crustaceans and cephalopods. Each station carries 29 variables:

| Group | Variables |
|---|---|
| **Spatial** | Bathymetry (`bat`), distance from coast (`dist`), slope |
| **Environmental** | Temperature and salinity, quarterly minima and maxima |
| **Biological** | `Spawners` (mature adults), `Recruits` (juveniles), both kg/km² |

Two properties of this dataset shape every modelling decision that follows.

**Half the stations are empty.** Exactly 60 of 120 record zero biomass, in both years. On the original scale the distribution is heavily right-skewed with a spike at zero, which rules out a Gaussian likelihood and points toward a Tweedie. The route taken here is a log transform:

```r
# log(x+1) keeps the zeros at zero and leaves the non-zero values
# roughly symmetric, so the lognormal assumption becomes defensible
shrimp$log_tot <- log(shrimp$tot + 1)
```

**The obvious covariates are unusable.** `tot` is by construction `Spawners + Recruits`, so neither can enter the model without leaking the response. Predictors have to come from the environment instead.

<details>
<summary><b>Covariate selection via PCA</b> (click to expand)</summary>

<br>

The environmental variables are heavily collinear: the quarterly salinity measures correlate almost perfectly with each other, and so do the temperature ones. PCA on the standardised covariates puts roughly 68–70% of total variance in the first two components, and the highest-loading variables on those components become the model covariates.

The selection differs between years, which is itself a result worth noting:

| | 2002 | 2008 |
|---|---|---|
| **PC1** | `salinity.minq3` (0.2713) | `salinity.maxq3` (0.2565) |
| **PC2** | `dist` (0.2876), `bat` (0.2538), `temp.maxq3` (−0.4262) | `temp.maxq1` (0.3436), `slope` (0.2419), `bat` (−0.2288) |

Bathymetry loads strongly in both years. The temperature quarter that matters flips from Q3 to Q1.

</details>

A separate exploratory finding: 2002 shows medium-sized clusters of **recruits**, while 2008 shows large clusters of **spawners**. Overfishing pressure in 2002, predator pressure and life-cycle dynamics are all plausible explanations, none of which this dataset can separate. In both years biomass concentrates at 200–500 m depth.

---

## 🔬 Method

The spatial process is modelled as a mixed-effect model:

$$Z(s) = X(s)\beta + W(s) + \varepsilon(s)$$

| Term | Meaning |
|---|---|
| $X(s)\beta$ | Deterministic trend from the selected covariates |
| $W(s) \sim N_n(0, \sigma^2 H_{11}(\phi))$ | Spatially correlated random effect |
| $\varepsilon(s) \sim N_n(0, \tau^2)$ | Measurement error, the nugget |

### Maximum likelihood

Variogram parameters are estimated by likelihood, then prediction runs over the coastal grid:

```r
fit_exponential <- likfit(geodata, cov.model = "exponential",
                          ini.cov.pars = c(sigma2, phi), nugget = tau2)

kc <- krige.control(trend.d = "1st", trend.l = "1st", obj.model = fit_exponential)
krig  <- krige.conv(geodata, locations = grid, krige = kc)

lower <- krig$predict - 1.96 * sqrt(krig$krige.var)
upper <- krig$predict + 1.96 * sqrt(krig$krige.var)
```

### Bayesian

Non-informative priors throughout, with $\phi$ discretised over a range read off the empirical variogram:

```r
prior <- prior.control(
  beta.prior    = "normal",      # mean zero, large variance
  sigmasq.prior = "reciprocal",  # standard non-informative choice
  phi.prior     = "uniform",     # discretised: 25-50 (2002), 10-25 (2008)
  tausq.rel.prior = "uniform"
)

krig_bayes <- krige.bayes(geodata, locations = grid,
                          prior = prior, model = model.control(...))
```

The posteriors come out with clear peaks in both years, which is the check that the priors were wide enough to let the data speak.

### Validation

Both models are scored on 10 repeated train-test splits, using RMSE for point accuracy and the **interval score** for the quality of the uncertainty bands. The interval score rewards narrow intervals but penalises intervals that miss the true value:

```r
Int_Score <- function(Y, L, U, alpha = 0.05) {
  (U - L) + (2 / alpha) * ((L - Y) * (Y < L) + (Y - U) * (Y > U))
}
```

---

## 📈 Results

### Model comparison

Averages over 10 cross-validation repetitions. Lower is better on both metrics.

| Metric | ML kriging | Bayesian kriging |
|---|---|---|
| **RMSE 2002** | 2.2617 | **2.0759** |
| **RMSE 2008** | 2.3812 | **2.1481** |
| **Interval score 2002** | 9.731 | **8.384** |
| **Interval score 2008** | **9.708** | 9.950 |

The Bayesian model wins on point accuracy in both years and produces clearly better intervals in 2002. In 2008 the interval scores are effectively tied, with ML marginally ahead.

The honest reading: the Bayesian approach is the better of the two on this dataset, but the margin is small and it costs substantially more computation. What the priors bought here was a modest gain in accuracy, not a decisive one. On a dataset this size that trade is worth making; it would not obviously scale.

### Where the shrimp are

| Region | 2002 | 2008 |
|---|---|---|
| **Liguria** (north) | Low | Low |
| **Tuscany** (central) | Medium–high | Medium–high |
| **Lazio** (south) | High | Very high |

Both methods agree on the picture. Biomass concentrates along the Tuscany and Lazio coasts at 200–400 m depth, with hotspots around the Island of Elba and northern Lazio. Liguria stays poor throughout, which the bathymetry explains: the gradient there is steep, so the 200–400 m band that the species prefers is compressed into a narrow strip. Tuscany shows the most stable temperature and salinity profiles.

Between the two years the hotspots shift south-east and the high values spread over a wider area. The methods disagree most in 2008, where ML pushes the peak out toward the Lazio coast while the Bayesian fit keeps it concentrated centrally.

Prediction intervals are narrow in both years, and 2008 carries higher upper bounds, consistent with its heavier right tail.

---

## 🦐 About *Parapenaeus longirostris*

| Characteristic | Value |
|---|---|
| **Common name** | Deep-water rose shrimp |
| **Size** | 13–15 cm, females larger than males |
| **Depth** | 50–700 m, mainly 200–400 m |
| **Temperature** | 12–16 °C |
| **Salinity** | 35–39 PSU |
| **Substrate** | Soft, muddy or sandy-mud seabeds |
| **Why it matters** | Indicator species for climate change, habitat degradation and overfishing |

---

## ⚙️ Running the code

Requires R with `geoR`, `gstat`, `sp`, `sf`, `spdep`, `spatialreg`, `vegan`, `extRemes` and `ggplot2`.

> **Note on data.** Scripts 02, 04 and 05 load `shrimpsfull.RData` (the full multi-year survey panel) and `AllGrids.RData` (the prediction grids); script 09 loads `datiVE.RData`. These are course-provided datasets and are not redistributed here, so the shrimp scripts will not run end to end from a clone.
>
> What is included is enough to inspect the work: the two CSVs hold the 2002 and 2008 stations, and `ML.RData` and `Bayesian.rdata` hold the fitted kriging objects, so the maps can be reproduced without refitting. Scripts 01, 03, 07 and 08 run as they are, since their datasets ship with `geoR`, `spData` and `vegan`.

Scripts 04 and 05 are the slow ones. Bayesian kriging over the full grid takes a while.

---

## 📖 References

- **MEDITS Programme** — [sibm.it](https://www.sibm.it/SITO%20MEDITS/Medits%20programme.htm)
- Diggle, P.J. & Ribeiro Jr, P.J. (2007). *Model-based Geostatistics*. Springer.
- Ribeiro Jr, P.J. & Diggle, P.J. (2001). geoR: A package for geostatistical analysis. *R-NEWS*, 1(2), 15–18.
