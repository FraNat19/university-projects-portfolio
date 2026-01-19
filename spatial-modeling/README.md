# 🦐 Spatial Statistics Analysis of Mediterranean Shrimp Biomass

<p align="center">
  <img src="https://img.shields.io/badge/Language-R-blue?style=for-the-badge&logo=r" alt="R">
  <img src="https://img.shields.io/badge/Course-Spatial%20Statistics-red?style=for-the-badge" alt="Course">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Year-2024%2F2025-green?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>A Statistical Framework for Understanding Parapenaeus longirostris:<br>Bayesian and Maximum Likelihood Applications</b>
</p>

---

## 📋 Overview

This repository contains the complete work developed during the **Spatial Statistics and Statistical Tools for Environmental Data** course at Sapienza University of Rome (A.Y. 2024/2025). The project spans an entire semester of weekly assignments, focusing on the spatial analysis of **deep-water rose shrimp (*Parapenaeus longirostris*)** biomass along the Italian Tyrrhenian coast.

---

## Project Objectives

The main goals of this project were to:

1. **Explore spatial patterns** in Mediterranean shrimp biomass distribution
2. **Apply geostatistical methods** for spatial prediction and interpolation
3. **Compare frequentist and Bayesian approaches** to kriging
4. **Quantify uncertainty** in spatial predictions through confidence/credible intervals
5. **Investigate environmental drivers** affecting shrimp distribution (temperature, salinity, bathymetry)

---

## 📊 Dataset

### MEDITS Survey
The data comes from the **MEDITS (Mediterranean International Trawl Survey)** program, which conducts standardized trawl surveys across the Mediterranean Sea to monitor demersal species.

| Feature | Description |
|---------|-------------|
| **Study Area** | Tyrrhenian Sea (Liguria to Lazio coast) |
| **Years Analyzed** | 2002 and 2008 |
| **Target Species** | *Parapenaeus longirostris* (deep-water rose shrimp) |
| **Variables** | 29 including spatial, environmental, and biological data |
| **Response Variable** | Total biomass (kg/km²) |

### Key Covariates
- **Spatial**: Bathymetry, distance from coast, slope
- **Environmental**: Temperature (seasonal quartiles), Salinity (seasonal quartiles)
- **Biological**: Spawners (adults), Recruits (juveniles)

---

## 🔬 Methodology

### 1. Exploratory Data Analysis
- Principal Component Analysis (PCA) for covariate selection
- Correlation analysis and multivariate exploration
- Log-transformation of biomass data to handle zero-inflation and skewness

### 2. Variogram Modeling
- Empirical variogram estimation
- Model fitting: **Matérn**, **Exponential**, **Spherical**
- Model selection via cross-validation (RMSE)

### 3. Spatial Interpolation

#### Maximum Likelihood Kriging
```r
# Variogram model fitting
fit_exponential <- likfit(geodata, cov.model = "exponential", 
                          ini.cov.pars = c(sigma2, phi), nugget = tau2)

# Kriging interpolation
krig_result <- krige.conv(geodata, locations = grid, krige = krige.control(...))
```

#### Bayesian Kriging
```r
# Prior specification
prior <- list(
  beta.prior = "normal",
  sigmasq.prior = "reciprocal",
  phi.prior = "uniform",
  tausq.rel.prior = "uniform"
)

# Bayesian kriging
krige_bayes <- krige.bayes(geodata, locations = grid, prior = prior, model = model)
```

### 4. Uncertainty Quantification
- **Frequentist**: 95% Confidence Intervals
- **Bayesian**: 95% Credible Intervals from posterior distributions
- Interval Score comparison for model evaluation

---

## Key Results

### Biomass Distribution Patterns

| Region | 2002 | 2008 |
|--------|------|------|
| **Liguria (North)** | Low biomass | Low biomass |
| **Tuscany (Central)** | Medium-High | Medium-High |
| **Lazio (South)** | High | Very High |

### Model Performance Comparison

| Metric | MLE Kriging | Bayesian Kriging |
|--------|-------------|------------------|
| **RMSE 2002** | 2.26 | **2.08** ✓ |
| **RMSE 2008** | 2.38 | **2.15** ✓ |
| **Interval Score 2002** | 9.73 | **8.38** ✓ |
| **Interval Score 2008** | 9.71 | 9.95 |

> **Conclusion**: The Bayesian approach consistently outperforms MLE in terms of RMSE, demonstrating its ability to incorporate prior information and handle uncertainty more effectively.

### Spatial Predictions

<p align="center">
  <i>Interpolated shrimp biomass shows highest concentrations in the central-southern Tyrrhenian (Tuscany-Lazio coast), with optimal conditions at 200-400m depth.</i>
</p>


---

## Key Concepts

### Spatial Mixed Effect Model
The spatial process is modeled as:

$$Z(s) = X(s)\beta + W(s) + \varepsilon(s)$$

Where:
- $X(s)\beta$ : Deterministic trend (fixed effects)
- $W(s) \sim N_n(0, \sigma^2 H_{11}(\phi))$ : Spatial random effect
- $\varepsilon(s) \sim N_n(0, \tau^2)$ : Measurement error (nugget)

### Bayesian Framework
Prior distribution:
$$\pi(\beta, \omega) = \pi(\beta)\pi(\tau^2)\pi(\sigma^2, \phi)$$

Posterior distribution obtained via MCMC sampling for prediction at unobserved locations.

---

## About *Parapenaeus longirostris*

| Characteristic | Description |
|----------------|-------------|
| **Common Name** | Deep-water rose shrimp |
| **Habitat Depth** | 50-700m (optimal: 200-400m) |
| **Temperature Range** | 12-16°C |
| **Salinity Preference** | 35-39 PSU |
| **Substrate** | Soft, muddy seabeds |
| **Ecological Role** | Indicator species for climate change and overfishing |

---

## 📖 References

- **MEDITS Programme**: [www.sibm.it/MEDITS](https://www.sibm.it/SITO%20MEDITS/Medits%20programme.htm)
- Diggle, P.J. & Ribeiro Jr, P.J. (2007). *Model-based Geostatistics*. Springer.
- Ribeiro Jr, P.J. & Diggle, P.J. (2001). geoR: A package for geostatistical analysis. *R-NEWS*, 1(2), 15-18.

