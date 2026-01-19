# Multivariate Analysis with SAS (BSc Coursework)

<p align="center">
  <img src="https://img.shields.io/badge/Tools-SAS-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Methods-PCA%20%7C%20Clustering%20%7C%20Regression-blue?style=for-the-badge" alt="Methods">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Level-BSc-orange?style=for-the-badge" alt="Level">
</p>

<p align="center">
  <b>Hands-on multivariate workflows: exploration, dimensionality reduction, clustering, and modelling</b>
</p>

---

## Overview

This repository collects a set of **Bachelor-level assignments** developed in a Multivariate Data Analysis course using **SAS**, covering end-to-end workflows from descriptive statistics to multivariate modelling.

Each exercise:
- Starts from a real (or realistic) dataset
- Applies multivariate techniques
- Provides an interpretable statistical narrative

**Focus**: Relationships between variables and similarity patterns among statistical units.

---

## Case Studies exaples

### 1. European Diet Patterns (`dieta.pdf`)

**Dataset**: 16 European countries × 10 food consumption variables

| Variable | Unit | Description |
|----------|------|-------------|
| Cereali | kg | Cereals consumption |
| Riso | kg | Rice consumption |
| Patate | kg | Potatoes consumption |
| Zucchero | kg | Sugar consumption |
| Verdure | kg | Vegetables consumption |
| Vino | L | Wine consumption |
| Carne | kg | Meat consumption |
| Latte | L | Milk consumption |
| Burro | kg | Butter consumption |
| Uova | kg | Eggs consumption |

**Objectives:**
- Analyze and describe dietary habits across European countries
- Identify linear relationships between food consumption patterns
- Detect country similarities through dimensionality reduction

**Key Findings:**

| Food | Mean | Median | CV (%) | Interpretation |
|------|------|--------|--------|----------------|
| Milk | 128.4 L | 125.3 L | 35.8 | Most consumed, high variability |
| Vegetables | 95.6 kg | 82.5 kg | 57.5 | High heterogeneity (Greece outlier) |
| Wine | 25.9 L | 21.5 L | 75.9 | Highest variability across countries |
| Rice | 4.2 kg | 4.3 kg | 28.9 | Least consumed |

**PCA Results (4 components, 83.7% variance):**

| PC | Interpretation | Key Loadings |
|----|----------------|--------------|
| PC1 | Mediterranean diet | +Rice, +Vegetables, +Wine, −Sugar, −Milk |
| PC2 | Protein-rich diet | +Meat, +Eggs, +Butter, −Milk, −Rice |
| PC3 | Grain-based diet | +Cereals, +Vegetables, −Butter |
| PC4 | Starch indicator | +Potatoes |

**Country Profiles:**
- 🇮🇹🇪🇸🇵🇹 **Mediterranean cluster**: High PC1 scores (rice, vegetables, wine)
- 🇫🇷🇩🇪🇩🇰 **Central European cluster**: High PC2 scores (meat, eggs, butter)
- 🇮🇪 **Ireland**: Outlier on PC4 (highest potato consumption)

---

### 2. Life Expectancy Analysis (`countries_lifeexp.pdf`)

**Dataset**: 38 countries × 12 socio-economic indicators

| Variable | Description |
|----------|-------------|
| Region | Geographic region |
| Area | Surface area |
| Irrigated | % irrigated land |
| Pop | Population density |
| Under14 | % population under 14 |
| LifeExp | Life expectancy at birth |
| Literacy | Literacy rate |
| Unemployment | Unemployment rate |
| ISP | Social institutions per million |
| TV | TVs per person |
| Railways | Railway km / area |
| Airports | Airports / area |

**Objectives:**
1. Analyze life expectancy determinants via PCA and clustering
2. Build predictive regression model for life expectancy
3. Evaluate multicollinearity and influential observations

**Methods Applied:**
- Principal Component Analysis with Varimax rotation
- Hierarchical clustering (Ward's method)
- Multiple linear regression with diagnostics

---

### 3. Labor Force Analysis (`forza_lavoro.pdf`)

**Dataset**: Regional labor force statistics

**Objectives:**
- Study employment patterns across regions
- Identify structural differences in labor markets
- Cluster regions by workforce characteristics

---

## Techniques Applied

### Analysis Pipeline

```
┌─────────────────┐
│  Raw Data       │
└────────┬────────┘
         ↓
┌─────────────────┐
│ Descriptive     │  PROC MEANS, PROC UNIVARIATE
│ Statistics      │  PROC CORR
└────────┬────────┘
         ↓
┌─────────────────┐
│ Standardization │  PROC STANDARD
│ (if needed)     │
└────────┬────────┘
         ↓
┌─────────────────┐
│ Dimensionality  │  PROC PRINCOMP
│ Reduction (PCA) │  PROC FACTOR (+ Varimax)
└────────┬────────┘
         ↓
┌─────────────────┐
│ Clustering      │  PROC CLUSTER (Ward)
│                 │  PROC FASTCLUS (k-means)
└────────┬────────┘
         ↓
┌─────────────────┐
│ Regression      │  PROC REG
│ Modelling       │  (diagnostics, VIF, influence)
└─────────────────┘
```

### SAS Procedures Used

| Procedure | Purpose |
|-----------|---------|
| `PROC MEANS` | Descriptive statistics |
| `PROC UNIVARIATE` | Distribution analysis |
| `PROC CORR` | Correlation matrix |
| `PROC STANDARD` | Data standardization |
| `PROC PRINCOMP` | Principal Component Analysis |
| `PROC FACTOR` | Factor Analysis with rotation |
| `PROC CLUSTER` | Hierarchical clustering |
| `PROC FASTCLUS` | K-means clustering |
| `PROC REG` | Linear regression |
| `PROC SGPLOT` | Visualization |


---

## Key Methodological Concepts

### When to Standardize?

| Situation | Standardize? |
|-----------|--------------|
| Variables in different units (kg vs L vs %) | ✅ Yes |
| Variables with very different variances | ✅ Yes |
| All variables in same scale, similar variance | ❌ Optional |
| Clustering on mixed-scale data | ✅ Yes |

### PCA: How Many Components?

| Criterion | Rule |
|-----------|------|
| Kaiser | Eigenvalue > 1 |
| Scree plot | Elbow in variance curve |
| Cumulative variance | ≥ 70-80% explained |
| Interpretability | Components must make sense |

### Clustering: How Many Groups?

| Method | Criterion |
|--------|-----------|
| Dendrogram | Visual cut at large distance jump |
| CCC (Cubic Clustering Criterion) | Local peaks |
| Pseudo-F | Higher is better |
| Silhouette | Maximize average silhouette |

---

## What These Examples Demonstrate

| Skill | Application |
|-------|-------------|
| **Dimensionality Reduction** | Reduce 10 food variables to 4 interpretable components |
| **Pattern Discovery** | Identify Mediterranean vs Central European diets |
| **Clustering** | Group countries by consumption similarity |
| **Predictive Modelling** | Explain life expectancy from socio-economic indicators |
| **Diagnostics** | Detect outliers, leverage points, multicollinearity |

---

## 👤 Author

**Francesco Natali** (1945581)

---
  <i>Sapienza Università di Roma • Multivariate Data Analysis • BSc Program</i>
</p>
