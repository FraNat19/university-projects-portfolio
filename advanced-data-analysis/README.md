# Advances in Data Analysis & Statistical Modelling

<p align="center">
  <img src="https://img.shields.io/badge/Tools-MATLAB%20%7C%20Python-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Course-ADASM-blue?style=for-the-badge" alt="Course">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Year-2024%2F2025-orange?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>Advanced Multivariate Analysis & Model-Based Composite Indicators on Real Data</b>
</p>

---

## Overview

This repository collects coursework and projects developed for **Advances in Data Analysis and Statistical Modelling** (Sapienza University of Rome, A.Y. 2024/2025).

Work includes:
- A **group project** on statistical model-based composite indicators applied to real socio-economic + crime-related data from U.S. communities
- An **individual mini-project** on clustering and dimensionality reduction on an OECD-like wellbeing dataset

---

## Course Objectives

The course focuses on building advanced skills to analyze complex, high-dimensional data and to support evidence-based decisions.

**Key learning goals include:**
- Reorganizing multidimensional data for statistical analysis
- Applying advanced multivariate methodologies with interpretable strategies
- Extracting relevant information from large datasets (including "big data")

---

## Projects & Methodologies

### 1) Model-Based Composite Indicator (Group Project)

A **composite indicator (CI)** is built by aggregating multiple manifest indicators (MIs) into a single index using an underlying statistical model of a multidimensional concept.

The project discusses:
- Reflective vs formative measurement perspectives
- The role of normalization/polarity
- Model selection strategies (confirmatory, exploratory, mixed)

#### Pipeline (High Level)

```
Correlation Matrix → Factor Extraction → Rotation → Loading Interpretation
         ↓
    EFA / DFA
         ↓
Specific Composite Indicators (SCIs)
         ↓
General Composite Indicator (GCI)
```

**Steps:**

1. **Correlation-based modelling** and goodness-of-fit reasoning for how indicators relate to latent constructs

2. **Exploratory Factor Analysis (EFA)**: correlation matrix → factor extraction → rotation → interpretation of loadings

3. **Disjoint Factor Analysis (DFA)**: constrains each variable to belong to exactly one factor (block structure), improving interpretability for building Specific Composite Indicators (SCIs)

4. **Aggregation** from SCIs to a General Composite Indicator (GCI) via a weighting scheme (statistical/equal/subjective)

#### Dataset (Real Data)

| Source | Description |
|--------|-------------|
| 1990 US Census | Socio-economic community conditions |
| 1990 LEMAS | Law enforcement survey data |
| 1995 FBI UCR | Crime statistics |

> **Source**: Communities and Crime dataset (UCI Machine Learning Repository)

#### Output Highlight

A **"Strong Social and Environmental Dynamics"** general index computed for communities, with example top-ranked cities (e.g., Newark, Camden, Miami) and bottom-ranked locations reported in the slides.

---

### 2) Clustering & Dimensionality Reduction (Personal Mini-Project)

Study on the **Better Life Index 2022** dataset (41 countries, 24 variables) using preprocessing + PCA + multiple clustering approaches.

#### Workflow

```
Raw Data → Missing Value Imputation (kNN) → Standardization → PCA → Clustering → Evaluation
```

**Preprocessing:**
- Missing-value imputation (kNN mean of nearest neighbors)
- Standardization

**PCA Selection:**
- Kaiser rule applied
- **7 PCs ≈ 81.53% variance** explained

**Clustering Evaluation:**
- Pseudo-F statistic for optimal K selection
- Best result: **K = 2**

#### Methods Implemented

| Method | Description |
|--------|-------------|
| **PCA** | Component selection and variance explanation |
| **k-means** | Standard clustering on component scores |
| **RKM** | Reduced K-Means |
| **FKM** | Factorial K-Means |
| **CDPCA** | Clustering and Disjoint PCA |
| **DKM** | Double K-Means |

**Comparison via:**
- Contingency tables
- Explained variance metrics

---

## 🛠️ Tools & Stack

| Tool | Purpose |
|------|---------|
| **MATLAB** | Factor analysis, clustering pipelines, reproducible scripts |
| **Python** | Data handling and experimentation |

---
## 👤 Author

**Francesco Natali** (Matricola 1945581)

---

## 📖 References

- Jolliffe, I.T. (2002). *Principal Component Analysis*. Springer.
- Hair, J.F. et al. (2019). *Multivariate Data Analysis*. Cengage.
- OECD (2008). *Handbook on Constructing Composite Indicators*.
- UCI Machine Learning Repository - Communities and Crime Dataset.

