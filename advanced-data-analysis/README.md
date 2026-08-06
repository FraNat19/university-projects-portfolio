<div align="center">

# 🏙️ Advances in Data Analysis & Statistical Modelling

**Building a crime index for 1,994 US communities, and clustering 41 countries by quality of life**

<img src="https://img.shields.io/badge/MATLAB-0076A8?style=for-the-badge&logo=mathworks&logoColor=white" alt="MATLAB">
<img src="https://img.shields.io/badge/Factor%20Analysis-EFA%20%7C%20DFA-1a7f37?style=for-the-badge" alt="Methods">
<img src="https://img.shields.io/badge/Sapienza-8b1a1a?style=for-the-badge" alt="Sapienza">

<sub>MSc in Statistical Methods and Applications · Sapienza University of Rome · A.Y. 2024/2025</sub>

</div>

---

## 📋 Overview

Two projects on the same underlying question: when you have dozens of correlated measurements of something you cannot observe directly — how safe a neighbourhood is, how good a country is to live in — how do you compress them into a number without throwing away what matters?

| | Project | Data | Code |
|---|---|---|---|
| 1 | Model-based composite indicator for community crime | 1,994 US communities, 28 variables | [`CRIME.m`](CRIME.m), [`mapcrime.m`](mapcrime.m) |
| 2 | Clustering and dimensionality reduction | Better Life Index 2022, 41 countries, 24 variables | in the slides |

---

## 🔬 Project 1 — A composite indicator for community crime

### The data

The **Communities and Crime** dataset from the UCI repository, joining three sources:

| Source | Contributes |
|---|---|
| 1990 US Census | Socio-economic conditions: income, poverty, education, unemployment, urbanisation, demographics |
| 1990 LEMAS survey | Law enforcement resources and staffing |
| 1995 FBI UCR | Crime statistics |

[`NEWcrimedata (1).xlsx`](NEWcrimedata%20\(1\).xlsx) holds **1,994 communities across 28 variables**, including `racepctblack`, `racePctWhite`, `pctUrban`, `medIncome`, `pctWPubAsst`, `PctPopUnderPov`, `PctNotHSGrad` and `PctUnemployed`.

### The pipeline

```
standardise → drop collinear (|r| > 0.99) → PCA scree → Kaiser rule
      → EFA → Disjoint Factor Analysis → Cronbach's alpha
      → factor scores → weighted composite → normalise to [0,1] → map
```

**Collinearity is removed before anything else.** Any variable correlating above 0.99 with another is dropped, because factor analysis on a near-singular correlation matrix produces loadings that are numerically meaningless.

```matlab
threshold = 0.99;
redundant_vars = any(abs(correlation_matrix - eye(size(correlation_matrix))) > threshold, 1);
cleaned_data = standardized_data(:, ~redundant_vars);
```

**The number of factors is chosen, not assumed.** A scree plot is produced, then Kaiser's criterion counts eigenvalues above 1:

```matlab
[coeff, ~, latent] = pca(standardized_data);
eigenvalues = sort(eig(correlation_matrix), 'descend');
num_factors = sum(eigenvalues > 1);
[loadings, psi, stats] = factoran(cleaned_data, num_factors);
```

**Disjoint Factor Analysis is the methodological point.** Standard EFA lets every variable load on every factor, which makes interpretation a judgement call about where to draw the line. DFA constrains each variable to belong to **exactly one** factor, producing a block structure. The cost is fit; the gain is that each resulting Specific Composite Indicator has an unambiguous meaning.

**Reliability is checked before the index is trusted.** Cronbach's alpha is computed from the item variances against the variance of the summed score, so internal consistency is measured rather than assumed.

**The composite is a weighted sum of factor scores**, normalised to [0, 1] and written back alongside the community names.

### The map

[`mapcrime.m`](mapcrime.m) plots the resulting index on a US map, colour-coding each community by intensity with the `hot` colormap over state boundaries:

```matlab
geoshow(ax, states, 'FaceColor', [0.8 0.8 0.8]);
scatterm(latitude(i), longitude(i), 50, colors(colorIndices(i), :), 'filled', ...
         'DisplayName', communityname{i});
```

The general index is interpreted as **"Strong Social and Environmental Dynamics"**, with Newark, Camden and Miami among the top-ranked communities.

### Reproducibility notes

The two scripts do not run as they stand, for three reasons worth stating plainly:

| Missing | Used by | What it is |
|---|---|---|
| `CRIMEd.xlsx` | `CRIME.m` | The input table the script reads; the repository ships `NEWcrimedata (1).xlsx` instead |
| `DFA.m` | `CRIME.m` | The Disjoint Factor Analysis routine, provided by the course |
| `MapCRIME.xlsx` | `mapcrime.m` | The joined table of index, latitude and longitude |

There is also duplication to clean up: both scripts carry the same ~47 lines of hardcoded community names and coordinates, roughly 200 lines in total across the two files. That block belongs in the data file, not in the source.

---

## 🌍 Project 2 — Clustering countries on the Better Life Index

The OECD **Better Life Index 2022**: 41 countries scored on 24 variables covering housing, income, jobs, community, education, environment, civic engagement, health, life satisfaction, safety and work-life balance.

```
raw data → kNN imputation → standardise → PCA → clustering → validation
```

**Dimensionality.** Kaiser's rule retains **7 principal components, capturing 81.53% of total variance.**

**Choosing k.** The pseudo-F statistic across candidate partitions selects **k = 2**. Forty-one countries split into two groups rather than the four or five a regional intuition would suggest.

**Four clustering approaches compared**, which is the interesting part, because they disagree by construction:

| Method | What it optimises |
|---|---|
| **k-means** on component scores | Clusters, taking the PCA as given |
| **Reduced K-Means (RKM)** | Subspace and clusters jointly, favouring between-cluster separation |
| **Factorial K-Means (FKM)** | Subspace and clusters jointly, favouring within-cluster compactness |
| **Clustering and Disjoint PCA (CDPCA)** | Clusters plus a disjoint component structure, so components stay interpretable |

Running PCA first and clustering second assumes the directions of greatest variance are also the directions that separate groups. RKM, FKM and CDPCA drop that assumption and solve both problems at once, which is exactly the same idea as DFA in Project 1: constrain the model so the output means something.

---

## 📂 Files

| File | Contents |
|---|---|
| [`CRIME.m`](CRIME.m) | Full pipeline: standardisation, collinearity filter, PCA, EFA, DFA, Cronbach's alpha, composite index |
| [`mapcrime.m`](mapcrime.m) | Geographic visualisation of the index across US communities |
| [`NEWcrimedata (1).xlsx`](NEWcrimedata%20\(1\).xlsx) | 1,994 communities, 28 variables |
| `ADASM 24-25.pdf`, `ADASM_1945581 (1).pdf` | Course slides and project report |

Requires MATLAB with the Statistics and Machine Learning Toolbox (`pca`, `factoran`, `zscore`) and the Mapping Toolbox (`geoshow`, `scatterm`).

---

## 📖 References

- Vichi, M. & Saporta, G. (2009). Clustering and disjoint principal component analysis. *Computational Statistics & Data Analysis*, 53(8), 3194–3208.
- Vichi, M. & Kiers, H.A.L. (2001). Factorial k-means analysis for two-way data. *Computational Statistics & Data Analysis*, 37(1), 49–64.
- Jolliffe, I.T. (2002). *Principal Component Analysis*. Springer.
- OECD (2008). *Handbook on Constructing Composite Indicators*.
- Redmond, M. & Baveja, A. (2002). A data-driven software tool for enabling cooperative information sharing among police departments. *European Journal of Operational Research*, 141(3), 660–678.
