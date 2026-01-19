# Multiple Time Series Modelling (MTSM)

<p align="center">
  <img src="https://img.shields.io/badge/Tools-R%20%7C%20MATLAB-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Topics-VAR%20%7C%20IRF%20%7C%20Unit%20Roots%20%7C%20Cointegration-blue?style=for-the-badge" alt="Topics">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Year-2024%2F2025-orange?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>Hands-on assignments: stability, identification, spurious regressions, and long-run relationships</b>
</p>

---

## Overview

This repository contains coursework for **Multiple Time Series Modelling (MTSM)**, covering both simulation-based exercises on VAR systems and applied replication studies in econometrics/finance.

**Focus areas:**
- Understanding multivariate dynamics
- Diagnosing non-stationarity pitfalls
- Implementing standard tools for inference and interpretation in multi-series settings

---

## Course Objectives

| Objective | Description |
|-----------|-------------|
| **VAR Modelling** | Estimate and interpret Vector Autoregressive models |
| **Impulse Responses** | Analyze dynamic effects of shocks through IRFs |
| **Spurious Regression** | Understand and avoid pitfalls with non-stationary data |
| **Cointegration** | Test and model long-run equilibrium relationships |

---

## Topics & Methods

### 1. VAR Models & Identification

#### VAR(1) Representation

$$Y_t = c + A Y_{t-1} + u_t, \quad u_t \sim N(0, \Omega)$$

Where:
- $Y_t$ = $k \times 1$ vector of endogenous variables
- $A$ = $k \times k$ coefficient matrix
- $\Omega$ = $k \times k$ covariance matrix of reduced-form shocks

#### Stability Conditions

```
Stability Analysis via Eigenvalues of A
────────────────────────────────────────────────────

Let λ₁, λ₂, ..., λₖ be eigenvalues of A

┌─────────────────────────────────────────────────┐
│ |λᵢ| < 1  ∀i  →  STATIONARY (stable)           │
│ |λᵢ| = 1  ∃i  →  UNIT ROOT (non-stationary)    │
│ |λᵢ| > 1  ∃i  →  EXPLOSIVE (unstable)          │
└─────────────────────────────────────────────────┘
```

**Simulation scenarios tested:**

| Scenario | Eigenvalues | Behavior |
|----------|-------------|----------|
| Stationary | $|\lambda| < 1$ | Mean-reverting, stable variance |
| Unit root | $|\lambda| = 1$ | Random walk, growing variance |
| Explosive | $|\lambda| > 1$ | Divergent paths |

#### Impulse Response Functions (IRFs)

IRFs trace the effect of a one-unit shock to variable $j$ on variable $i$ over time:

$$IRF_{i,j}(h) = \frac{\partial Y_{i,t+h}}{\partial u_{j,t}}$$

**Cholesky Identification:**

$$\Omega = PP', \quad P = \text{lower triangular (Cholesky factor)}$$

```
Structural shocks: ε_t = P⁻¹ u_t

IRF computation:
─────────────────
h = 0:  Ψ₀ = P
h = 1:  Ψ₁ = A·P
h = 2:  Ψ₂ = A²·P
  ⋮
h = H:  Ψₕ = Aʰ·P
```

**Important caveat**: Cholesky decomposition imposes a **recursive ordering**—the first variable responds contemporaneously only to its own shock. Different orderings yield different IRFs!

#### Rotation Invariance Problem

Multiple decompositions satisfy $\Omega = PP'$:

$$\Omega = PP' = (PQ)(PQ)' \quad \text{for any orthogonal } Q$$

**Implication**: Without additional economic restrictions, structural identification is not unique.

---

### 2. Spurious Regressions & Non-Stationarity

#### The Granger-Newbold (1974) Problem

When regressing **unrelated non-stationary** series:

```
Monte Carlo Setup (Granger-Newbold replication)
───────────────────────────────────────────────
Generate independent random walks:
  x_t = x_{t-1} + ε_t,  ε_t ~ N(0,1)
  y_t = y_{t-1} + η_t,  η_t ~ N(0,1)

Regress: y_t = α + β·x_t + u_t

SPURIOUS RESULTS:
  • High R² (often > 0.5)
  • Low Durbin-Watson (near 0)
  • Highly significant t-statistics
  • BUT β should be 0 (no true relationship!)
```

**Diagnostic symptoms of spurious regression:**

| Symptom | Typical Value | Indicates |
|---------|---------------|-----------|
| $R^2$ | High (> 0.5) | Apparent explanatory power |
| DW | Very low (< 1) | Severe autocorrelation in residuals |
| t-stat | Very high | Inflated significance |
| Residual ACF | Persistent | Non-stationary residuals |

#### Solution: First Differencing

$$\Delta y_t = y_t - y_{t-1}$$

Differencing converts I(1) series to I(0):
- Removes spurious correlation
- Restores valid inference
- But loses long-run information!

---

### 3. Cointegration & Long-Run Relationships

#### Concept

Two I(1) series are **cointegrated** if a linear combination is I(0):

$$y_t \sim I(1), \quad x_t \sim I(1)$$
$$y_t - \beta x_t = u_t \sim I(0)$$

**Interpretation**: Short-run deviations from equilibrium, but long-run relationship exists.

#### Testing Procedures

**a) Engle-Granger (1987) Two-Step:**

```
Step 1: Estimate cointegrating regression
        y_t = α + β·x_t + u_t

Step 2: Test residuals for stationarity
        ADF test on û_t
        
H₀: No cointegration (residuals are I(1))
H₁: Cointegration (residuals are I(0))
```

**b) Phillips-Ouliaris:**

Similar to EG but uses Phillips-Perron test statistics (robust to heteroskedasticity/autocorrelation).

**c) Johansen (1988) Trace Test:**

System-based approach for multivariate cointegration:

$$\Delta Y_t = \Pi Y_{t-1} + \sum_{i=1}^{p-1} \Gamma_i \Delta Y_{t-i} + u_t$$

Where $\Pi = \alpha \beta'$ and $rank(\Pi) = r$ = number of cointegrating vectors.

```
Johansen Trace Test
───────────────────
H₀: r = 0  vs  H₁: r > 0   (test for at least 1 CI vector)
H₀: r ≤ 1  vs  H₁: r > 1   (test for at least 2 CI vectors)
  ⋮
```

---

### 4. Applied Analysis: MSCI World Indices

**Dataset**: MSCI World and sector indices (daily/monthly)

**Analysis pipeline:**

```
┌──────────────────┐
│ Raw Price Data   │
└────────┬─────────┘
         ↓
┌──────────────────┐
│ Unit Root Tests  │  ADF, PP, KPSS
│ (Levels)         │  → Confirm I(1)
└────────┬─────────┘
         ↓
┌──────────────────┐
│ First Differences│  Compute returns
│ (Returns)        │  → Confirm I(0)
└────────┬─────────┘
         ↓
┌──────────────────┐
│ Cointegration    │  Pairwise tests
│ Tests            │  Johansen system test
└────────┬─────────┘
         ↓
┌──────────────────┐
│ VECM if          │  Error correction model
│ cointegrated     │  if CI found
└──────────────────┘
```

---

## Key Results Summary

### VAR(1) Simulation

| Setting | Eigenvalues | Stability | IRF Behavior |
|---------|-------------|-----------|--------------|
| Baseline | (0.5, 0.3) | Stationary | Decaying to zero |
| Near unit root | (0.95, 0.8) | Stationary | Slow decay |
| Unit root | (1.0, 0.5) | Non-stationary | Persistent |

### Spurious Regression (Monte Carlo)

| Metric | True Value | Spurious Result |
|--------|------------|-----------------|
| $\beta$ | 0 | Often significant |
| $R^2$ | ~0 | 0.3 - 0.8 |
| DW | ~2 | < 0.5 |
| Rejection rate (5%) | 5% | > 70% |

### Cointegration Tests (Shiller Data)

| Test | Null Hypothesis | Result |
|------|-----------------|--------|
| Engle-Granger | No cointegration | Reject at 5% |
| Phillips-Ouliaris | No cointegration | Reject at 5% |
| Johansen (r=0) | No CI vectors | Reject |
| Johansen (r≤1) | At most 1 CI vector | Fail to reject |

---

## Tools & Stack

| Tool | Purpose |
|------|---------|
| **R** | Econometric analysis, testing, visualization |
| `vars` | VAR estimation and IRFs |
| `urca` | Unit root and cointegration tests |
| `tseries` | Time series utilities |
| `MASS` | Multivariate simulation |
| **MATLAB** | Alternative for simulations |
---

## Key concepts

| Concept | Takeaway |
|---------|----------|
| **VAR Stability** | Eigenvalues inside unit circle ensure stationarity |
| **IRF Identification** | Cholesky ordering matters—economic theory should guide |
| **Spurious Regression** | Never regress I(1) series in levels without checking |
| **Cointegration** | Allows modeling long-run relationships between non-stationary series |
| **Testing Strategy** | Always test for unit roots before regression analysis |

---

## 📚 References

- Hamilton, J.D. (1994). *Time Series Analysis*. Princeton University Press.
- Lütkepohl, H. (2005). *New Introduction to Multiple Time Series Analysis*. Springer.
- Granger, C.W.J. & Newbold, P. (1974). Spurious Regressions in Econometrics. *Journal of Econometrics*.
- Engle, R.F. & Granger, C.W.J. (1987). Co-Integration and Error Correction. *Econometrica*.
- Johansen, S. (1988). Statistical Analysis of Cointegration Vectors. *Journal of Economic Dynamics and Control*.

  <br>
  <i>Sapienza Università di Roma • Multiple Time Series Modelling • A.Y. 2024/2025</i>
</p>
