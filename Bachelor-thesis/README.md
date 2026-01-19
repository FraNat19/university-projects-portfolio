# 🎓 Bachelor's Thesis — Advanced Statistical Inference

<p align="center">
  <img src="https://img.shields.io/badge/Tools-R%20%7C%20LaTeX-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Topic-Binomial%20Inference%20%7C%20Confidence%20Intervals-blue?style=for-the-badge" alt="Topic">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/A.Y.-2022%2F2023-orange?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>Approximate confidence intervals for the binomial model: new proposals and empirical comparisons</b>
</p>

---

## Overview

This repository contains the work from my **Bachelor's thesis**, focused on statistical inference for the binomial model through the construction and evaluation of **approximate confidence intervals** for the success probability $p$.

The thesis:
- Revisits the classical **Wald interval** and investigates why it can behave poorly in finite samples
- Studies **alternative intervals** (Wilson, Agresti-Coull, Andersson-Nerman)
- Explores improved approaches motivated by refined approximate pivotal quantities

---

## 👤 Author

**Francesco Natali** (1945581)

---

## 🏛️ Academic Context

| | |
|---|---|
| **Degree** | B.Sc. in Statistics, Economics, Finance and Insurance |
| **University** | Sapienza University of Rome |

---

## 🎯Research Goal

The main goal is to **study and improve confidence interval estimation** for a binomial parameter $p$—a classic but still subtle problem due to:

- **Discreteness** of the binomial distribution
- **Finite-sample effects** that asymptotic theory ignores

**Key Focus**: Understanding how "easy" asymptotic constructions (especially Wald-type intervals) may lead to **unstable coverage**, and how alternative approximations can yield more reliable inference.

```
Research Question:
─────────────────────────────────────────────────────────────
Can we construct confidence intervals for binomial p that:
  ✓ Achieve nominal coverage more reliably?
  ✓ Avoid the "lucky/unlucky n" phenomenon?
  ✓ Maintain reasonable expected length?
─────────────────────────────────────────────────────────────
```

---

## 🔬 Methods & Content

### The Problem: Binomial Confidence Intervals

Given $X \sim \text{Binomial}(n, p)$, construct an interval $[L(X), U(X)]$ such that:

$$P_p(L(X) \leq p \leq U(X)) \geq 1 - \alpha \quad \forall p \in (0,1)$$

### 1. Wald Interval

The most common "textbook" interval:

$$\hat{p} \pm z_{\alpha/2} \sqrt{\frac{\hat{p}(1-\hat{p})}{n}}$$

Where $\hat{p} = X/n$ is the MLE.

**Problems identified:**

| Issue | Description |
|-------|-------------|
| **Coverage oscillations** | Actual coverage varies wildly with $p$ |
| **Lucky/unlucky $n$** | Some sample sizes perform much better than others |
| **Boundary issues** | Degenerates when $\hat{p} = 0$ or $\hat{p} = 1$ |
| **Undercoverage** | Often falls below nominal level |



### 2. Wilson (Score) Interval

Derived by inverting the score test:

$$\frac{\hat{p} + \frac{z^2}{2n} \pm z\sqrt{\frac{\hat{p}(1-\hat{p})}{n} + \frac{z^2}{4n^2}}}{1 + \frac{z^2}{n}}$$

**Advantages:**
- Never degenerates at boundaries
- More stable coverage
- Respects $(0,1)$ parameter space

### 3. Agresti-Coull Interval

Simple adjustment using "pseudo-counts":

$$\tilde{p} = \frac{X + 2}{n + 4}$$

Then apply Wald formula with $\tilde{p}$ and $\tilde{n} = n + 4$.

**Interpretation**: Add 2 successes and 2 failures before computing the interval.

**Advantages:**
- Easy to compute and explain
- Excellent coverage properties
- "Add 2 successes and 2 failures" rule

### 4. Andersson-Nerman Approach

Motivated by reducing dependence between $\hat{p}$ and its estimated variance through a **correlation-adjustment** idea.

**Key insight**: The correlation between $\hat{p}$ and $\widehat{Var}(\hat{p})$ causes distortions in the pivotal quantity.

$$\rho(\hat{p}, \hat{\sigma}^2) \neq 0 \Rightarrow \text{Pivot is non-standard}$$

**Correction**: Adjust the pivot to account for this correlation.

---

## Evaluation Criteria

### 1. Coverage Probability

$$CP(p) = P_p(p \in CI) = \sum_{x=0}^{n} \mathbf{1}[p \in CI(x)] \binom{n}{x} p^x (1-p)^{n-x}$$

**Desideratum**: $CP(p) \approx 1 - \alpha$ for all $p$.

### 2. Non-Coverage Components

| Component | Definition |
|-----------|------------|
| **Mesial** | $P(p < L)$ — interval too far right |
| **Distal** | $P(p > U)$ — interval too far left |

Balanced non-coverage: $P(p < L) \approx P(p > U) \approx \alpha/2$

### 3. Expected Length

$$EL(p) = E_p[U(X) - L(X)]$$

**Trade-off**: Better coverage often requires wider intervals.

### 4. Pivotal Diagnostics

For pivot $Z = \frac{\hat{p} - p}{\hat{\sigma}}$, compare:

| Moment | Target (Normal) | Actual |
|--------|-----------------|--------|
| Mean | 0 | $E[Z]$ |
| Variance | 1 | $Var(Z)$ |
| Skewness | 0 | $\gamma_1(Z)$ |

---

## 🧪 Computational Experiments (R)

### Simulation Design

```r
# Parameters
n_values <- c(10, 20, 30, 50, 100, 200)
p_grid <- seq(0.001, 0.999, by = 0.001)
alpha <- 0.05
n_sim <- 10000  # For Monte Carlo coverage
```

### Coverage Computation

```r
# Exact coverage probability (no simulation needed for binomial)
coverage_exact <- function(n, p, ci_function, alpha = 0.05) {
  x_values <- 0:n
  probs <- dbinom(x_values, n, p)
  
  covered <- sapply(x_values, function(x) {
    ci <- ci_function(x, n, alpha)
    p >= ci[1] & p <= ci[2]
  })
  
  sum(probs * covered)
}
```

### Interval Functions

```r
# Wald interval
ci_wald <- function(x, n, alpha = 0.05) {
  p_hat <- x / n
  z <- qnorm(1 - alpha/2)
  se <- sqrt(p_hat * (1 - p_hat) / n)
  c(max(0, p_hat - z * se), min(1, p_hat + z * se))
}

# Wilson interval
ci_wilson <- function(x, n, alpha = 0.05) {
  z <- qnorm(1 - alpha/2)
  p_hat <- x / n
  denom <- 1 + z^2 / n
  center <- (p_hat + z^2 / (2*n)) / denom
  margin <- z * sqrt(p_hat*(1-p_hat)/n + z^2/(4*n^2)) / denom
  c(center - margin, center + margin)
}

# Agresti-Coull interval
ci_agresti_coull <- function(x, n, alpha = 0.05) {
  z <- qnorm(1 - alpha/2)
  n_tilde <- n + z^2
  p_tilde <- (x + z^2/2) / n_tilde
  se <- sqrt(p_tilde * (1 - p_tilde) / n_tilde)
  c(max(0, p_tilde - z * se), min(1, p_tilde + z * se))
}
```

### Example Results

**Coverage curves (n = 30):**

| Method | Min Coverage | Mean Coverage | Max Coverage |
|--------|--------------|---------------|--------------|
| Wald | 0.812 | 0.921 | 0.987 |
| Wilson | 0.923 | 0.954 | 0.982 |
| Agresti-Coull | 0.931 | 0.958 | 0.989 |
| Andersson-Nerman | 0.928 | 0.951 | 0.978 |

---

## Key Findings

| Interval | Coverage Stability | Simplicity | Recommended Use |
|----------|-------------------|------------|-----------------|
| **Wald** | ❌ Poor | ✅ Simple | Avoid in practice |
| **Wilson** | ✅ Good | ⚠️ Moderate | General use |
| **Agresti-Coull** | ✅ Very Good | ✅ Simple | Teaching & practice |
| **Andersson-Nerman** | ✅ Good | ⚠️ Complex | Research |

### Main Conclusions

1. **The Wald interval should be avoided** in practice despite its ubiquity in textbooks
2. **Wilson and Agresti-Coull** provide substantial improvements with minimal added complexity
3. **The "add 2 successes and 2 failures" rule** is an excellent practical recommendation
4. **Correlation-adjusted pivots** (Andersson-Nerman) offer theoretical insights and competitive performance

---

## 🛠️ Tools & Stack

| Tool | Purpose |
|------|---------|
| **R** | Statistical computing & simulations |
| **ggplot2** | Visualization |
| **LaTeX** | Thesis typesetting |
| **BibTeX** | Bibliography management |


  <br>
  <i>Sapienza Università di Roma • B.Sc. Thesis • A.Y. 2022/2023</i>
</p>
