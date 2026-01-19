# Sample Survey Project: Adaptive Web Sampling (AWS)

<p align="center">
  <img src="https://img.shields.io/badge/Tools-R%20%7C%20Python-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Methods-Adaptive%20Sampling%20%7C%20MCMC%20%7C%20Rao--Blackwell-blue?style=for-the-badge" alt="Methods">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Year-2024%2F2025-orange?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>From design definition to simulation study (without focusing on results)</b>
</p>

---

## Overview

This repository contains the exam project for the **Sample Survey** course, developed as a group project on **Adaptive Web Sampling (AWS)**—a flexible class of adaptive link-tracing designs for network and spatial populations.

The work covers the full pipeline:
- Problem framing
- Formal design specification
- Estimator definitions
- Simulation study to compare design/estimation choices


---

## Project Goals

| Objective | Description |
|-----------|-------------|
| **Present AWS** | Alternative to classic link-tracing approaches (snowball sampling, adaptive cluster sampling) |
| **Flexibility** | Highlight AWS's control over sample size and exploration vs link-following balance |
| **Describe Strategy 1** | Initial sample, evolving active set, mixture-based transitions with dampening parameter $d$ |
| **Simulation Framework** | Study estimator behavior, variance reduction via Rao–Blackwellization, MCMC approximation |

---

## Concepts Covered

### Network & Spatial Sampling Setup

In AWS, the population consists of units (nodes) with:
- A variable of interest $y_i$
- Link/adjacency information $w_{ij}$ observed only through sampling

```
Population U = {1, 2, ..., N}

For each unit i:
  • y_i = variable of interest
  • N_i = {j : w_ij = 1} = neighborhood (links)
  
Links are revealed only when a unit is sampled.
```

### AWS Design Mechanics

**Key Innovation**: Mixture of link-tracing and random jumps controlled by dampening parameter $d$.

```
┌─────────────────────────────────────────────────────────────┐
│                    AWS Strategy 1                           │
├─────────────────────────────────────────────────────────────┤
│                                                             │
│  Step 0: Draw initial sample S₀ (size n₀) via SRS/other    │
│          Initialize Active Set A₀ = S₀                      │
│                                                             │
│  Step t: For t = n₀+1, ..., n:                              │
│                                                             │
│          With probability (1-d):                            │
│            → Select from neighbors of Active Set            │
│            → (Link-tracing / exploration)                   │
│                                                             │
│          With probability d:                                │
│            → Random jump to any unsampled unit              │
│            → (Prevents trapping in clusters)                │
│                                                             │
│          Update Active Set Aₜ                               │
│                                                             │
└─────────────────────────────────────────────────────────────┘
```

**Dampening Parameter $d$**:

| Value of $d$ | Behavior |
|--------------|----------|
| $d = 0$ | Pure link-tracing (like snowball) |
| $d = 1$ | Pure random sampling (no adaptation) |
| $0 < d < 1$ | Hybrid: balances exploration and coverage |

### Transition Mechanism

At each step, the selection probability is a **mixture**:

$$P(\text{select } j \mid S_{t-1}, A_{t-1}) = (1-d) \cdot \pi_j^{(link)} + d \cdot \pi_j^{(random)}$$

Where:
- $\pi_j^{(link)}$ = probability based on links from active set
- $\pi_j^{(random)}$ = uniform probability over unsampled units

---

## 📐 Estimation in AWS

### Challenge

AWS is an **adaptive design** → standard Horvitz-Thompson estimators need careful handling of inclusion probabilities.

### Estimators Implemented

| Estimator | Description |
|-----------|-------------|
| **Horvitz-Thompson (HT)** | $\hat{\mu}_{HT} = \frac{1}{N} \sum_{i \in S} \frac{y_i}{\pi_i}$ |
| **Hájek (Ratio)** | $\hat{\mu}_{H} = \frac{\sum_{i \in S} y_i / \pi_i}{\sum_{i \in S} 1 / \pi_i}$ |
| **Rao-Blackwell (RB)** | Conditional expectation given reduced data |

### Rao-Blackwellization

**Key Idea**: Improve estimators by conditioning on **minimal sufficient statistics** (the reduced data).

$$\hat{\mu}_{RB} = E[\hat{\mu} \mid D_{red}]$$

Where $D_{red}$ = reduced data (what was observed, not the order of selection).

**Variance Reduction**:
$$Var(\hat{\mu}_{RB}) \leq Var(\hat{\mu})$$

---

## 🔄 MCMC for Rao-Blackwellization

Computing $E[\hat{\mu} \mid D_{red}]$ exactly is often intractable → use **MCMC approximation**.

### Metropolis-Hastings Resampling

```
Algorithm: MCMC Rao-Blackwellization
─────────────────────────────────────
Input: Observed sample S, reduced data D_red, number of iterations M

1. Initialize: current path π⁽⁰⁾ = observed selection order

2. For m = 1, ..., M:
   
   a. Propose: π' = permute(π⁽ᵐ⁻¹⁾)  [swap two adjacent selections]
   
   b. Check validity: Is π' consistent with D_red?
      - If NO: reject, π⁽ᵐ⁾ = π⁽ᵐ⁻¹⁾
      - If YES: continue to (c)
   
   c. Compute acceptance ratio:
      α = min(1, P(π') / P(π⁽ᵐ⁻¹⁾))
   
   d. Accept/Reject:
      - With prob α: π⁽ᵐ⁾ = π'
      - Otherwise: π⁽ᵐ⁾ = π⁽ᵐ⁻¹⁾
   
   e. Compute θ̂⁽ᵐ⁾ = estimator under path π⁽ᵐ⁾

3. Return: θ̂_RB ≈ (1/M) Σ θ̂⁽ᵐ⁾
```

---

## 🧪 Simulation Study (High Level)

### Objective

Compare estimator performance under different AWS configurations.

### Design Parameters

| Parameter | Symbol | Values Tested |
|-----------|--------|---------------|
| Population size | $N$ | Fixed |
| Initial sample size | $n_0$ | 10, 20, 30 |
| Final sample size | $n$ | 50, 100, 150 |
| Dampening parameter | $d$ | 0.1, 0.3, 0.5, 0.7 |
| Replications | $R$ | 1000 |

### Simulation Workflow

```
┌──────────────────┐
│ Generate/Load    │
│ Population       │
└────────┬─────────┘
         ↓
┌──────────────────┐
│ For r = 1 to R:  │
│  • Draw AWS      │
│    sample        │
│  • Compute all   │
│    estimators    │
│  • Store results │
└────────┬─────────┘
         ↓
┌──────────────────┐
│ Aggregate:       │
│  • Bias          │
│  • Variance      │
│  • MSE           │
│  • Coverage      │
└────────┬─────────┘
         ↓
┌──────────────────┐
│ Compare across   │
│ (n₀, n, d)       │
│ configurations   │
└──────────────────┘
```

### Performance Metrics

| Metric | Formula |
|--------|---------|
| Bias | $\frac{1}{R} \sum_{r=1}^{R} (\hat{\mu}^{(r)} - \mu)$ |
| Variance | $\frac{1}{R-1} \sum_{r=1}^{R} (\hat{\mu}^{(r)} - \bar{\hat{\mu}})^2$ |
| MSE | $Bias^2 + Variance$ |
| Coverage | $\frac{1}{R} \sum_{r=1}^{R} \mathbf{1}[\mu \in CI^{(r)}]$ |

### Population Types

| Type | Description |
|------|-------------|
| **Network** | Social network structure with clustering |
| **Spatial** | Lattice with spatial autocorrelation |

---

## 📚 Key References

- **Thompson, S.K.** (2006). Adaptive Web Sampling. *Biometrics*, 62(4), 1224-1234.
- **Thompson, S.K. & Seber, G.A.F.** (1996). *Adaptive Sampling*. Wiley.
- **Rao, J.N.K.** (1965). On two simple schemes of unequal probability sampling without replacement. *JASA*.
- **Särndal, C.E., Swensson, B., & Wretman, J.** (1992). *Model Assisted Survey Sampling*. Springer.

---

## 📊 Comparison with Other Designs

| Design | Sample Size Control | Link-Tracing | Random Jumps |
|--------|---------------------|--------------|--------------|
| SRS | ✅ Fixed | ❌ | ✅ |
| Snowball | ❌ Variable | ✅ | ❌ |
| Adaptive Cluster | ❌ Variable | ✅ | ❌ |
| **AWS** | ✅ Fixed | ✅ | ✅ |

**AWS Advantage**: Combines benefits of adaptive designs (efficient for clustered populations) with sample size control (important for budget/logistics).

  <i>Sapienza Università di Roma • Sample Survey • A.Y. 2024/2025</i>
</p>
