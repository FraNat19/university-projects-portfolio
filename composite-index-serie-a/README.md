# ⚽ Advanced Composite Indicators (Mini-course Project)

<p align="center">
  <img src="https://img.shields.io/badge/Tools-Python%20%7C%20Excel-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Methods-Normalization%20%7C%20Weighting%20%7C%20Sub--Indices-blue?style=for-the-badge" alt="Methods">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Project-Composite%20Index-orange?style=for-the-badge" alt="Project">
</p>

<p align="center">
  <b>Final project: design and build an advanced composite indicator from scratch</b>
</p>

---

## 📋 Overview

This repository contains the final project developed for a **mini-course focused on advanced composite indicator construction**, where the exam consisted of designing and building a composite index on a topic of choice.

The project delivers a **full composite-index pipeline**:
- Goal definition
- Data collection
- Indicator design
- Normalization
- Weighting
- Aggregation into sub-indices
- Time-series tracking

---

## Project: TCCI (Title Contention Composite Index)

The final project proposes the **TCCI – Title Contention Composite Index**, designed to assess a football team's likelihood of winning the league title by combining performance and contextual dimensions into a single score.

### Objectives

| Goal | Description |
|------|-------------|
| 🏆 **Identify Contenders** | Rank teams by title-winning potential |
| 📈 **Track Changes** | Monitor evolution across seasons |
| ⚠️ **Highlight Under-performance** | Detect gaps between expectations and outcomes |

```
TCCI Score → Higher = More likely title contender
```

---

## 🧾 Data Sources

| Source | Data Type | Coverage |
|--------|-----------|----------|
| **Opta / Optascore** | Match statistics, xG, shots, possession | 2015/16 – 2022/23 |
| **Transfermarkt** | Market values, squad spending, financials | 2015/16 – 2022/23 |

**Time-series perspective**: Index rankings are compared against actual league-table outcomes across 8 seasons.

---

## Methodology

### Pipeline Overview

```
┌─────────────┐    ┌─────────────┐    ┌─────────────┐    ┌─────────────┐
│    RAW      │ → │ NORMALIZED  │ → │ SUB-INDICES │ → │    TCCI     │
│    DATA     │    │ INDICATORS  │    │   (4)       │    │   SCORE    │
└─────────────┘    └─────────────┘    └─────────────┘    └─────────────┘
     ↑                   ↑                  ↑                  ↑
   Opta            AMPI-style          Weighted           Weighted
 Transfermarkt     goalposts           average            average
```

### 1. Normalization: AMPI Method

**AMPI (Adjusted Mazziotta-Pareto Index)** style rescaling using goalposts (min/max bounds) to allow **absolute comparisons over time**.

$$x_{i,t}^{norm} = \frac{x_{i,t} - GP_{min}}{GP_{max} - GP_{min}} \times 100$$

Where:
- $GP_{min}$ = goalpost minimum (fixed reference)
- $GP_{max}$ = goalpost maximum (fixed reference)

**Advantage**: Enables cross-season comparisons (not relative to single-year distribution).

### 2. Index Structure

**TCCI** is a weighted aggregation of **4 sub-indices**:

```
                         ┌──────────────────┐
                         │      TCCI        │
                         │  (Final Score)   │
                         └────────┬─────────┘
                                  │
        ┌─────────────┬───────────┼───────────┬─────────────┐
        │             │           │           │             │
   ┌────▼────┐   ┌────▼────┐ ┌────▼────┐ ┌────▼────┐
   │OFFENSIVE│   │DEFENSIVE│ │HISTORICAL│ │ MARKET │
   │  INDEX  │   │  INDEX  │ │  INDEX  │ │  INDEX │
   │  (0.3)  │   │  (0.3)  │ │  (0.2)  │ │  (0.2) │
   └─────────┘   └─────────┘ └─────────┘ └─────────┘
```

| Sub-Index | Weight | Focus |
|-----------|--------|-------|
| 🔴 **Offensive** | 30% | Attacking performance |
| 🔵 **Defensive** | 30% | Defensive solidity |
| 🟡 **Historical** | 20% | Past achievements |
| 🟢 **Market** | 20% | Financial strength |

### 3. Indicators by Sub-Index

#### 🔴 Offensive Index (0.3)

| Indicator | Description | Polarity |
|-----------|-------------|----------|
| xG (Expected Goals) | Attacking quality | + |
| xG Overperformance | Goals - xG gap | + |
| Shots per match | Volume of attacks | + |
| Possession % | Ball control | + |
| Pass accuracy | Build-up quality | + |

#### 🔵 Defensive Index (0.3)

| Indicator | Description | Polarity |
|-----------|-------------|----------|
| xGA (Expected Goals Against) | Defensive exposure | - |
| Clean sheets | Matches without conceding | + |
| Goals conceded | Actual defensive output | - |
| xGA Overperformance | Conceded - xGA gap | - |

#### 🟡 Historical Index (0.2)

| Indicator | Description | Polarity |
|-----------|-------------|----------|
| League titles (last 10y) | Championship pedigree | + |
| Average league position | Consistency | + |
| Points per season (avg) | Historical strength | + |
| European qualifications | Top-tier finishes | + |

#### 🟢 Market Index (0.2)

| Indicator | Description | Polarity |
|-----------|-------------|----------|
| Squad market value | Overall squad worth | + |
| Net transfer spending | Investment capacity | + |
| Wage bill | Financial commitment | + |
| Revenue | Club financial health | + |

### 4. Aggregation Formula

**Sub-index calculation:**
$$SI_k = \sum_{j=1}^{m_k} w_j \cdot x_j^{norm}$$

**Final TCCI:**
$$TCCI = 0.3 \cdot SI_{off} + 0.3 \cdot SI_{def} + 0.2 \cdot SI_{hist} + 0.2 \cdot SI_{mkt}$$


## 🛠️ Tools & Stack

| Tool | Purpose |
|------|---------|
| **R** | Data processing & analysis |
| **Excel** | Initial data exploration |


---

## 📈 Key Insights

### Methodological Takeaways

1. **Goalposts matter**: Fixed min/max bounds enable meaningful cross-season comparisons
2. **Balance is key**: 60% performance (off+def) vs 40% context (history+market)
3. **xG metrics**: Expected goals provide stable signal vs raw goals (less noise)
4. **Financial dimension**: Market value correlates strongly with title contention

### Validation

| Metric | Value |
|--------|-------|
| Correlation (TCCI vs Final Points) | ~0.85 |
| Correct Top-3 Prediction | 7/8 seasons |
| Title Winner in Top-2 TCCI | 8/8 seasons |

---

## 📚 References

- OECD (2008). *Handbook on Constructing Composite Indicators: Methodology and User Guide*.
- Mazziotta, M. & Pareto, A. (2013). Methods for Constructing Composite Indices. *Rivista Italiana di Economia Demografia e Statistica*.
- Opta Sports: https://www.optasports.com/
- Transfermarkt: https://www.transfermarkt.com/


  <br>
  <i>Sapienza Università di Roma • Advanced Composite Indicators • A.Y. 2024/2025</i>
</p>
