# 📊 Data-Driven Decision Making

<p align="center">
  <img src="https://img.shields.io/badge/Tools-Excel%20%7C%20Python-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Course-Decision%20Analysis-blue?style=for-the-badge" alt="Course">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Year-2024%2F2025-orange?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>Multi-Criteria Decision Analysis: From Theory to Real-World Applications</b>
</p>

---

## Overview

This repository contains coursework and implementations developed during the **Data-Driven Decision Making** course at Sapienza University of Rome (A.Y. 2024/2025). The projects span multiple decision analysis methodologies, applying various optimization and ranking models to real-world scenarios such as smartphone selection and strategic business planning.


---

## Course Objectives

The course covers fundamental and advanced techniques for:

1. **Multi-Criteria Decision Analysis (MCDA)** - Evaluating alternatives across multiple conflicting criteria
2. **Preference Modeling** - Capturing and formalizing decision-maker preferences
3. **Strategic Planning Tools** - Implementing frameworks for organizational decision-making
4. **Optimization Models** - Solving multi-objective optimization problems



## Methodologies Implemented

### 1. Multi-Attribute Utility Theory (MAUT)

MAUT provides a systematic approach to evaluate alternatives based on multiple attributes by computing a utility score.

**Application**: Smartphone Ranking (8 devices compared)

```
Utility Score = Σ (weight_i × utility_i)
```

**Key Findings**:
| Ranking Type | Best Choice | Worst Choice |
|--------------|-------------|--------------|
| Personal Preference | iPhone 16 Pro Max | LG W31 |
| MAUT Model | Galaxy A16 5G (430 pts) | LG W31 (211 pts) |

> **Insight**: 18 out of 28 pairwise comparisons were consistent between personal and data-driven rankings, revealing the impact of cognitive biases in decision-making.

---

### 2. Multi-Attribute Attitude (MAA) Model

A scoring model that combines attribute importance weights with performance ratings.

**Attributes Evaluated**:
- 📶 Network capabilities
- 📱 Display quality
- ⚡ Platform performance
- 🔋 Battery life
- 📸 Camera specifications
- 💾 Memory & Storage

**Formula**:
```
Score = Σ (Importance_i × Rating_i)
```

---

### 3. Analytic Hierarchy Process (AHP)

AHP structures complex decisions into a hierarchy of criteria and alternatives, using pairwise comparisons to derive priority weights.

**Process**:
1. Define the goal and criteria hierarchy
2. Perform pairwise comparisons
3. Calculate priority vectors
4. Check consistency ratio (CR < 0.1)
5. Synthesize final rankings

**Consistency Check**:
```
CR = CI / RI

Where:
- CI = (λmax - n) / (n - 1)
- RI = Random Index (tabulated)
```

---

### 4. Balanced Scorecard (BSC)

Strategic management framework applied to **Urban Greening & Smart Technology Integration**.

**Four Perspectives**:

| Perspective | Key Objectives | Targets |
|-------------|----------------|---------|
| 💰 **Financial** | Reduce operation costs, optimize scaling | 51% cost optimization |
| 👥 **Customer** | Improve air quality, green area access | 70% population access |
| ⚙️ **Internal Processes** | Scale IT technologies, improve incident response | 83% IT scalability |
| 📚 **Learning & Growth** | Encourage innovation, AI governance | Continuous R&D |

---

### 5. PROMETHEE

**Preference Ranking Organization METHod for Enrichment Evaluation**

Applied to car selection with preference functions for each criterion.

**Key Concepts**:
- Preference functions (V-shape, U-shape, Gaussian)
- Positive and negative outranking flows
- Net flow for final ranking

---

### 6. Multi-Objective Linear Programming (MOLP)

Optimization problems with multiple conflicting objectives using Clex of IBM.

**Bi-Objective Analysis**:
- Pareto frontier identification
- Trade-off analysis between objectives
- Compromise solutions



## Tools & Software

| Tool | Purpose |
|------|---------|
| **Microsoft Excel** | MAUT, MAA, AHP, BSC, PROMETHEE implementations |
| **Python** | Data processing and analysis |
| **Word** | Documentation and reports |

---
### 👤 Author
**Francesco Natali** 

---

## References

- Saaty, T.L. (1980). *The Analytic Hierarchy Process*. McGraw-Hill.
- Keeney, R.L. & Raiffa, H. (1976). *Decisions with Multiple Objectives*. Wiley.
- Kaplan, R.S. & Norton, D.P. (1996). *The Balanced Scorecard*. Harvard Business Press.
- Brans, J.P. & Vincke, P. (1985). A Preference Ranking Organisation Method. *Management Science*, 31(6).

---

