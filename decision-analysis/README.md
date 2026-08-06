<div align="center">

# ⚖️ Multi-Criteria Decision Analysis

**Six methods for choosing between alternatives when the criteria disagree**

<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
<img src="https://img.shields.io/badge/Excel-217346?style=for-the-badge&logo=microsoftexcel&logoColor=white" alt="Excel">
<img src="https://img.shields.io/badge/CPLEX-0530AD?style=for-the-badge" alt="CPLEX">
<img src="https://img.shields.io/badge/Sapienza-8b1a1a?style=for-the-badge" alt="Sapienza">

<sub>Data-Driven Decision Making · MSc · Sapienza University of Rome · A.Y. 2024/2025</sub>

</div>

---

## 📋 Overview

Every method here answers the same question in a different way: given several alternatives and several criteria that point in different directions, which do you pick, and how much do you trust the answer?

The models were built in Excel during the course. [`mcda.py`](mcda.py) reimplements three of them in Python so the spreadsheets can be **checked rather than trusted** — and one of the checks found a real error.

```bash
python mcda.py
```

| Method | Decision problem | Files |
|---|---|---|
| **PROMETHEE II** | Rank 5 smartphones on price, weight, camera, screen | [`Template Promethee (cars).xlsx`](Template%20Promethee%20\(cars\).xlsx) |
| **AHP** | Where to put campus parking | [`AHP_model_Naali_1945581 (1).xlsx`](AHP_model_Naali_1945581%20\(1\).xlsx), `AHPNatali.docx` |
| **MAUT** | Rank 8 smartphones, compared against gut feeling | `MAUT_NATALI_1945581.docx` |
| **Fishbein MAA** | Score the same 8 phones on 9 weighted attributes | [`multi-attribute-attitude-model.xlsx`](multi-attribute-attitude-model.xlsx) |
| **Topological sorting** | Complete an incomplete preference order over 7 EVs | `DM_preferences_Natali_1945581.docx` |
| **Bi-objective MOLP** | Knapsack trading safety against value, ε-constraint | `Biobj_1945581.docx`, `MolpNatali (1).pdf` |
| **Balanced Scorecard** | Strategy maps for a brewery and for urban greening | `BSC_Group_AGATE_NATALI.xlsx`, `BSC_AGATENATALI.docx` |

---

## 🔎 PROMETHEE II, and an error in the spreadsheet

Five phones, four criteria, linear preference functions with indifference threshold *q* and preference threshold *p*.

| Criterion | Direction | Weight | q | p |
|---|:-:|---:|---:|---:|
| Price | min | 0.35 | 100 | 250 |
| Weight | min | 0.15 | 4 | 9 |
| Videocamera | max | 0.25 | 1 | 2 |
| Screen | max | 0.25 | 0 | 0.04 |

The net flow can be computed two ways, and in correct PROMETHEE II they must agree: as φ⁺ − φ⁻ from the aggregated preference matrix, or as the weighted sum of the unicriterion net flows. **The spreadsheet reports both, and they disagree.**

| Phone | Excel: φ⁺ − φ⁻ | Excel: weighted sum | Python |
|---|---:|---:|---:|
| Xiaomi Redmi A3 Pro | 0.1415 | 0.0228 | **0.1415** |
| Honor 200 Smart | 0.1375 | 0.2625 | **0.1375** |
| iPhone 16 Pro Max | 0.0938 | 0.1313 | **0.0938** |
| Galaxy A16 5G | 0.0022 | −0.0227 | **0.0022** |
| Huawei Pura 70 Ultra | −0.3750 | −0.5250 | **−0.3750** |

The Python implementation computes both routes and confirms they coincide to floating-point precision:

```
max |phi+ - phi-  minus  weighted unicriterion sum| = 5.55e-17
the two routes to the net flow agree
```

It reproduces the φ⁺ − φ⁻ column exactly, to four decimals on all five alternatives. **The second table in the spreadsheet is the faulty one**, and it matters: it puts Honor 200 Smart first, where the correct computation ranks Xiaomi Redmi A3 Pro ahead of it.

The winner is a phone that loses badly on camera (−1.00, the worst possible unicriterion flow) and wins on price. With price weighted at 0.35 against the camera's 0.25, that trade is exactly what the weights were asking for.

---

## 🏛️ AHP with consistency checking

A two-level hierarchy for a campus parking decision between three alternatives: a large area outside the campus, a garage, and a small area inside.

```
GOAL
├── CUSTOMER (0.581)
│   ├── SERVICE      → 0.3382 global
│   ├── COST         → 0.1797
│   └── CONVENIENCE  → 0.0637
├── CONTROL          → 0.3090
└── INVESTMENT       → 0.1095
```

Priorities come from the geometric mean of each pairwise comparison matrix, and every matrix is checked before its weights are used:

```python
def consistency(matrix):
    m = np.asarray(matrix, dtype=float)
    n = m.shape[0]
    w = priority_vector(m)
    lam = float((m @ w / w).mean())
    ci = (lam - n) / (n - 1)
    cr = ci / RANDOM_INDEX[n]
    return {"lambda_max": lam, "CI": ci, "CR": cr, "acceptable": cr < 0.10}
```

| Matrix | λ_max | CI | CR | Acceptable |
|---|---:|---:|---:|:-:|
| GOAL | 3.0037 | 0.0018 | 0.0032 | ✔ |
| CUSTOMER | 3.0037 | 0.0018 | 0.0032 | ✔ |
| CONTROL | 3.0649 | 0.0324 | 0.0559 | ✔ |
| SERVICE, COST, CONVENIENCE, INVESTMENT | 3.0000 | 0.0000 | 0.0000 | ✔ |

Every ratio sits below the 0.10 threshold, so the judgements hold together. `CONTROL` is the loosest, which makes sense: it is the only matrix where a 7 and a 5 appear together, and strong judgements are where inconsistency creeps in.

**Result:** Garage 0.5938, Small area inside 0.2067, Large area outside 0.1995. The spreadsheet gives 0.5898 / 0.2055 / 0.1988; the 0.004 gap comes from the rounded intermediate priorities stored in the sheet.

The garage wins because it dominates the two heaviest criteria at once, taking 0.73 on `CONTROL` and 0.71 on `INVESTMENT`. It is not close.

---

## 📱 MAUT versus intuition

The same eight phones ranked twice: once by gut feeling before looking at specifications, once by the MAUT model.

| | Personal ranking | MAUT model |
|---|---|---|
| Best | iPhone 16 Pro Max | **Galaxy A16 5G** (430) |
| Second | Galaxy A16 5G | iPhone 16 Pro Max (429) |
| Worst | LG W31 | LG W31 (211) |

Across all C(8,2) = 28 pairwise comparisons, **10 disagree — a 35.7% error rate.** More than a third of the intuitive preferences do not survive contact with the attribute data.

The top two are separated by a single point out of 430, which is the more honest reading of the result: the model does not really prefer the Galaxy, it declares a tie that the ranking is then forced to break.

### Fishbein multi-attribute attitude model

A separate scoring exercise on the same phones, nine attributes weighted 1–9 by importance, score = Σ (importance × rating):

| Rank | Phone | Score |
|---:|---|---:|
| 1 | iPhone 16 Pro Max | 380 |
| 2 | Galaxy A16 5G | 374 |
| 3 | Oppo K12 Plus | 355 |
| 4 | Honor 200 Smart | 354 |
| … | | |
| 8 | LG W31 | 215 |

Two models, two different winners, same loser. The Fishbein weights put price at 9 and display at 8, which is enough to flip the top of the ranking relative to MAUT while leaving the bottom untouched.

---

## 🔀 Completing an incomplete preference order

A decision maker gives preferences over seven electric vehicles, but only some pairs:

```
Tesla Model 3     ≻ Nissan Leaf, Ford Mustang Mach-E, Hyundai Kona Electric
Polestar 2        ≻ Audi e-tron, Volkswagen ID.4
Audi e-tron       ≻ Nissan Leaf
Ford Mustang Mach-E ≻ Hyundai Kona Electric
Volkswagen ID.4   ≻ Hyundai Kona Electric
```

These arcs form a directed acyclic graph, and every topological ordering of that graph is a complete ranking consistent with what the decision maker actually said. Enumerating them shows how much freedom the missing comparisons leave, and which questions would need answering to narrow it down.

---

## 🎯 Bi-objective optimisation

A 0-1 knapsack with two competing objectives, safety **S** and total value **T**, solved under three weightings and analysed with the ε-constraint method to trace the Pareto frontier.

| Scenario | Weights | Safety | Total value | Investments selected |
|---|---|---:|---:|---|
| 1 | Original | 133 | **154** | 1, 3, 5, 7, 9, 14 |
| 2 | Increased | **286** | 145 | 3, 7, 8, 14 |
| 3 | Decreased | 66 | 76 | 1, 3, 5, 7, 9, 14 |

Scenario 2 is the interesting one. Raising both weights more than doubles safety, from 133 to 286, at a cost of nine points of value, and it does so by holding only four investments instead of six. Concentration buys safety here; diversification buys value.

Scenario 3 selects the same portfolio as Scenario 1 at lower objective values, which is the expected behaviour when weights scale without changing their ratio.

---

## 📈 Balanced Scorecard

Two distinct exercises.

**Built from scratch:** a strategy map for a **brewery**, across the four standard perspectives. Objectives are paired with measurable targets rather than intentions, which is the part most scorecards get wrong.

| Perspective | Example objective | Target |
|---|---|---|
| Customer | Introduce a takeaway service | 10% of total sales from takeaway |
| Internal process | Data-driven loyalty programme | 20% of revenue from loyalty members |
| Economic | Branded merchandise line | Additional revenue stream, brand promotion |
| Learning & growth | Test new beers on a small scale first | Sales per batch before wider launch |

**Analysed, not built:** a scorecard for sustainable **urban greening through smart technology**, reviewed across its four pillars — 51% cost optimisation and 75% unit-level cost-efficiency on the financial side, air quality and green-area access on the customer side, IT scalability and incident response internally, and AI governance under learning and growth.

---

## ⚙️ Running the code

```bash
pip install numpy pandas
python mcda.py
```

`mcda.py` implements PROMETHEE II with linear preference functions, AHP over an arbitrary two-level hierarchy with consistency ratios, and the Fishbein model. Each run prints its results next to the spreadsheet figures so any divergence is visible immediately.

---

## 📖 References

- Saaty, T.L. (1980). *The Analytic Hierarchy Process*. McGraw-Hill.
- Keeney, R.L. & Raiffa, H. (1976). *Decisions with Multiple Objectives*. Wiley.
- Brans, J.P. & Vincke, P. (1985). A preference ranking organisation method. *Management Science*, 31(6), 647–656.
- Fishbein, M. (1963). An investigation of the relationships between beliefs about an object and the attitude toward that object. *Human Relations*, 16(3), 233–239.
- Kaplan, R.S. & Norton, D.P. (1996). *The Balanced Scorecard*. Harvard Business Press.
- Haimes, Y.Y., Lasdon, L.S. & Wismer, D.A. (1971). On a bicriterion formulation of the problems of integrated system identification. *IEEE Transactions on Systems, Man, and Cybernetics*, 1(3), 296–297.
