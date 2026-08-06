<div align="center">

# 📐 Multivariate Analysis with SAS

**Five case studies: description, dimensionality reduction, clustering, regression**

<img src="https://img.shields.io/badge/SAS-1a5f7a?style=for-the-badge" alt="SAS">
<img src="https://img.shields.io/badge/PCA%20·%20Clustering%20·%20Regression-1a7f37?style=for-the-badge" alt="Methods">
<img src="https://img.shields.io/badge/Sapienza-8b1a1a?style=for-the-badge" alt="Sapienza">

<sub>Multivariate Data Analysis · BSc · Sapienza University of Rome</sub>

</div>

---

## 📋 Overview

Five independent case studies, each following the same arc: describe the data honestly, find the linear structure inside it, reduce it to something interpretable, then group the units that behave alike.

The datasets are deliberately unglamorous — fruit harvests, washing machines, exam grades — which is the point. The techniques are the same ones used on data that matters, and small tidy datasets make it obvious when a method is being misapplied.

| Exercise | File | Units × variables | Question |
|---|---|---|---|
| 1 · Fruit | [`frutta.pdf`](frutta.pdf) | 20 Italian regions × 7 crops | Which regions specialise in what, and which behave alike? |
| 3 · Students | [`student.pdf`](student.pdf) | 88 students × 5 exams | Do open-book exams produce different outcomes from closed-book? |
| 5 · Washing machines | [`lavatrici.pdf`](lavatrici.pdf) | 16 models × 6 features | What actually drives the price? |
| · Diet | [`dieta.pdf`](dieta.pdf) | European countries × 10 foods | What dietary patterns exist across Europe? |
| · Cereals, Labour force | [`Cereali.docx`](Cereali.docx), [`FORZA LAVORO.docx`](FORZA%20LAVORO.docx) | — | Further case studies |

---

## 🍇 Exercise 1 — Fruit production across Italian regions

Seven woody crops measured across all 20 regions: table grapes, wine grapes, wine produced, oranges, mandarins, clementines, lemons.

| Variable | Mean | Median | CV |
|---|---:|---:|---:|
| Wine grapes | 3,677.35 | 2,364.50 | 94.48 |
| Wine produced | 2,746.35 | 1,618.50 | 97.10 |
| Oranges | 924.25 | 0.00 | 289.32 |
| Table grapes | 539.80 | 20.00 | 261.13 |
| Lemons | 331.70 | 0.00 | 382.56 |

**Every mean sits far above its median, and four of the seven medians are zero.** More than half of Italy's regions produce no citrus at all — Lombardy, Piedmont, Valle d'Aosta, Trentino Alto Adige and Veneto grow none — so the distributions are not skewed by outliers so much as split into producers and non-producers.

That structure is the finding. Citrus and table grapes separate the south; wine grapes are the one crop grown, at varying intensity, everywhere.

---

## 🎓 Exercise 3 — Do open-book exams change outcomes?

Grades out of 100 for 88 students across five exams. **Rational Mechanics and Vector Algebra are closed-book; Algebra, Calculus and Statistics are open-book.**

| Exam | Format | Mean | Median | SD |
|---|---|---:|---:|---:|
| Algebra | open | 50.60 | 50.00 | 10.56 |
| Vector Algebra | closed | 50.59 | 51.00 | 13.07 |
| Calculus | open | 46.68 | 49.00 | 14.76 |
| Statistics | open | 42.31 | 40.00 | 17.16 |
| Rational Mechanics | closed | 38.95 | 41.50 | 17.39 |

The design invites the conclusion that open-book exams are easier. The data does not support it cleanly: the two highest means are one open-book and one closed-book exam, separated by 0.01 points. What separates the exams is **dispersion**, not level. Algebra has a standard deviation of 10.56 while Mechanics and Statistics both exceed 17, so the open-book format is not raising grades so much as the individual exams differ in how sharply they sort students.

---

## 🔌 Exercise 5 — What drives a washing machine's price?

Sixteen models, six characteristics: price, spin speed, water consumption, depth, energy use, load capacity.

| Variable | Mean | Median | CV | Range |
|---|---:|---:|---:|---:|
| Price ($) | 624.88 | 568.50 | 45.98 | 309–1,469 |
| Spin (rpm) | 1,065.63 | 1,000.00 | 31.77 | 600–1,600 |
| Water (l) | 54.25 | 53.50 | 18.84 | 39–75 |
| Depth (cm) | 56.19 | 55.00 | 5.49 | 51–60 |
| Energy (kWh) | 102.13 | 95.00 | 11.42 | 94–133 |
| Load (kg) | 53.44 | 50.00 | 10.83 | 50–70 |

**Price varies eight times more than depth.** The coefficient of variation ranges from 5.49 on depth to 45.98 on price, which says machines that are nearly identical physically are sold at wildly different prices. Depth is standardised by the cabinets they slot into; price is not standardised by anything.

The median load equals the minimum, so more than half the models sit at exactly 50 kg — another variable carrying almost no information for a regression.

---

## 🥗 Diet patterns across Europe

Ten food categories per country: cereals, rice, potatoes, sugar, vegetables, wine, meat, milk, butter, eggs.

| Food | Mean | Median | CV |
|---|---:|---:|---:|
| Milk | 128.4 L | 125.3 L | 35.8 |
| Vegetables | 95.6 kg | 82.5 kg | 57.5 |
| Wine | 25.9 L | 21.5 L | 75.9 |
| Rice | 4.2 kg | 4.3 kg | 28.9 |

**PCA retains four components explaining 83.7% of variance:**

| Component | Interpretation | Dominant loadings |
|---|---|---|
| PC1 | Mediterranean diet | +rice, +vegetables, +wine, −sugar, −milk |
| PC2 | Protein-rich diet | +meat, +eggs, +butter, −milk, −rice |
| PC3 | Grain-based diet | +cereals, +vegetables, −butter |
| PC4 | Starch | +potatoes |

PC1 is the one worth reading: it does not describe how much people eat, it describes a *trade-off*. Rice, vegetables and wine load positively while sugar and milk load negatively, so the component separates southern from northern European eating rather than ranking countries on quantity.

---

## 🧭 Method notes

**Choosing how many components to keep** — four criteria applied together rather than any one alone:

| Criterion | Rule |
|---|---|
| Kaiser | Eigenvalue > 1 |
| Scree plot | Elbow in the variance curve |
| Cumulative variance | 70–80% explained |
| Interpretability | The components have to mean something |

The last one is doing real work here. A component that explains 8% of variance and loads on three unrelated variables is noise with a number attached, and no amount of eigenvalue arithmetic makes it a construct.

**Reading the coefficient of variation first** is the habit that runs through all five exercises. Comparing standard deviations across variables measured in litres, kilograms, centimetres and dollars means nothing; the CV makes them comparable, and it is what exposes the zero-inflation in the fruit data and the standardised depth in the appliance data before any model is fitted.

---

<sub>Written in SAS. Each PDF contains the full output: descriptive tables, correlation matrices, component loadings, dendrograms and regression diagnostics.</sub>
