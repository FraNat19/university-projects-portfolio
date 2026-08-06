<div align="center">

# 🎓 University Projects

**Statistics · Machine learning · Data engineering**

<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
<img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R">
<img src="https://img.shields.io/badge/SQL-4479A1?style=for-the-badge&logo=mysql&logoColor=white" alt="SQL">
<img src="https://img.shields.io/badge/Stata-1a5f7a?style=for-the-badge" alt="Stata">
<img src="https://img.shields.io/badge/MATLAB-0076A8?style=for-the-badge&logo=mathworks&logoColor=white" alt="MATLAB">
<img src="https://img.shields.io/badge/SAS-0766D1?style=for-the-badge" alt="SAS">
<img src="https://img.shields.io/badge/Power%20BI-F2C811?style=for-the-badge&logo=powerbi&logoColor=black" alt="Power BI">

**Francesco Natali** · MSc Statistical Methods and Applications · Sapienza University of Rome · 2022–2025

</div>

---

## ⭐ Start here

Three projects that show the range, if you only have a few minutes.

<table>
<tr>
<td width="33%" valign="top">

### [⚖️ GraphRAG for Safety Law](llm-compliance-agent/)
**Master's thesis**

A retrieval system that pairs vector search with a Neo4j graph of Italian legislation, so it can tell you the law it just cited was repealed in 2008.

<sub>Python · Neo4j · Qdrant · Qwen 2.5 72B · SLURM</sub>

</td>
<td width="33%" valign="top">

### [🧭 Graph Algorithms](algorithms/)
**Measured, not assumed**

Routing on real street networks, minimum-cut hierarchies, and an [interactive dashboard](algorithms/benchmark_dashboard.html) of the benchmarks.

<sub>Python · NetworkX · OSMnx</sub>

</td>
<td width="33%" valign="top">

### [⏳ Discount Rates](empirical-economics-SHIW-data/)
**A framing effect**

Ordered probit on 14,000 Bank of Italy survey responses. The best predictor is not income or age — it is how the question was asked.

<sub>Stata</sub>

</td>
</tr>
</table>

---

## 🤖 AI and machine learning

| Project | What it does | Stack |
|---|---|---|
| [**GraphRAG for Workplace Safety**](llm-compliance-agent/) <sub>MSc thesis</sub> | Hybrid retrieval over 1,000+ INAIL and EU-OSHA documents plus a knowledge graph of 50+ laws and 3,000+ articles. The graph tracks repeals and amendments, so answers built on a superseded statute get flagged instead of quoted. | Python · Neo4j · Qdrant · Docling · Ollama |
| [**Machine Learning & Big Data**](machine-learning/) <sub>2023/24</sub> | Income prediction on 32,561 census records across logistic regression, random forest, gradient boosting and a 16-configuration neural network sweep. Plus computer vision: a CNN trained from scratch on 320 images reaches 28% validation accuracy where a pretrained DenseNet121 reaches 97%. | TensorFlow · scikit-learn · BeautifulSoup |

---

## 📐 Statistical modelling and inference

| Project | What it does | Stack |
|---|---|---|
| [**Spatial Prediction of Shrimp Biomass**](spatial-modeling/) <sub>2024/25</sub> | Kriging over 120 MEDITS trawl stations per year, half of them recording zero catch. Maximum-likelihood against Bayesian, compared on RMSE and interval score. | R · geoR · gstat |
| [**Intertemporal Discount Rates**](empirical-economics-SHIW-data/) <sub>2024</sub> | Ordered probit on the Bank of Italy SHIW 2010 survey. Adding one dummy for question order triples the model's pseudo R². | Stata |
| [**Bachelor Thesis**](Bachelor-thesis/) <sub>2023</sub> | Improved approximate confidence intervals for the binomial proportion, with simulation study. | R · LaTeX |
| [**Multiple Time Series**](multiple-time-series/) | VAR stability and impulse responses, spurious regression, unit roots, cointegration. | R |
| [**Advances in Data Analysis**](advanced-data-analysis/) <sub>2024/25</sub> | Composite crime index for 1,994 US communities via disjoint factor analysis; clustering 41 countries on the Better Life Index. | MATLAB |
| [**Multivariate Analysis**](stat-mult/) <sub>BSc</sub> | Five case studies: PCA, clustering and regression on regional agriculture, exam grades and appliance pricing. | SAS |
| [**Adaptive Web Sampling**](adaptive-web-sampling/) | Sampling designs for hidden and networked populations. | R |

---

## 💻 Computer science and data engineering

| Project | What it does | Stack |
|---|---|---|
| [**Graph Algorithms**](algorithms/) <sub>2024/25</sub> | Four projects: OpenStreetMap routing, social network analysis, empirical complexity benchmarking, and Gomory-Hu cut trees applied to dependency-network fragility. Every implementation verified against a reference before being benchmarked. | Python · NetworkX · OSMnx |
| [**Cardiology Ward Database**](sql-hospital/) <sub>2021/22</sub> | A MySQL schema designed from an interview with a cardiology resident. Twenty-one tables, and a validator that runs the whole thing without a MySQL server. | MySQL · Python |
| [**Country Economic Dashboards**](Country-Economic-Dashboard-Automation/) <sub>2025</sub> | Upload a country CSV to SharePoint, receive an email with a Power BI link already filtered to that country. | Power BI · Power Automate |
| [**ECB Lending Pipeline**](ECB-Data-Automation-Pipeline/) <sub>2025</sub> | The same pattern over 273 monthly observations of euro-area corporate borrowing cost. Logs before it works, so failed runs stay visible. | Power BI · Power Automate |

---

## 🎯 Decision science

| Project | What it does | Stack |
|---|---|---|
| [**Multi-Criteria Decision Analysis**](decision-analysis/) <sub>2024/25</sub> | PROMETHEE II, AHP with consistency ratios, MAUT, Fishbein scoring, topological sorting of incomplete preferences, bi-objective optimisation. The Python reimplementation found an inconsistency in one spreadsheet that changed which alternative wins. | Python · Excel · CPLEX |
| [**Serie A Competitiveness Index**](composite-index-serie-a/) | Four AMPI sub-indices — historical, defensive, offensive, market — across 12 clubs and 8 seasons, with a penalty term limiting compensation between dimensions. | R · Compind |

---

## 🗂️ Reading this repository

Each folder has its own README covering what the project does, what the results were, and what does not run without files that cannot be redistributed.

The results reported are the ones the code actually produced, including the ones that did not work: an overfitted CNN, a refinement that turned out to be unnecessary, a spreadsheet whose two calculations disagreed.

Where a project depends on data I cannot share — course archives, Bank of Italy microdata, third-party corpora — the README says which files are missing and why.

---

## 📬 Contact

**Francesco Natali** · [franknatali01@gmail.com](mailto:franknatali01@gmail.com)
