<div align="center">

# 🎓 University Projects

**Statistics, machine learning and data engineering · 2022–2025**

<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
<img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R">
<img src="https://img.shields.io/badge/SQL-4479A1?style=for-the-badge&logo=mysql&logoColor=white" alt="SQL">
<img src="https://img.shields.io/badge/Stata-1a5f7a?style=for-the-badge" alt="Stata">
<img src="https://img.shields.io/badge/MATLAB-0076A8?style=for-the-badge&logo=mathworks&logoColor=white" alt="MATLAB">
<img src="https://img.shields.io/badge/Power%20BI-F2C811?style=for-the-badge&logo=powerbi&logoColor=black" alt="Power BI">

**Francesco Natali** · MSc Statistical Methods and Applications · Sapienza University of Rome

</div>

---

Work from my Bachelor's and Master's degrees, collected here with the code, the data where it can be shared, and the reports. Each folder has its own README explaining what the project does, what the results were, and what does not run without files I cannot redistribute.

---

## 🚀 Main projects

### [LLM Agent for Workplace Safety Compliance](llm-compliance-agent/) · 2025
Master's thesis. A GraphRAG system that answers questions on Italian workplace safety law by combining vector search over 2,000+ INAIL and EU-OSHA documents with a Neo4j knowledge graph of the legislation itself. The graph is what makes it useful: it knows which laws have been repealed, so the system flags an answer built on a regulation that is no longer in force instead of quoting it as current.
<sub>Python · Docling · Qdrant · Neo4j · Qwen 2.5 72B</sub>

### [Spatial Prediction of Shrimp Biomass](spatial-modeling/) · 2024/25
Kriging models predicting deep-water rose shrimp biomass along the Tyrrhenian coast from 120 MEDITS trawl stations per survey year, half of which record zero catch. Maximum-likelihood and Bayesian approaches compared on RMSE and interval score across 2002 and 2008.
<sub>R · geoR · gstat · spdep</sub>

### [Graph Algorithms](algorithms/) · 2024/25
Four projects: routing on real OpenStreetMap street networks, social network analysis, empirical complexity benchmarking, and minimum-cut hierarchies. Includes an [interactive benchmark dashboard](algorithms/benchmark_dashboard.html) and extensions written afterwards to explain why A\* came out slower than Dijkstra on a small graph. It expands 3–4× fewer nodes at every size; the heuristic just costs more than it saves until the graph grows.
<sub>Python · NetworkX · OSMnx · Folium</sub>

### [Machine Learning & Big Data Analytics](machine-learning/) · 2023/24
Income prediction on 32,561 census records comparing logistic regression, random forest, gradient boosting and a tuned neural network, plus computer vision and web scraping. A CNN trained from scratch on 320 images reaches 28% validation accuracy where a pretrained DenseNet121 reaches 97%, which is the whole argument for transfer learning in one comparison.
<sub>Python · TensorFlow · scikit-learn · BeautifulSoup · Selenium</sub>

### [Multi-Criteria Decision Analysis](decision-analysis/) · 2024/25
PROMETHEE II, AHP with consistency checking, MAUT, Fishbein scoring, topological sorting of incomplete preferences, and bi-objective optimisation with the ε-constraint method. The Python reimplementation found an inconsistency in one of the spreadsheets that changed which alternative wins.
<sub>Python · Excel · CPLEX</sub>

### [Cardiology Ward Database](sql-hospital/) · 2021/22
A MySQL database designed from a structured interview with a cardiology resident at Policlinico Gemelli. Twenty-one tables covering the clinical path from symptom to therapy and the documentary path from exam to archive, with a validation script that runs the whole schema and every query without a MySQL server.
<sub>MySQL · Python</sub>

### [Intertemporal Discount Rates](empirical-economics-SHIW-data/) · 2024
Ordered probit on the Bank of Italy SHIW 2010 survey, estimating how steeply people discount the future. The strongest predictor turns out not to be income, age or education but **the order in which the question was asked** — adding that single dummy triples the model's pseudo R².
<sub>Stata</sub>

### [Bachelor Thesis: Confidence Intervals for the Binomial](Bachelor-thesis/) · 2023
Study and proposal of improved approximate confidence intervals for the binomial proportion, with simulations in R and the full thesis in LaTeX.
<sub>R · LaTeX</sub>

---

## 📊 Data pipelines and dashboards

### [Country Economic Dashboard Automation](Country-Economic-Dashboard-Automation/) · 2025
Upload `Italy_Economy.csv` to SharePoint and an email arrives with a Power BI link already filtered to Italy. World Bank indicators for 6 countries over 25 years, with productivity measures derived against working-age population rather than total population.
<sub>Power BI · Power Automate · SharePoint</sub>

### [ECB Corporate Lending Pipeline](ECB-Data-Automation-Pipeline/) · 2025
The same pattern applied to 273 monthly observations of euro-area corporate borrowing cost, 2003 to 2024. Logs before it works rather than after, so failed runs stay visible.
<sub>Power BI · Power Automate · SharePoint</sub>

---

## 📚 Statistical modelling

| Project | What it does | Tools |
|---|---|---|
| [Advances in Data Analysis](advanced-data-analysis/) | Composite crime index for 1,994 US communities via disjoint factor analysis; clustering 41 countries on the Better Life Index | MATLAB |
| [Serie A Competitiveness Index](composite-index-serie-a/) | Four AMPI sub-indices — historical, defensive, offensive, market — over 12 clubs and 8 seasons | R · Compind |
| [Multiple Time Series](multiple-time-series/) | VAR stability and impulse responses, spurious regression, unit roots, cointegration | R |
| [Multivariate Statistics](stat-mult/) | PCA, factor analysis, clustering and regression case studies | SAS |
| [Adaptive Web Sampling](adaptive-web-sampling/) | Sampling designs for hidden and networked populations | R |

---

## 🗂️ How to read this repository

Each folder contains its code, whatever data can be redistributed, and a report or slide deck. The READMEs state the results the code actually produced, including the ones that did not work: an overfitted CNN, a refinement that turned out to be unnecessary, a spreadsheet whose two calculations disagreed.

Where a project needs data I cannot share — course archives, Bank of Italy microdata, image datasets — the README says so and lists exactly which files are missing.

---

## 📬 Contact

**Francesco Natali** · [franknatali01@gmail.com](mailto:franknatali01@gmail.com)
