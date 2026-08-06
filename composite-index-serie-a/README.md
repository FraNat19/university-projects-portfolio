<div align="center">

# ⚽ Team Competitiveness Composite Index

**Ranking a decade of Serie A with the Adjusted Mazziotta-Pareto Index**

<img src="https://img.shields.io/badge/R-276DC3?style=for-the-badge&logo=r&logoColor=white" alt="R">
<img src="https://img.shields.io/badge/Compind-AMPI-1a7f37?style=for-the-badge" alt="Compind">
<img src="https://img.shields.io/badge/Sapienza-8b1a1a?style=for-the-badge" alt="Sapienza">

<sub>Composite Indicators · MSc in Statistical Methods and Applications · Sapienza University of Rome</sub>

</div>

---

## 📋 The idea

League tables rank teams by points. Points measure one season of results and nothing else — not the squad behind them, not the finances holding it up, not whether the team has been good for a decade or got lucky in May.

This project builds a **composite index of competitiveness** for twelve Serie A clubs across eight seasons, aggregating four dimensions that a table cannot show.

| Sub-index | What it captures | Indicators |
|---|---|---:|
| 🏆 **Historical** | Sustained standing over time | 7 |
| 🛡️ **Defensive** | Solidity without the ball | 8 |
| ⚡ **Offensive** | Creation and conversion | 6 |
| 💰 **Market** | Financial and structural base | 8 |

Teams: Atalanta, Bologna, Fiorentina, Inter, Juventus, Lazio, Milan, Napoli, Roma, Sassuolo, Torino, Udinese.

---

## 🔧 Why AMPI

A composite index has to solve a problem a simple average ignores: **should a team be allowed to compensate a terrible defence with a brilliant attack?**

Weighted averages say yes. The **Adjusted Mazziotta-Pareto Index** says only partly. It normalises each indicator against fixed goalposts, then subtracts a penalty proportional to the variability across indicators. A team that is uniformly good scores better than a team with the same mean but wild swings between dimensions.

```r
goals1 <- apply(df1[df1$Year == "2016/17", c(3:9)], 2, median)

CI1 <- ci_ampi(df1,
               indic_col = c(3:9),
               gp        = goals1,
               time      = df1$Year,
               polarity  = c("NEG","POS","POS","POS","POS","POS","NEG"),
               penalty   = "POS")
```

Two design choices in there are worth naming.

**Fixed goalposts from a reference season.** The medians of 2016/17 anchor the scale, so scores are comparable *across years* rather than being renormalised each season. Without this, a team could improve while its index falls, simply because the league improved faster.

**Explicit polarity per indicator.** `Position` and `Ncoach` are `NEG`: finishing lower and burning through coaches both hurt. Everything else is `POS`. Getting this vector wrong silently inverts an indicator's meaning, which is the classic way composite indices go wrong without anyone noticing.

---

## 📊 The four sub-indices

<details>
<summary><b>Indicator lists and polarities</b></summary>

<br>

**Historical** — `Position`, `PPG` (points per game), `DR` (goal difference), `CupWins`, `EuropePart`, `Trophies`, `Ncoach`
<sub>NEG, POS, POS, POS, POS, POS, NEG</sub>

**Defensive** — `XGS` (expected goals against), `GS` (goals conceded), `PercPar` (save percentage), `PortInv` (clean sheets), `Amm` (bookings), `Fouls`, `Int` (interceptions), `TackW` (tackles won)
<sub>POS, NEG, POS, POS, NEG, NEG, POS, POS</sub>

**Offensive** — `Bposs` (possession), `XG`, `GF` (goals for), `Shoot`, `PassCompl`, `Dribbling`
<sub>POS, NEG, POS, POS, POS, POS</sub>

**Market** — `DIRITTI` (broadcast rights), `COSTOROSA` (squad cost), `VALKEYPLAYERS`, `DEBIT` (debt), `TRANSF` (transfer balance), `STADIUM`, `OLD` (average age), `MANAGER`
<sub>POS, POS, POS, NEG, POS, POS, NEG, POS</sub>

</details>

Each sub-index gets its own goalposts and penalty term, then a trend plot across seasons, so a club's trajectory is visible rather than just its current standing.

---

## ⚠️ A weighting bug worth fixing

The four sub-indices are combined at line 213:

```r
CI <- CI1$ci_ampi_est*0.3 + CI2$ci_ampi_est*0.3 + CI3$ci_ampi_est*0.3 + CI4$ci_ampi_est*0.2
```

**Those weights sum to 1.1, not 1.**

AMPI produces indices centred on 100 by construction. A weighted sum whose weights total 1.1 therefore centres the final TCCI near 110, on a scale that no longer means what the sub-indices mean. The ordering between clubs is unaffected, since the combination is linear and fixed, but the absolute value cannot be read against the 100 baseline and comparing the composite to any sub-index is misleading.

The intended split was presumably 0.3 / 0.3 / 0.2 / 0.2, giving the three sporting dimensions three quarters of the weight and finance one quarter.

Two smaller issues sit in the badge-plotting block:

```r
img <- png::readPNG(Images)
df  <- cbind(df, Images)              # df does not exist yet
df  <- read.csv(".../CompInd.csv")    # it is created here, four lines later
```

`df` is used before assignment, and `readPNG` receives a 96-element vector where it expects one path. Cosmetic in intent, but they stop the script if it is run top to bottom.

Paths are absolute and local throughout:

```r
sPath = "C:/Users/frank/OneDrive/Desktop/Composite indicators/"
```

---

## 📂 Files

| File | Contents |
|---|---|
| [`HistInd.R`](HistInd.R) | Full pipeline: four AMPI sub-indices, aggregation, trend plots |
| [`CompInd.xlsx`](CompInd.xlsx) | Sporting indicators, 12 clubs across the seasons |
| [`CompInd2.xlsx`](CompInd2.xlsx) | Financial and structural indicators |
| [`composit_index (1).pdf`](composit_index%20\(1\).pdf) | Report with the resulting rankings and trend charts |

**Running it.** Needs R with `Compind`, `dplyr`, `ggplot2`, `reshape2`, `ggimage` and `png`. The script reads `CompInd.csv`; the workbooks here hold the same data in Excel form. The club badge PNGs referenced by `sPath` are not included, so the `geom_image` blocks will not run.

---

## 📖 References

- Mazziotta, M. & Pareto, A. (2016). On a generalized non-compensatory composite index for measuring socio-economic phenomena. *Social Indicators Research*, 127(3), 983–1003.
- Fusco, E., Vidoli, F. & Sahoo, B.K. (2018). Spatial heterogeneity in composite indicators. *Social Indicators Research*, 137(2), 635–658.
- OECD (2008). *Handbook on Constructing Composite Indicators: Methodology and User Guide*.
