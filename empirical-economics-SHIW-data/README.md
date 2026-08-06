<div align="center">

# ⏳ How Much Is a Year of Waiting Worth?

**Intertemporal discount rates, and a framing effect that beats every demographic**

<img src="https://img.shields.io/badge/Stata-1a5f7a?style=for-the-badge" alt="Stata">
<img src="https://img.shields.io/badge/SHIW%202010-Banca%20d'Italia-8b1a1a?style=for-the-badge" alt="SHIW">
<img src="https://img.shields.io/badge/Ordered%20probit-survey%20weighted-1a7f37?style=for-the-badge" alt="Method">

<sub>Empirical Economics · MSc in Statistical Methods and Applications · Sapienza University of Rome · June 2024</sub>

</div>

---

## 📋 The question

You have just won a lottery prize equal to a full year of your income. You can have it in twelve months, or you can have it right now if you give up part of it. How much would you give up?

The Bank of Italy asked exactly that in the 2010 **Survey on Household Income and Wealth**, question C33. The answer measures an individual's *intertemporal discount rate*: how steeply they discount the future.

The elicitation is a bracketed sequence with a skip pattern, and the survey randomises which end it starts from:

| Group | First offer | Then |
|---|---|---|
| Half the sample | Give up **10%**? | If yes → 20%. If no → 5%, then 2% |
| Other half | Give up **20%**? | If no → 10%, then 5%, then 2% |

Whichever discount they accept becomes their `SCONTO`, on a five-point scale: 0, 2, 5, 10, 20.

```stata
gen SCONTO = .
replace SCONTO = 10 if sconto1 == 1
replace SCONTO = 20 if sconto1 == 1 & sconto2 == 1
replace SCONTO = 5  if sconto1 == 2 & sconto3 == 1
replace SCONTO = 2  if sconto1 == 2 & sconto2 == 2 & sconto3 == 2 & sconto4 == 1
replace SCONTO = 0  if sconto1 == 2 & sconto2 == 2 & sconto3 == 2 & sconto4 == 2
```

That randomisation is the whole point of what follows.

---

## 📊 What people answered

**N = 13,996** household heads, merged across seven SHIW files and weighted by the survey design.

| Discount given up | Share |
|---:|---:|
| 0% — waits for the full amount | 20.32% |
| 2% | 12.09% |
| 5% | 23.27% |
| 10% | 23.23% |
| 20% — most impatient | 21.09% |

A fifth of respondents refuse to give up anything at all. Another fifth give up a full 20%. The middle is thin at 2%, which the report reads carefully: it may reflect genuine preferences, or simply that respondents who accept a discount stop the sequence early and never reach the 2% question.

```stata
svyset nquest [pweight=pesofit]
```

Everything downstream uses the survey design, so estimates are population-representative rather than sample-representative.

---

## 🎯 The result

An ordered probit with sixteen covariates — age, sex, education, region, work status, marital status, risk attitude, income, savings, precautionary reserves, ability to reach the end of the month, debt, property, savings plan — explains a **pseudo R² of 0.0380**.

Then one dummy is added. `id_group` records nothing about the respondent. It records only **which version of the question they happened to be asked**.

| Model | Covariates | LR χ² | Pseudo R² |
|---|---|---:|---:|
| Respondent characteristics | 32 | 239.66 | 0.0380 |
| **Plus the question order** | 33 | **730.45** | **0.1158** |

| Variable | Coefficient | z | P>&#124;z&#124; |
|---|---:|---:|---:|
| `id_group` | **−1.2488** | **−21.83** | **0.000** |
| age | 0.0234 | 6.44 | 0.000 |
| income (high) | −0.3585 | −1.44 | 0.149 |
| sex | 0.0111 | 0.19 | 0.851 |

**Pseudo R² triples, and the added variable is the largest coefficient in the model by a factor of five.** How the question was framed predicts the answer better than income, education, age, wealth, debt and risk attitude put together.

This is a framing effect, and it is measurable here precisely because the survey randomised the starting point. Anchoring on 20% first and walking down produces systematically different discount rates than anchoring on 10% and walking up — from people drawn from the same population.

The practical reading: a discount rate elicited this way is partly a property of the instrument, not only of the person. Any downstream use of these numbers has to carry that caveat.

---

## 📈 What still holds after that

The respondent characteristics that survive as significant, from the model without `id_group`:

| Variable | Coefficient | Direction |
|---|---:|---|
| **Age** | +0.0222 *** | Older respondents are more impatient |
| **Farm or non-farm land** | +1.0742 *** | Property owners accept larger discounts |
| **All types of debt** | +1.1809 ** | Debt pushes hard toward immediate cash |
| **No debt** | −0.2981 ** | And its absence pulls the other way |
| **Self-employed** | −0.2984 *** | Less urgency for liquidity than employees |
| **Income** (all bands) | −0.20 to −0.41 *** | Higher income, more patience |
| **High precautionary savings** | −0.31 to −0.45 ** | A safety net buys the ability to wait |
| **Has a savings plan** | −0.2054 ** | Planning correlates with patience |
| **University degree** | +0.1961 ** | Counterintuitively, *more* impatient |
| Sex | +0.0094 | Not significant |
| Risk attitude | all n.s. | Not significant |

<sub>*** p<0.01, ** p<0.05</sub>

Two of these are worth pausing on.

**Sex drops out.** The descriptive tables show a clear gender gap — 25.65% of men choose 5% against 19.78% of women, while women more often accept 10% and 20%. In the regression that difference vanishes entirely (p = 0.851). It was never about gender; it was about the income, work status and education that differ by gender in this sample.

**Risk attitude is irrelevant.** None of the four risk categories reach significance. Impatience and risk aversion are commonly treated as the same underlying trait, and here they behave as separate things.

**Education runs the wrong way.** Degree holders are *more* willing to take the discount, not less. The report's reading is that they have somewhere to put the money — reinvestment, career acceleration — which makes immediate liquidity worth more to them, not less.

---

## 📂 Files

| File | Contents |
|---|---|
| [`exam2024.do`](exam2024.do) | The full Stata pipeline: variable construction, merges, survey design, nine nested ordered probits |
| [`Exam2024_Empirical.pdf`](Exam2024_Empirical.pdf) | Report with all frequency tables, boxplots and regression output |

The do-file merges seven SHIW 2010 archives on `nquest` (`q10c2`, `carcom10`, `q10a`, `lavoro`, `allb1`, `q10d`, `q10e`) and keeps household heads only (`nord==1`). The models are built up one block at a time, so each addition can be judged against the previous fit.

**Running it.** The SHIW microdata are distributed by the Bank of Italy and are not redistributed here. The path at the top of the do-file points at a university drive and needs changing to a local copy:

```stata
cd "X:\APPLIEDECONOMICS\EMPEC24\SHIW10\"
```

The observation count also drops from 13,996 in the descriptive tables to 2,003 in the regressions, because Stata deletes listwise across all sixteen covariates. The estimates are precise on that subsample, but it is a subsample.

---

## 📖 Source

Banca d'Italia, *Survey on Household Income and Wealth*, 2010 wave. Question C33 for the discount sequence, C34 for precautionary reserves, E13 for end-of-month financial condition, B1 question 7 for net payroll income.
