<div align="center">

# 🌍 Country Economic Dashboard Automation

**Upload a CSV, get a personalised dashboard in your inbox**

<img src="https://img.shields.io/badge/Power%20BI-F2C811?style=for-the-badge&logo=powerbi&logoColor=black" alt="Power BI">
<img src="https://img.shields.io/badge/Power%20Automate-0066FF?style=for-the-badge&logo=powerautomate&logoColor=white" alt="Power Automate">
<img src="https://img.shields.io/badge/SharePoint-036C70?style=for-the-badge&logo=microsoftsharepoint&logoColor=white" alt="SharePoint">
<img src="https://img.shields.io/badge/World%20Bank-WDI-1a7f37?style=for-the-badge" alt="World Bank">

</div>

---

## 📋 What it does

Drop `Italy_Economy.csv` into a SharePoint folder. Ten seconds later an email arrives with a link that opens a Power BI dashboard already filtered to Italy.

Nobody clicks through a report and applies a slicer. The filter travels in the URL, generated from the filename.

```
Upload "Italy_Economy.csv"
        ↓
Power Automate extracts "Italy" from the filename
        ↓
Creates a SharePoint log row: Country=Italy, Records=25, Status=Completed
        ↓
Sends "Italy Dashboard Ready" with a pre-filtered link
        ↓
Sets Email_Sent = Yes
```

The detail that makes it feel finished is the filter in the query string:

```
dashboard_url?filter=All_Countries_Economy/Country eq 'Italy'
```

Power BI applies it on load, so the recipient never sees the other five countries. One dashboard serves every country without being duplicated.

---

## 📊 The data

[`All_Countries_Economy.csv`](All_Countries_Economy.csv) — **150 rows: 6 countries × 25 years**, from the World Bank World Development Indicators.

| | |
|---|---|
| **Countries** | Italy, Germany, France, Spain, United Kingdom, United States |
| **Period** | 2000–2024 |
| **Format** | Semicolon-separated, comma decimal separator |

Eight measures per country-year, five pulled directly and three derived:

| Column | Source |
|---|---|
| `GDP_Current_USD` | WDI |
| `GDP_Growth_Pct` | WDI |
| `Inflation_CPI_Pct` | WDI |
| `Population_15_64` | WDI |
| `Unemployment_Pct` | WDI |
| `GDP_per_capita_15_64` | derived |
| `GDP_per_employed_15_64` | derived |
| `Real_Growth_minus_Inflation` | derived |

The derived columns are the ones worth having. Dividing GDP by working-age population rather than by total population strips out demographic structure, so Italy and Germany become comparable on productivity instead of on age profile. And carrying nominal growth minus inflation as a column keeps the dashboard from reporting a boom that is only price level.

---

## 📈 The dashboard

[`Country_Economic_Dashboard.pbix`](Country_Economic_Dashboard.pbix) — open in Power BI Desktop.

- GDP trend 2000–2024, area chart
- Unemployment rate evolution
- GDP growth against inflation on a dual axis
- Detailed data table
- Country slicer, which the emailed link presets

The dual-axis view is where 2008 and 2020 separate: countries whose growth and inflation move together behave differently from those where the two diverge.

---

## ⚙️ The flow

[`Country_Economic_Report_Generator_20251130214959.zip`](Country_Economic_Report_Generator_20251130214959.zip) is the exported Power Automate solution, importable into another tenant.

| Step | Action |
|---|---|
| 1 | Trigger: file created in SharePoint Documents |
| 2 | Parse the country name from the filename, before `_Economy.csv` |
| 3 | Create a tracking row: country, timestamp, record count, status |
| 4 | Send email containing the filtered dashboard URL |
| 5 | Set `Email_Sent = Yes` on the tracking row |

The tracking list is what makes this operable rather than a demo. Every upload leaves a row, so a failed run shows up as a record stuck at `Pending` instead of disappearing silently.

**To reuse it** you need a SharePoint document library and a tracking list, a published Power BI report whose dataset exposes a `Country` column, and your own report URL in step 4.

---

## 📖 Source

World Bank, [World Development Indicators](https://data.worldbank.org). Series used: GDP (current US$), GDP growth (annual %), inflation consumer prices (annual %), population ages 15–64, unemployment (% of total labour force).
