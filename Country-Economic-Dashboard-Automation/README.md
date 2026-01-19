# Country Economic Dashboard Automation

**Fully automated pipeline for country-specific economic analysis with Power BI, Power Automate, and SharePoint**

##  What This Project Does

Built an **end-to-end automated system** that:

1. Detects when a country CSV file is uploaded to SharePoint (e.g., `Italy_Economy.csv`)
2. Extracts the country name from the filename
3. Creates a tracking record in SharePoint list
4. Sends email notification with personalized dashboard link
5. Dashboard opens **automatically filtered** for that country

**Result:** Upload file → 10 seconds → Receive email with ready-to-view dashboard

## Architecture

Upload "Italy_Economy.csv"
↓
Power Automate extracts "Italy" from filename
↓
Creates log: Country=Italy, Records=25, Status=Completed
↓
Sends email: "Italy Dashboard Ready" + filtered link
↓
Dashboard shows ONLY Italy data (2000-2024)


## 🛠️ Tech Stack

- **Power BI** - Interactive dashboard (GDP, unemployment, inflation trends)
- **Power Automate** - 5-step automated workflow
- **SharePoint** - File storage + tracking list
- **World Bank Data** - 6 countries, 25 years, 5 economic indicators


## Dashboard Features

**Data:** 150 records (6 countries × 25 years)

**Visualizations:**
- 📈 GDP Trend (2000-2024) - Area chart
- 📉 Unemployment Rate Evolution
- 💹 GDP Growth vs Inflation (dual-axis)
- 📋 Detailed data table
- 🌍 Country filter (Italy, Germany, France, Spain, UK, USA)

## Automation Workflow

### Trigger:
File uploaded to SharePoint Documents

### Steps:
1. **Extract country name** from filename (`Italy_Economy.csv` → `Italy`)
2. **Create SharePoint record** with country, date, status
3. **Send email** with dashboard link filtered by country: dashboard_url?filter=All_Countries_Economy/Country eq 'Italy'
4. **Update status** Email_Sent = Yes

## Data Source

**World Bank - World Development Indicators**
- 6 countries: Italy, Germany, France, Spain, UK, USA
- 25 years: 2000-2024
- 5 indicators: GDP, GDP growth, inflation, unemployment, population
- [Download data](https://data.worldbank.org)
