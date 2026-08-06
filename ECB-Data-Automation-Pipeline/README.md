<div align="center">

# 🏦 ECB Corporate Lending Pipeline

**A tracked, notified, hands-off path from CSV drop to published dashboard**

<img src="https://img.shields.io/badge/Power%20BI-F2C811?style=for-the-badge&logo=powerbi&logoColor=black" alt="Power BI">
<img src="https://img.shields.io/badge/Power%20Automate-0066FF?style=for-the-badge&logo=powerautomate&logoColor=white" alt="Power Automate">
<img src="https://img.shields.io/badge/SharePoint-036C70?style=for-the-badge&logo=microsoftsharepoint&logoColor=white" alt="SharePoint">
<img src="https://img.shields.io/badge/ECB-SDW-8b1a1a?style=for-the-badge" alt="ECB">

</div>

---

## 📋 What it does

Drop a CSV into SharePoint. Everything after that is automatic: the file is logged, processed, marked complete, and a notification goes out with a link to the refreshed dashboard.

```
File uploaded (CSV)
        ↓
Create log entry — Status: Pending
        ↓
Delay 5s (processing placeholder)
        ↓
Update Status: Completed
        ↓
Send email with dashboard link
        ↓
Set Email_Sent: Yes
```

Six steps, one manual action: choosing which file to upload.

The interesting choice is the **log-first ordering**. The row is written with `Pending` *before* any work happens, then updated to `Completed` afterwards. If the flow dies halfway, the evidence survives as a stuck row. A flow that only logs on success cannot tell you about the runs that failed.

The five-second delay stands in for a real transformation step. In this build it does nothing, which is worth saying rather than dressing up: it marks the place where parsing, validation or an incremental load would go.

---

## 📊 The data

[`ECB Data Portal_20251129191514.csv`](ECB%20Data%20Portal_20251129191514.csv) — **273 monthly observations, January 2003 to 2024.**

| | |
|---|---|
| **Source** | ECB Statistical Data Warehouse |
| **Series** | `MIR.M.U2.B.A2I.AM.R.A.2240.EUR.N` |
| **Meaning** | Cost of borrowing for corporations, euro area, monthly |
| **First value** | 4.54% in January 2003 |

Reading the series key is the fastest way to know what you are looking at: `M` monthly, `U2` euro area aggregate, `EUR` denominated, `MIR` the MFI interest rate domain. This is a **single euro-area aggregate**, not a country breakdown, so the dashboard shows one line through time rather than a comparison across member states.

Twenty-two years of monthly observations covers the 2008 crisis, the sovereign debt crisis, the negative-rate period and the 2022 tightening. Corporate borrowing cost is the transmission channel: it shows how far and how fast policy rates actually reach firms.

---

## 📈 The dashboard

[`ECB_Corporate_Lending_Dashboard.pbix`](ECB_Corporate_Lending_Dashboard.pbix) — open in Power BI Desktop.

A time series of borrowing cost against date, with filtering on the period. Deliberately plain: the point of this project is the pipeline behind it, not the chart on top.

---

## ⚙️ The flow

[`ECBDataAuto-ProcessingPipeline_20251130140453.zip`](ECBDataAuto-ProcessingPipeline_20251130140453.zip) is the exported Power Automate solution.

| Step | Action | Why |
|---|---|---|
| 1 | Trigger on file creation in SharePoint | No polling, no scheduled job |
| 2 | Create log row, `Status: Pending` | Evidence exists before work starts |
| 3 | Delay 5 seconds | Placeholder for real processing |
| 4 | Update `Status: Completed` | Closes the loop |
| 5 | Send email with dashboard link | Notification carries the destination |
| 6 | Set `Email_Sent: Yes` | Prevents duplicate notification on re-run |

**To reuse it** you need a SharePoint library and a tracking list with `Status` and `Email_Sent` columns, a mail connection, and your own dashboard URL in step 5.

A companion project, [`Country-Economic-Dashboard-Automation`](../Country-Economic-Dashboard-Automation), takes the same pattern further: it parses the filename and emails a link filtered to that country.

---

## 📖 Source

European Central Bank, [Statistical Data Warehouse](https://data.ecb.europa.eu). Series `MIR.M.U2.B.A2I.AM.R.A.2240.EUR.N`, cost of borrowing for corporations.
