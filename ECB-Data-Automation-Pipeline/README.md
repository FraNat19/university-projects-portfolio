# 📊 ECB Data Automation Pipeline

**Automated workflow for processing European Central Bank interest rate data**


## 🎯 What This Project Does

Built an **automated data pipeline** using Microsoft Power Platform that:

1. ✅ **Detects** when a CSV file is uploaded to SharePoint
2. ✅ **Creates** a log entry with processing status
3. ✅ **Waits** 5 seconds (simulates data processing)
4. ✅ **Updates** the status to "Completed"
5. ✅ **Sends** an email notification with dashboard link
6. ✅ **Marks** the email as sent in the log

**Key Achievement:** Automated workflow with only manual dataset selection


## 🛠️ Technologies

- **Power Automate** - 6-step automated workflow
- **SharePoint Online** - File storage + tracking list
- **Power BI** - Dashboard visualization
- **Gmail API** - Email notifications


## 🔄 The Automation Flow

File Upload (CSV)
↓
Create Log Entry (Status: Pending)
↓
Delay 5 seconds
↓
Update Status (Completed)
↓
Send Email with Link
↓
Update Email Flag (Sent: Yes)


## 📈 Power BI Dashboard

Simple visualization of **ECB interest rates over time**:
- Date vs Interest Rate line chart
- Basic filtering by country
- Data: 274 monthly observations (2003-2024)


## 📊 Dataset

**Source:** ECB Statistical Data Warehouse  
**Type:** Monthly interest rates  
**Size:** 274 records  
**Period:** 2003-2024

