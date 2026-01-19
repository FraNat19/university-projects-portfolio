# Database Project: Cardiology Ward Management System

<p align="center">
  <img src="https://img.shields.io/badge/Tools-MySQL%20Workbench%20%7C%20Python-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Focus-ER%20Model%20%7C%20SQL%20%7C%20Implementation-blue?style=for-the-badge" alt="Focus">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Year-2021%2F2022-orange?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>From requirements elicitation to MySQL implementation: Cardiology Ward Database</b>
</p>

---

## Overview

This repository contains a group project for the **Database Systems** course during my Bachelor's degree. The project involved designing and implementing a realistic database covering the full lifecycle: from requirements gathering to physical implementation and testing.

**Modeled domain:** clinical and administrative workflow of a **Cardiology Ward**.

```
Paziente → Sintomo → Visita/Esami → Diagnosi → Cura → Referti/Archivio
```

---

## Project Goals

| Phase | Description |
|-------|-------------|
| **Requirements Gathering** | Structured interviews with a cardiology resident (Policlinico Gemelli) |
| **Design** | Conceptual ER model → Logical relational schema → Physical MySQL implementation |
| **Validation** | Test data, representative SQL queries, Python interface for interaction |

---

## Design & Implementation

### Domain Analysis

The modeled clinical workflow covers the full patient path:

```
┌─────────┐    ┌─────────┐    ┌─────────┐    ┌─────────┐    ┌─────────┐
│ TRIAGE  │ → │ PRESA   │ → │ANAMNESI │ → │  ESAMI  │ → │DIAGNOSI │
│         │    │IN CARICO│    │         │    │ (ECG,..)│    │         │
└─────────┘    └─────────┘    └─────────┘    └─────────┘    └─────────┘
                                                                  ↓
┌─────────┐    ┌─────────┐    ┌─────────┐                  ┌─────────┐
│ARCHIVIO │ ← │ REFERTO │ ← │FOLLOW-UP│ ←──────────────── │ TERAPIA │
│         │    │         │    │         │                  │         │
└─────────┘    └─────────┘    └─────────┘                  └─────────┘
```


### ER Model (Conceptual)

**Main Entities:**

| Entity | Description | Key Attributes |
|--------|-------------|----------------|
| `PATIENT` | Patient demographic info | ID, FirstName, LastName, BirthDate |
| `DOCTOR` | Medical staff | StaffID, Specialty |
| `VISIT` | Clinical encounter | ID, Date, Time, Type |
| `EXAM` | Diagnostic test | ID, Type (ECG, Echo, etc.), Result |
| `DIAGNOSIS` | Identified condition | ID, ICD Code, Description |
| `TREATMENT` | Therapy plan | ID, Type (Medication/Procedure) |
| `MEDICATION` | Prescribed drugs | ATC Code, Name, Dosage |
| `REPORT` | Clinical document | ID, Date, Content |

**Relationships:**

```
PAZIENTE ──(1,N)── VISITA ──(N,1)── MEDICO
    │                 │
    │                 ├──(1,N)── ESAME
    │                 │
    └──(1,N)── DIAGNOSI ──(1,N)── TERAPIA ──(N,M)── MEDICINALE
```


---

### Estimated Volumes & Operations

**Annual Volumes:**

| Entity/Relation | Estimated Volume |
|-----------------|-----------------|
| Patient | 5,000 |
| Visit | 15,000 |
| Exam | 30,000 |
| Diagnosis | 8,000 |
| Report | 20,000 |
| Doctor | 50 |

**Operations Table:**

| Operation | Type | Frequency |
|-----------|------|-----------|
| Add new visit | Interactive | 50/day |
| Retrieve patient history | Interactive | 100/day |
| Monthly diagnosis report | Batch | 1/month |
| Report archive backup | Batch | 1/week |

---

### Physical Implementation

**Database:** MySQL  

```sql
-- Example: create PATIENT table
CREATE TABLE Patient (
    ID VARCHAR(16) PRIMARY KEY,
    FirstName VARCHAR(50) NOT NULL,
    LastName VARCHAR(50) NOT NULL,
    BirthDate DATE NOT NULL,
    Gender ENUM('M', 'F') NOT NULL,
    Phone VARCHAR(15),
    Email VARCHAR(100),
    Address VARCHAR(200)
);

-- Example Foreign Key
ALTER TABLE Visit
ADD CONSTRAINT fk_visit_patient
FOREIGN KEY (PatientID) REFERENCES Patient(ID)
ON DELETE RESTRICT ON UPDATE CASCADE;

**Indexes:**
# 🗄️ Database Project (Real-World DB)

<p align="center">
  <img src="https://img.shields.io/badge/Tools-MySQL%20Workbench%20%7C%20Python-green?style=for-the-badge" alt="Tools">
  <img src="https://img.shields.io/badge/Focus-ER%20Model%20%7C%20SQL%20%7C%20Implementation-blue?style=for-the-badge" alt="Focus">
  <img src="https://img.shields.io/badge/University-Sapienza-darkred?style=for-the-badge" alt="Sapienza">
  <img src="https://img.shields.io/badge/Year-2021%2F2022-orange?style=for-the-badge" alt="Year">
</p>

<p align="center">
  <b>From requirements elicitation to MySQL implementation: Cardiology ward database</b>
</p>

---

## 📋 Overview

Repository for the **Database** exam project, developed as a group assignment: building a "real-world" database following the complete lifecycle, from requirements gathering to implementation and testing.

**Domain modeled**: organization and clinical-documentary pathway of a **Cardiology ward**.

```
Patient → Symptom → Visit/Exams → Diagnosis → Treatment → Reports/Archive
```

---

## 👥 Team

| Member | Role |
|--------|------|
| **Francesco Natali** | Design & Implementation |
| **Davide Di Brango** | Design & Implementation |
| **Davide Anello** | Design & Implementation |
| **Leonardo Agate** | Design & Implementation |

---

## 🎯 Project Objectives

| Phase | Description |
|-------|-------------|
| 📝 **Requirements Gathering** | Structured interview with a cardiovascular disease resident (Policlinico Gemelli) |
| 📐 **Design** | Conceptual ER schema → Logical relational schema → Physical MySQL implementation |
| ✅ **Validation** | Test data, representative SQL queries, Python interface |

---

## 📁 Repository Structure

```
cardiology-database/
│
├── 📂 docs/
│   ├── Project_Report.pdf          # Complete documentation
│   ├── Requirements_Interview.pdf  # Interview transcript
│   └── Volume_Operations_Table.xlsx
│
├── 📂 diagrams/
│   ├── ER_Conceptual_Schema.png    # ER diagram
│   ├── ER_Conceptual_Schema.mwb    # MySQL Workbench file
│   └── Logical_Schema.png          # Relational schema
│
├── 📂 sql/
│   ├── 01_create_database.sql      # DB and tables creation
│   ├── 02_constraints.sql          # PK, FK, constraints
│   ├── 03_indexes.sql              # Optimization indexes
│   ├── 04_insert_test_data.sql     # ~80 test records
│   └── 05_queries.sql              # Validation queries
│
├── 📂 python/
│   ├── db_interface.ipynb          # Colab notebook with interface
│   └── db_connector.py             # MySQL connection module
│
├── 📂 data/
│   └── db62_dump.sql               # Complete database dump
│
└── README.md
```

---

## 🧱 Design & Implementation

### Domain Analysis

The modeled clinical pathway covers the entire flow of a cardiology patient:

```
┌─────────┐    ┌─────────┐    ┌─────────┐    ┌─────────┐    ┌─────────┐
│ TRIAGE  │ → │ PATIENT │ → │ MEDICAL │ → │  EXAMS  │ → │DIAGNOSIS│
│         │    │ INTAKE  │    │ HISTORY │    │ (ECG,..)│    │         │
└─────────┘    └─────────┘    └─────────┘    └─────────┘    └─────────┘
                                                                  ↓
┌─────────┐    ┌─────────┐    ┌─────────┐                  ┌─────────┐
│ ARCHIVE │ ← │ REPORT  │ ← │FOLLOW-UP│ ←──────────────── │TREATMENT│
│         │    │         │    │         │                  │         │
└─────────┘    └─────────┘    └─────────┘                  └─────────┘
```

### ER Schema (Conceptual)

**Main Entities:**

| Entity | Description | Key Attributes |
|--------|-------------|----------------|
| `PATIENT` | Patient registry | SSN, FirstName, LastName, BirthDate |
| `PHYSICIAN` | Healthcare staff | ID, Specialization |
| `VISIT` | Clinical encounter | ID, Date, Time, Type |
| `EXAM` | Diagnostic test | ID, Type (ECG, Echo, ...), Result |
| `DIAGNOSIS` | Identified pathology | ID, ICD Code, Description |
| `TREATMENT` | Therapeutic plan | ID, Type (drug/procedure) |
| `MEDICATION` | Prescribed drug | ATC Code, Name, Dosage |
| `REPORT` | Clinical document | ID, Date, Content |

**Relationships:**

```
PATIENT ──(1,N)── VISIT ──(N,1)── PHYSICIAN
    │                 │
    │                 ├──(1,N)── EXAM
    │                 │
    └──(1,N)── DIAGNOSIS ──(1,N)── TREATMENT ──(N,M)── MEDICATION
```

### Volume & Operations Estimation

**Volume Table (annual):**

| Entity/Relationship | Estimated Volume |
|---------------------|------------------|
| Patient | 5,000 |
| Visit | 15,000 |
| Exam | 30,000 |
| Diagnosis | 8,000 |
| Report | 20,000 |
| Physician | 50 |

**Operations Table:**

| Operation | Type | Frequency |
|-----------|------|-----------|
| Insert new visit | Interactive | 50/day |
| Search patient history | Interactive | 100/day |
| Monthly diagnosis report | Batch | 1/month |
| Archive reports backup | Batch | 1/week |

### Physical Implementation

**Database**: `db62` (MySQL 8.0)

```sql
-- Example: PATIENT table creation
CREATE TABLE Patient (
    SSN VARCHAR(16) PRIMARY KEY,
    FirstName VARCHAR(50) NOT NULL,
    LastName VARCHAR(50) NOT NULL,
    BirthDate DATE NOT NULL,
    Gender ENUM('M', 'F') NOT NULL,
    Phone VARCHAR(15),
    Email VARCHAR(100),
    Address VARCHAR(200)
);

-- Foreign Key example
ALTER TABLE Visit
ADD CONSTRAINT fk_visit_patient
FOREIGN KEY (PatientSSN) REFERENCES Patient(SSN)
ON DELETE RESTRICT ON UPDATE CASCADE;
```

**Indexes created:**

```sql
-- Indexes for frequent searches
CREATE INDEX idx_patient_lastname ON Patient(LastName);
CREATE INDEX idx_visit_date ON Visit(Date);
CREATE INDEX idx_exam_type ON Exam(Type);
```

**Test data**: ~80 records inserted to simulate realistic clinical pathways.

---

## Queries & Validation

### Test SQL Queries

**1. Patient-Exams Join (clinical history):**
```sql
SELECT p.SSN, p.LastName, p.FirstName, e.Type, e.Date, e.Result
FROM Patient p
JOIN Visit v ON p.SSN = v.PatientSSN
JOIN Exam e ON v.ID = e.VisitID
WHERE p.SSN = 'RSSMRA80A01H501Z'
ORDER BY e.Date DESC;
```

**2. Average age aggregation by symptom:**
```sql
SELECT s.Description AS Symptom,
       AVG(YEAR(CURDATE()) - YEAR(p.BirthDate)) AS AvgAge,
       COUNT(*) AS PatientCount
FROM Patient p
JOIN Visit v ON p.SSN = v.PatientSSN
JOIN Symptom s ON v.ID = s.VisitID
GROUP BY s.Description
ORDER BY PatientCount DESC;
```

**3. Prescriptions Physician → Patient → Medication:**
```sql
SELECT ph.LastName AS Physician,
       p.LastName AS Patient,
       med.Name AS Drug,
       t.Dosage,
       t.Duration
FROM Physician ph
JOIN Visit v ON ph.ID = v.PhysicianID
JOIN Patient p ON v.PatientSSN = p.SSN
JOIN Treatment t ON v.ID = t.VisitID
JOIN Medication med ON t.MedicationID = med.Code
WHERE ph.ID = 'MED001';
```

**4. Diagnosis count by pathology:**
```sql
SELECT d.ICDCode, d.Description, COUNT(*) AS Occurrences
FROM Diagnosis d
GROUP BY d.ICDCode, d.Description
ORDER BY Occurrences DESC
LIMIT 10;
```

---

## Python Interface

Colab notebook with MySQL connection and functions for frequent queries.

### Database Connection

```python
import mysql.connector
import pandas as pd

def connect_db():
    return mysql.connector.connect(
        host="localhost",
        user="root",
        password="password",
        database="db62"
    )

def execute_query(query):
    conn = connect_db()
    df = pd.read_sql(query, conn)
    conn.close()
    return df
```

### Interactive Menu

```python
def menu():
    print("=" * 50)
    print("  CARDIOLOGY INFORMATION SYSTEM - db62")
    print("=" * 50)
    print("1. Search patient by SSN")
    print("2. Patient exam history")
    print("3. Visit list by date")
    print("4. Monthly diagnosis report")
    print("5. Prescriptions by physician")
    print("0. Exit")
    print("=" * 50)
    return input("Select option: ")
```

### Example Output

```
┌────────────┬──────────┬────────┬────────────┬─────────┐
│ SSN        │ LastName │ First  │ Exam Type  │ Result  │
├────────────┼──────────┼────────┼────────────┼─────────┤
│ RSSMRA80.. │ Rossi    │ Mario  │ ECG        │ Normal  │
│ RSSMRA80.. │ Rossi    │ Mario  │ Echocardio │ Mild..  │
│ RSSMRA80.. │ Rossi    │ Mario  │ Holter     │ Normal  │
└────────────┴──────────┴────────┴────────────┴─────────┘
```

---

## Tools & Stack

| Tool | Purpose |
|------|---------|
| **MySQL 8.0** | Relational DBMS |
| **MySQL Workbench** | ER design & administration |
| **Python 3.x** | Query interface |
| **pandas** | Tabular result display |
| **mysql-connector-python** | Connection driver |
| **Google Colab** | Notebook environment |




