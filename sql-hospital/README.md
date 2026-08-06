<div align="center">

# 🫀 Cardiology Ward Database

**From a bedside interview to a normalised MySQL schema**

<img src="https://img.shields.io/badge/MySQL-4479A1?style=for-the-badge&logo=mysql&logoColor=white" alt="MySQL">
<img src="https://img.shields.io/badge/Python-3776AB?style=for-the-badge&logo=python&logoColor=white" alt="Python">
<img src="https://img.shields.io/badge/3NF-normalised-1a7f37?style=for-the-badge" alt="3NF">
<img src="https://img.shields.io/badge/Sapienza-8b1a1a?style=for-the-badge" alt="Sapienza">

<sub>Database Systems · BSc · Sapienza University of Rome · January 2022</sub>

</div>

---

## 📋 Overview

The requirements for this database did not come from a specification sheet. They came from a structured interview with **Claudio, a cardiology resident at Policlinico Universitario Agostino Gemelli**, who was asked to describe what actually happens to a patient who walks into the ward with chest pain.

Two connected worlds came out of that conversation, and they shaped everything downstream:

- **The clinical path** — symptom, triage, admission, history, exams, diagnosis, therapy, follow-up.
- **The documentary path** — how each of those events becomes a report, and how reports are filed so a patient's history can be reconstructed years later.

```
Sintomo → Triage → Visita → Esame → Diagnosi → Cura → Referto → Archivio
```

| | |
|---|---|
| **Domain** | Cardiology ward: 25 rooms, 2 beds each, ~60 medical staff |
| **Method** | Interview → E-R model → logical schema → MySQL implementation → test queries → Python interface |
| **Design strategy** | Inside-out, starting from the concepts central to the interview |
| **Normalisation** | Third normal form on all principal entities |
| **Database** | `db62`, built in MySQL Workbench |

---

## 📂 What is in this folder

| File | Contents |
|---|---|
| [`schema.sql`](schema.sql) | 21 tables, primary and foreign keys, check constraints, indexes |
| [`seed_data.sql`](seed_data.sql) | Sample records reproducing typical cardiology pathways |
| [`queries.sql`](queries.sql) | The validation queries, one per informational requirement |
| [`validate.py`](validate.py) | Runs all three without needing a MySQL server |
| [`report_conclusivo_DB.pdf`](report_conclusivo_DB.pdf) | Full report: interview, E-R diagram, logical schema, physical model |

```
$ python validate.py

schema.sql: 37 statements executed
seed_data.sql: 21 statements executed
tables created: 21
foreign key violations: 0
all statements executed, referential integrity holds
```

`validate.py` translates the MySQL dialect to SQLite in memory, loads the schema and sample data with foreign keys enforced, then runs every query and prints its output. A broken reference or a typo in a join fails the run.

---

## 📊 Sizing the database

Volumes were estimated from the interview, annually, and used to justify the choice of keys and indexes.

| Entity | Type | Annual volume | | Entity | Type | Annual volume |
|---|:-:|---:|---|---|:-:|---:|
| Referto | E | 30,000 | | Stanza | E | 25 |
| Diagnosi | E | 30,000 | | Infermiere | E | 22 |
| Esame | E | 18,000 | | OSS | E | 12 |
| Sintomo | E | 15,000 | | Letto | E | 50 |
| Paziente | E | 12,000 | | Turno | E | 4 |
| Visita | E | 12,000 | | Composizione | R | 3 |
| Cura | E | 11,000 | | Archivio | E | 1 |
| Medico | E | 60 <sub>(20 strutturati)</sub> | | Rappresentante | E | 1 |

Two relations carry volume of their own: `Svolgimento_Visita` and `Esaminazione_Medico`, both at 12,000.

The shape of that table drove the design. **Referti and diagnosi are the heaviest tables, at more than twice the patient count**, because one admission generates many documents. Meanwhile rooms, beds and shifts are tiny and effectively static. Indexing effort went where the rows are: `Referto(paziente, data_referto)`, `Diagnosi(paziente)`, `Esame(paziente, data_esame)`.

### Operations

| # | Operation | Type | Annual frequency |
|---|---|:-:|---:|
| 1 | Register patient demographics at first access or admission | I | 12,000 |
| 2 | File a new report in the clinical archive | B | 30,000 |
| 3 | Carry out a visit, updating diagnosis and therapy | I | 12,000 |
| 4 | Prescribe a drug or a specific therapy | I | 18,000 |
| 5 | Issue a materials order through the representative | I | 200 |
| 6 | Report on an exam and attach it to the patient | B | 30,000 |
| 7 | Record a new symptom reported by the patient | I | 15,000 |
| 8 | Perform a diagnostic exam and record its outcome | I | 18,000 |

<sub>I = interactive, B = batch.</sub>

The two batch operations are also the two highest-volume ones, which is convenient: the heaviest work is the work that does not need to answer in real time.

---

## 🗂️ Schema

Twenty-one tables in four groups.

**Clinical path** — `Paziente`, `Sintomo`, `Visita`, `Esame`, `Diagnosi`, `Cura`

**Documents** — `Referto`, `Archivio`

**Staff and structure** — `Medico`, `Specializzando`, `Infermiere`, `OSS`, `Turno`, `Stanza`, `Letto`

**Pharmacy** — `Medicinale`, `Prescrizione`, `Rappresentante`, `Ordine`

**Associative** — `Composizione`, `Esaminazione_Medico`

```sql
CREATE TABLE Referto (
    id            INT AUTO_INCREMENT PRIMARY KEY,
    paziente      CHAR(16)     NOT NULL,
    id_archivio   INT          NOT NULL,
    visita        INT          NULL,
    esame         INT          NULL,
    data_referto  DATETIME     NOT NULL,
    contenuto     TEXT         NOT NULL,
    CONSTRAINT fk_ref_paziente FOREIGN KEY (paziente) REFERENCES Paziente(cf)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT chk_ref_origine CHECK (visita IS NOT NULL OR esame IS NOT NULL)
) ENGINE=InnoDB;
```

Three modelling decisions worth pointing at:

**A report must come from somewhere.** `Referto` links optionally to a visit *and* optionally to an exam, but `chk_ref_origine` requires at least one of them. A document with no clinical event behind it cannot be filed.

**Clinical history is never deleted.** Foreign keys from `Visita`, `Esame`, `Diagnosi` and `Referto` back to `Paziente` use `ON DELETE RESTRICT`. A patient with any recorded history cannot be removed. `Sintomo` is the exception and cascades, since a symptom has no meaning detached from the patient who reported it.

**Prescriptions are traceable to a doctor.** `Prescrizione` carries a foreign key to `Medico` alongside the one to `Cura`, so the question "who prescribed this" is answerable directly rather than inferred from the visit.

---

## 🔍 Validation queries

Eight queries, each covering an informational requirement from the analysis. All of them run in `validate.py`.

| | Question |
|---|---|
| Q1 | Demographics of every patient filed under a given archive |
| Q2 | Record a new myocardial infarction diagnosis |
| Q3 | Every exam performed on each patient, with outcomes |
| Q4 | Mean patient age per recorded symptom |
| Q5 | Youngest patient carrying a heart failure diagnosis |
| Q6 | Which drugs each doctor prescribed, and to which patients |
| Q7 | Bed occupancy by room |
| Q8 | One patient's complete path, symptom to therapy |

Q6 is the one that exercises the schema hardest, crossing five tables to connect a prescribing doctor to the patient who ends up taking the drug:

```sql
SELECT m.cognome AS medico, med.nome AS medicinale, p.cognome AS paziente, pr.posologia
FROM Prescrizione pr
JOIN Medico m       ON m.cf = pr.medico
JOIN Medicinale med ON med.codice_atc = pr.medicinale
JOIN Cura c         ON c.id = pr.cura
JOIN Diagnosi d     ON d.id = c.diagnosi
JOIN Paziente p     ON p.cf = d.paziente
ORDER BY m.cognome, p.cognome;
```

```
medico   | medicinale  | paziente | posologia
Marroni  | Ramipril    | Bianco   | 5 mg once daily
Marroni  | Bisoprololo | Bianco   | 2.5 mg once daily
Marroni  | Furosemide  | Marino   | 25 mg twice daily
Verdi    | Rivaroxaban | Greco    | 20 mg once daily
```

Q7 answers a logistics question rather than a clinical one, which is the point of having modelled rooms and beds at all:

```
stanza | tipo                | letti | occupati | liberi
ST-103 | osservazione        | 2     | 0        | 2
ST-201 | preospedalizzazione | 2     | 0        | 2
ST-101 | degenza             | 2     | 1        | 1
ST-102 | degenza             | 2     | 2        | 0
```

---

## ⚙️ Running it

Against a real MySQL server:

```bash
mysql -u root -p < schema.sql
mysql -u root -p < seed_data.sql
mysql -u root -p db62 < queries.sql
```

Without one, using only the Python standard library:

```bash
python validate.py
```

The original project also included a Python interface built on `mysql.connector` and tested in Colab, exposing the queries above through a small menu so the database could be explored without writing SQL by hand. That notebook is described in the report but is not part of this folder.

---

## 👥 Team

| Name | Role |
|--------|------|
| **Francesco Natali** | Design & Implementation |
| **Davide Di Brango** | Design & Implementation |
| **Davide Anello** | Design & Implementation |
| **Leonardo Agate** | Design & Implementation |
