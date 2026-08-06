-- Validation queries for db62
-- One query per informational requirement listed in report_conclusivo_DB.pdf

USE db62;

-- Q1 - Demographics of every patient filed under a given clinical archive.
SELECT p.cf, p.nome, p.cognome, p.data_nascita, p.sesso, a.descrizione AS archivio
FROM Paziente p
JOIN Archivio a ON a.id = p.id_archivio
WHERE p.id_archivio = 1
ORDER BY p.cognome, p.nome;

-- Q2 - Record a new diagnosis of myocardial infarction for a patient.
INSERT INTO Diagnosi (paziente, esame, nome, descrizione, data_diagnosi)
VALUES ('RSSMRA80A01H501U', 3, 'Infarto miocardico acuto',
        'STEMI anteriore confermato da coronarografia', '2022-01-14');

-- Q3 - Every exam performed on each patient, with its outcome.
SELECT p.cognome, p.nome, e.tipo, e.data_esame, e.esito
FROM Paziente p
JOIN Esame e ON e.paziente = p.cf
ORDER BY p.cognome, e.data_esame;

-- Q4 - Mean patient age per recorded symptom.
SELECT s.descrizione AS sintomo,
       COUNT(DISTINCT s.paziente) AS pazienti,
       ROUND(AVG(TIMESTAMPDIFF(YEAR, p.data_nascita, s.data_comparsa)), 1) AS eta_media
FROM Sintomo s
JOIN Paziente p ON p.cf = s.paziente
GROUP BY s.descrizione
ORDER BY eta_media DESC;

-- Q5 - Youngest patient carrying a diagnosis of heart failure.
SELECT MIN(TIMESTAMPDIFF(YEAR, p.data_nascita, d.data_diagnosi)) AS eta_minima
FROM Diagnosi d
JOIN Paziente p ON p.cf = d.paziente
WHERE d.nome LIKE '%scompenso cardiaco%';

-- Q6 - Which drugs each doctor prescribed, and to which patients.
SELECT m.cognome AS medico, med.nome AS medicinale, med.principio_attivo,
       p.cognome AS paziente, pr.posologia, pr.data_prescr
FROM Prescrizione pr
JOIN Medico m       ON m.cf = pr.medico
JOIN Medicinale med ON med.codice_atc = pr.medicinale
JOIN Cura c         ON c.id = pr.cura
JOIN Diagnosi d     ON d.id = c.diagnosi
JOIN Paziente p     ON p.cf = d.paziente
ORDER BY m.cognome, p.cognome;

-- Q7 - Bed occupancy by room, for the ward logistics view.
SELECT st.codice AS stanza, st.tipo,
       COUNT(l.codice) AS letti,
       SUM(l.occupato) AS occupati,
       COUNT(l.codice) - SUM(l.occupato) AS liberi
FROM Stanza st
LEFT JOIN Letto l ON l.stanza = st.codice
GROUP BY st.codice, st.tipo
ORDER BY liberi DESC, st.codice;

-- Q8 - Complete clinical path of one patient, symptom to therapy.
SELECT p.cognome, s.descrizione AS sintomo, s.codice_triage,
       v.data_visita, e.tipo AS esame, d.nome AS diagnosi, c.nome AS cura
FROM Paziente p
LEFT JOIN Sintomo  s ON s.paziente = p.cf
LEFT JOIN Visita   v ON v.paziente = p.cf
LEFT JOIN Esame    e ON e.visita = v.id
LEFT JOIN Diagnosi d ON d.esame = e.id
LEFT JOIN Cura     c ON c.diagnosi = d.id
WHERE p.cf = 'RSSMRA80A01H501U'
ORDER BY v.data_visita;
