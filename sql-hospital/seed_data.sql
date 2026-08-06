-- Sample records for db62, reproducing typical cardiology pathways
USE db62;

INSERT INTO Stanza (codice, piano, tipo) VALUES
('ST-101', 1, 'degenza'), ('ST-102', 1, 'degenza'),
('ST-103', 1, 'osservazione'), ('ST-201', 2, 'preospedalizzazione');

INSERT INTO Letto (codice, stanza, occupato) VALUES
('LT-101-A','ST-101',TRUE), ('LT-101-B','ST-101',FALSE),
('LT-102-A','ST-102',TRUE), ('LT-102-B','ST-102',TRUE),
('LT-103-A','ST-103',FALSE), ('LT-103-B','ST-103',FALSE),
('LT-201-A','ST-201',FALSE), ('LT-201-B','ST-201',FALSE);

INSERT INTO Turno (id, nome, ora_inizio, ora_fine) VALUES
(1,'mattina','07:00:00','14:00:00'),
(2,'pomeriggio','14:00:00','21:00:00'),
(3,'notte','21:00:00','07:00:00'),
(4,'riposo',NULL,NULL);

INSERT INTO Medico (cf, nome, cognome, specializzazione, strutturato) VALUES
('BNCLCU70M12H501A','Luca','Bianchi','Cardiologia interventistica',TRUE),
('VRDGNN75F55H501B','Giovanna','Verdi','Elettrofisiologia',TRUE),
('MRRPLA68B20H501C','Paolo','Marroni','Cardiologia clinica',TRUE),
('CNTSRA82D44H501D','Sara','Conti','Ecocardiografia',FALSE);

INSERT INTO Specializzando (cf, nome, cognome, anno_corso, tutor) VALUES
('FRRCLD95T10H501E','Claudio','Ferrari',4,'BNCLCU70M12H501A'),
('GLLMRA96A41H501F','Marta','Galli',2,'VRDGNN75F55H501B');

INSERT INTO Infermiere (cf, nome, cognome, turno) VALUES
('RSITRS85C45H501G','Teresa','Rosi',1),
('LMBGRG88E03H501H','Giorgio','Lombardi',2),
('PLLNNA90H62H501I','Anna','Pellegrini',3);

INSERT INTO OSS (cf, nome, cognome, turno) VALUES
('DMCLGU92L15H501J','Luigi','De Micheli',1),
('SNTRTA87P51H501K','Rita','Santoro',2);

INSERT INTO Archivio (id, descrizione, data_apertura) VALUES
(1,'Archivio clinico Cardiologia - degenza','2021-09-01'),
(2,'Archivio clinico Cardiologia - ambulatoriale','2021-09-01');

INSERT INTO Paziente (cf, nome, cognome, data_nascita, sesso, telefono, id_archivio, letto) VALUES
('RSSMRA80A01H501U','Mario','Rossi','1958-01-01','M','3331112221',1,'LT-101-A'),
('BNCLRA65E41H501V','Laura','Bianco','1965-05-01','F','3332223332',1,'LT-102-A'),
('GRECRL72T22H501W','Carlo','Greco','1972-12-22','M','3333334443',1,'LT-102-B'),
('FNTNNA90M50H501X','Anna','Fontana','1990-08-10','F','3334445554',2,NULL),
('MRNGPP48S30H501Y','Giuseppe','Marino','1948-11-30','M','3335556665',2,NULL);

INSERT INTO Sintomo (paziente, descrizione, data_comparsa, codice_triage) VALUES
('RSSMRA80A01H501U','Dolore toracico','2022-01-12','rosso'),
('BNCLRA65E41H501V','Dispnea','2022-01-13','giallo'),
('GRECRL72T22H501W','Palpitazioni','2022-01-15','verde'),
('FNTNNA90M50H501X','Palpitazioni','2022-01-18','verde'),
('MRNGPP48S30H501Y','Dispnea','2022-01-20','arancione');

INSERT INTO Visita (paziente, medico, data_visita, tipo) VALUES
('RSSMRA80A01H501U','BNCLCU70M12H501A','2022-01-12 08:30:00','urgenza'),
('BNCLRA65E41H501V','MRRPLA68B20H501C','2022-01-13 10:00:00','prima_visita'),
('GRECRL72T22H501W','VRDGNN75F55H501B','2022-01-15 11:15:00','prima_visita'),
('FNTNNA90M50H501X','CNTSRA82D44H501D','2022-01-18 09:45:00','controllo'),
('MRNGPP48S30H501Y','MRRPLA68B20H501C','2022-01-20 08:00:00','urgenza');

INSERT INTO Esame (paziente, visita, tipo, data_esame, esito) VALUES
('RSSMRA80A01H501U',1,'ECG','2022-01-12 08:45:00','Sopraslivellamento ST in sede anteriore'),
('RSSMRA80A01H501U',1,'esami_ematochimici','2022-01-12 09:10:00','Troponina 4.2 ng/mL, elevata'),
('RSSMRA80A01H501U',1,'coronarografia','2022-01-12 11:00:00','Occlusione IVA prossimale'),
('BNCLRA65E41H501V',2,'ecocardiogramma','2022-01-13 10:30:00','FE 38%, ipocinesia diffusa'),
('GRECRL72T22H501W',3,'holter','2022-01-15 12:00:00','Fibrillazione atriale parossistica'),
('FNTNNA90M50H501X',4,'ECG','2022-01-18 10:00:00','Ritmo sinusale, nella norma'),
('MRNGPP48S30H501Y',5,'ecocardiogramma','2022-01-20 08:40:00','FE 30%, dilatazione ventricolare');

INSERT INTO Diagnosi (paziente, esame, nome, descrizione, data_diagnosi) VALUES
('BNCLRA65E41H501V',4,'Scompenso cardiaco','Disfunzione sistolica moderata','2022-01-13'),
('GRECRL72T22H501W',5,'Fibrillazione atriale','Forma parossistica','2022-01-15'),
('MRNGPP48S30H501Y',7,'Scompenso cardiaco','Disfunzione sistolica severa','2022-01-20');

INSERT INTO Cura (diagnosi, nome, tipo, data_inizio, data_fine) VALUES
(1,'Terapia con ACE-inibitore e betabloccante','farmacologica','2022-01-13',NULL),
(2,'Anticoagulazione orale','farmacologica','2022-01-15',NULL),
(3,'Impianto di resincronizzatore cardiaco','dispositivo','2022-01-22','2022-01-22');

INSERT INTO Rappresentante (id, nome, cognome, azienda, telefono) VALUES
(1,'Elena','Ricci','CardioPharma SpA','3401234567');

INSERT INTO Medicinale (codice_atc, nome, principio_attivo, rappresentante) VALUES
('C09AA05','Ramipril','Ramipril',1),
('C07AB07','Bisoprololo','Bisoprololo',1),
('B01AF01','Rivaroxaban','Rivaroxaban',1),
('C03CA01','Furosemide','Furosemide',1);

INSERT INTO Prescrizione (cura, medicinale, medico, posologia, data_prescr) VALUES
(1,'C09AA05','MRRPLA68B20H501C','5 mg once daily','2022-01-13'),
(1,'C07AB07','MRRPLA68B20H501C','2.5 mg once daily','2022-01-13'),
(2,'B01AF01','VRDGNN75F55H501B','20 mg once daily','2022-01-15'),
(3,'C03CA01','MRRPLA68B20H501C','25 mg twice daily','2022-01-20');

INSERT INTO Ordine (rappresentante, medicinale, quantita, data_ordine) VALUES
(1,'C09AA05',200,'2022-01-05'),
(1,'B01AF01',150,'2022-01-05');

INSERT INTO Referto (paziente, id_archivio, visita, esame, data_referto, contenuto) VALUES
('RSSMRA80A01H501U',1,1,3,'2022-01-12 12:00:00','Referto coronarografia: occlusione IVA prossimale, angioplastica con stent medicato'),
('BNCLRA65E41H501V',1,2,4,'2022-01-13 11:00:00','Referto ecocardiogramma: FE 38%, indicazione a terapia medica ottimizzata'),
('GRECRL72T22H501W',1,3,5,'2022-01-15 13:00:00','Referto holter: FA parossistica, indicazione ad anticoagulazione'),
('MRNGPP48S30H501Y',2,5,7,'2022-01-20 09:30:00','Referto ecocardiogramma: FE 30%, valutazione per CRT');

INSERT INTO Composizione (visita, specializzando, infermiere, oss) VALUES
(1,'FRRCLD95T10H501E','RSITRS85C45H501G','DMCLGU92L15H501J'),
(3,'GLLMRA96A41H501F','LMBGRG88E03H501H','SNTRTA87P51H501K');

INSERT INTO Esaminazione_Medico (esame, medico, ruolo) VALUES
(1,'BNCLCU70M12H501A','esecutore'),
(3,'BNCLCU70M12H501A','esecutore'),
(3,'BNCLCU70M12H501A','refertatore'),
(4,'CNTSRA82D44H501D','esecutore'),
(5,'VRDGNN75F55H501B','refertatore'),
(7,'CNTSRA82D44H501D','esecutore');
