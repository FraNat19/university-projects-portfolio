-- Cardiology ward database (db62)
-- Physical schema derived from the E-R and logical models in report_conclusivo_DB.pdf
-- Target: MySQL 8.0, InnoDB, utf8mb4

DROP DATABASE IF EXISTS db62;
CREATE DATABASE db62 CHARACTER SET utf8mb4 COLLATE utf8mb4_unicode_ci;
USE db62;

-- ---------------------------------------------------------------- structure

CREATE TABLE Stanza (
    codice        CHAR(6)      PRIMARY KEY,
    piano         TINYINT      NOT NULL,
    tipo          ENUM('degenza','preospedalizzazione','osservazione') NOT NULL DEFAULT 'degenza'
) ENGINE=InnoDB;

CREATE TABLE Letto (
    codice        CHAR(8)      PRIMARY KEY,
    stanza        CHAR(6)      NOT NULL,
    occupato      BOOLEAN      NOT NULL DEFAULT FALSE,
    CONSTRAINT fk_letto_stanza FOREIGN KEY (stanza) REFERENCES Stanza(codice)
        ON DELETE RESTRICT ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE INDEX idx_letto_stanza ON Letto(stanza);

CREATE TABLE Turno (
    id            TINYINT      PRIMARY KEY,
    nome          ENUM('mattina','pomeriggio','notte','riposo') NOT NULL UNIQUE,
    ora_inizio    TIME         NULL,
    ora_fine      TIME         NULL
) ENGINE=InnoDB;

-- ---------------------------------------------------------------- staff

CREATE TABLE Medico (
    cf            CHAR(16)     PRIMARY KEY,
    nome          VARCHAR(50)  NOT NULL,
    cognome       VARCHAR(50)  NOT NULL,
    specializzazione VARCHAR(80) NOT NULL,
    strutturato   BOOLEAN      NOT NULL DEFAULT TRUE
) ENGINE=InnoDB;

CREATE TABLE Specializzando (
    cf            CHAR(16)     PRIMARY KEY,
    nome          VARCHAR(50)  NOT NULL,
    cognome       VARCHAR(50)  NOT NULL,
    anno_corso    TINYINT      NOT NULL CHECK (anno_corso BETWEEN 1 AND 5),
    tutor         CHAR(16)     NULL,
    CONSTRAINT fk_spec_tutor FOREIGN KEY (tutor) REFERENCES Medico(cf)
        ON DELETE SET NULL ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE TABLE Infermiere (
    cf            CHAR(16)     PRIMARY KEY,
    nome          VARCHAR(50)  NOT NULL,
    cognome       VARCHAR(50)  NOT NULL,
    turno         TINYINT      NOT NULL,
    CONSTRAINT fk_inf_turno FOREIGN KEY (turno) REFERENCES Turno(id)
        ON DELETE RESTRICT ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE TABLE OSS (
    cf            CHAR(16)     PRIMARY KEY,
    nome          VARCHAR(50)  NOT NULL,
    cognome       VARCHAR(50)  NOT NULL,
    turno         TINYINT      NOT NULL,
    CONSTRAINT fk_oss_turno FOREIGN KEY (turno) REFERENCES Turno(id)
        ON DELETE RESTRICT ON UPDATE CASCADE
) ENGINE=InnoDB;

-- ---------------------------------------------------------------- patients

CREATE TABLE Archivio (
    id            INT          PRIMARY KEY,
    descrizione   VARCHAR(120) NOT NULL,
    data_apertura DATE         NOT NULL
) ENGINE=InnoDB;

CREATE TABLE Paziente (
    cf            CHAR(16)     PRIMARY KEY,
    nome          VARCHAR(50)  NOT NULL,
    cognome       VARCHAR(50)  NOT NULL,
    data_nascita  DATE         NOT NULL,
    sesso         ENUM('M','F') NOT NULL,
    telefono      VARCHAR(20)  NULL,
    id_archivio   INT          NOT NULL,
    letto         CHAR(8)      NULL,
    CONSTRAINT fk_paz_archivio FOREIGN KEY (id_archivio) REFERENCES Archivio(id)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT fk_paz_letto FOREIGN KEY (letto) REFERENCES Letto(codice)
        ON DELETE SET NULL ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE INDEX idx_paz_archivio ON Paziente(id_archivio);
CREATE INDEX idx_paz_cognome  ON Paziente(cognome, nome);

CREATE TABLE Sintomo (
    id            INT AUTO_INCREMENT PRIMARY KEY,
    paziente      CHAR(16)     NOT NULL,
    descrizione   VARCHAR(120) NOT NULL,
    data_comparsa DATE         NOT NULL,
    codice_triage ENUM('bianco','verde','giallo','arancione','rosso') NOT NULL,
    CONSTRAINT fk_sint_paziente FOREIGN KEY (paziente) REFERENCES Paziente(cf)
        ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE INDEX idx_sint_paziente ON Sintomo(paziente);
CREATE INDEX idx_sint_descr    ON Sintomo(descrizione);

-- ---------------------------------------------------------------- clinical path

CREATE TABLE Visita (
    id            INT AUTO_INCREMENT PRIMARY KEY,
    paziente      CHAR(16)     NOT NULL,
    medico        CHAR(16)     NOT NULL,
    data_visita   DATETIME     NOT NULL,
    tipo          ENUM('prima_visita','controllo','follow_up','urgenza') NOT NULL,
    CONSTRAINT fk_vis_paziente FOREIGN KEY (paziente) REFERENCES Paziente(cf)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT fk_vis_medico FOREIGN KEY (medico) REFERENCES Medico(cf)
        ON DELETE RESTRICT ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE INDEX idx_vis_paziente ON Visita(paziente, data_visita);
CREATE INDEX idx_vis_medico   ON Visita(medico);

CREATE TABLE Esame (
    id            INT AUTO_INCREMENT PRIMARY KEY,
    paziente      CHAR(16)     NOT NULL,
    visita        INT          NULL,
    tipo          ENUM('ECG','ecocardiogramma','holter','RX_torace',
                       'coronarografia','esami_ematochimici','angioplastica') NOT NULL,
    data_esame    DATETIME     NOT NULL,
    esito         TEXT         NULL,
    CONSTRAINT fk_esa_paziente FOREIGN KEY (paziente) REFERENCES Paziente(cf)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT fk_esa_visita FOREIGN KEY (visita) REFERENCES Visita(id)
        ON DELETE SET NULL ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE INDEX idx_esa_paziente ON Esame(paziente, data_esame);
CREATE INDEX idx_esa_tipo     ON Esame(tipo);

CREATE TABLE Diagnosi (
    id            INT AUTO_INCREMENT PRIMARY KEY,
    paziente      CHAR(16)     NOT NULL,
    esame         INT          NULL,
    nome          VARCHAR(120) NOT NULL,
    descrizione   TEXT         NULL,
    data_diagnosi DATE         NOT NULL,
    CONSTRAINT fk_dia_paziente FOREIGN KEY (paziente) REFERENCES Paziente(cf)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT fk_dia_esame FOREIGN KEY (esame) REFERENCES Esame(id)
        ON DELETE SET NULL ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE INDEX idx_dia_paziente ON Diagnosi(paziente);
CREATE INDEX idx_dia_nome     ON Diagnosi(nome);

CREATE TABLE Cura (
    id            INT AUTO_INCREMENT PRIMARY KEY,
    diagnosi      INT          NOT NULL,
    nome          VARCHAR(120) NOT NULL,
    tipo          ENUM('farmacologica','interventistica','dispositivo') NOT NULL,
    data_inizio   DATE         NOT NULL,
    data_fine     DATE         NULL,
    CONSTRAINT fk_cura_diagnosi FOREIGN KEY (diagnosi) REFERENCES Diagnosi(id)
        ON DELETE CASCADE ON UPDATE CASCADE,
    CONSTRAINT chk_cura_periodo CHECK (data_fine IS NULL OR data_fine >= data_inizio)
) ENGINE=InnoDB;

CREATE INDEX idx_cura_diagnosi ON Cura(diagnosi);

-- ---------------------------------------------------------------- documents

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
    CONSTRAINT fk_ref_archivio FOREIGN KEY (id_archivio) REFERENCES Archivio(id)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT fk_ref_visita FOREIGN KEY (visita) REFERENCES Visita(id)
        ON DELETE SET NULL ON UPDATE CASCADE,
    CONSTRAINT fk_ref_esame FOREIGN KEY (esame) REFERENCES Esame(id)
        ON DELETE SET NULL ON UPDATE CASCADE,
    CONSTRAINT chk_ref_origine CHECK (visita IS NOT NULL OR esame IS NOT NULL)
) ENGINE=InnoDB;

CREATE INDEX idx_ref_paziente ON Referto(paziente, data_referto);
CREATE INDEX idx_ref_archivio ON Referto(id_archivio);

-- ---------------------------------------------------------------- drugs

CREATE TABLE Rappresentante (
    id            INT          PRIMARY KEY,
    nome          VARCHAR(50)  NOT NULL,
    cognome       VARCHAR(50)  NOT NULL,
    azienda       VARCHAR(80)  NOT NULL,
    telefono      VARCHAR(20)  NULL
) ENGINE=InnoDB;

CREATE TABLE Medicinale (
    codice_atc    VARCHAR(10)  PRIMARY KEY,
    nome          VARCHAR(80)  NOT NULL,
    principio_attivo VARCHAR(80) NOT NULL,
    rappresentante INT         NULL,
    CONSTRAINT fk_med_rappr FOREIGN KEY (rappresentante) REFERENCES Rappresentante(id)
        ON DELETE SET NULL ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE TABLE Prescrizione (
    id            INT AUTO_INCREMENT PRIMARY KEY,
    cura          INT          NOT NULL,
    medicinale    VARCHAR(10)  NOT NULL,
    medico        CHAR(16)     NOT NULL,
    posologia     VARCHAR(80)  NOT NULL,
    data_prescr   DATE         NOT NULL,
    CONSTRAINT fk_pre_cura FOREIGN KEY (cura) REFERENCES Cura(id)
        ON DELETE CASCADE ON UPDATE CASCADE,
    CONSTRAINT fk_pre_medicinale FOREIGN KEY (medicinale) REFERENCES Medicinale(codice_atc)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT fk_pre_medico FOREIGN KEY (medico) REFERENCES Medico(cf)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT uq_pre UNIQUE (cura, medicinale)
) ENGINE=InnoDB;

CREATE INDEX idx_pre_medico     ON Prescrizione(medico);
CREATE INDEX idx_pre_medicinale ON Prescrizione(medicinale);

CREATE TABLE Ordine (
    id            INT AUTO_INCREMENT PRIMARY KEY,
    rappresentante INT         NOT NULL,
    medicinale    VARCHAR(10)  NOT NULL,
    quantita      INT          NOT NULL CHECK (quantita > 0),
    data_ordine   DATE         NOT NULL,
    CONSTRAINT fk_ord_rappr FOREIGN KEY (rappresentante) REFERENCES Rappresentante(id)
        ON DELETE RESTRICT ON UPDATE CASCADE,
    CONSTRAINT fk_ord_medicinale FOREIGN KEY (medicinale) REFERENCES Medicinale(codice_atc)
        ON DELETE RESTRICT ON UPDATE CASCADE
) ENGINE=InnoDB;

-- ---------------------------------------------------------------- teams

CREATE TABLE Composizione (
    visita        INT          NOT NULL,
    specializzando CHAR(16)    NULL,
    infermiere    CHAR(16)     NULL,
    oss           CHAR(16)     NULL,
    PRIMARY KEY (visita, specializzando, infermiere, oss),
    CONSTRAINT fk_comp_visita FOREIGN KEY (visita) REFERENCES Visita(id)
        ON DELETE CASCADE ON UPDATE CASCADE,
    CONSTRAINT fk_comp_spec FOREIGN KEY (specializzando) REFERENCES Specializzando(cf)
        ON DELETE CASCADE ON UPDATE CASCADE,
    CONSTRAINT fk_comp_inf FOREIGN KEY (infermiere) REFERENCES Infermiere(cf)
        ON DELETE CASCADE ON UPDATE CASCADE,
    CONSTRAINT fk_comp_oss FOREIGN KEY (oss) REFERENCES OSS(cf)
        ON DELETE CASCADE ON UPDATE CASCADE
) ENGINE=InnoDB;

CREATE TABLE Esaminazione_Medico (
    esame         INT          NOT NULL,
    medico        CHAR(16)     NOT NULL,
    ruolo         ENUM('esecutore','refertatore') NOT NULL,
    PRIMARY KEY (esame, medico, ruolo),
    CONSTRAINT fk_esm_esame FOREIGN KEY (esame) REFERENCES Esame(id)
        ON DELETE CASCADE ON UPDATE CASCADE,
    CONSTRAINT fk_esm_medico FOREIGN KEY (medico) REFERENCES Medico(cf)
        ON DELETE RESTRICT ON UPDATE CASCADE
) ENGINE=InnoDB;
