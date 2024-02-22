-- phpMyAdmin SQL Dump
-- version 4.9.1
-- https://www.phpmyadmin.net/
--
-- Host: wibi-guardian.helmholtz-hzi.de
-- Erstellungszeit: 11. Okt 2021 um 13:45
-- Server-Version: 10.1.48-MariaDB
-- PHP-Version: 7.0.33

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET AUTOCOMMIT = 0;
START TRANSACTION;
SET time_zone = "+00:00";

USE ;
/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;

--
-- Datenbank: `agress_id_mapping`
--

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `UNIPROT`
--

CREATE TABLE `UNIPROT` (
  `Uniprot_Ac` varchar(16) COLLATE utf8_unicode_ci NOT NULL,
  `Uniprot_Id` varchar(32) COLLATE utf8_unicode_ci DEFAULT NULL,
  `RefSeq` varchar(16) COLLATE utf8_unicode_ci DEFAULT NULL,
  `RefSeq_NT` varchar(16) COLLATE utf8_unicode_ci DEFAULT NULL,
  `Sequence` BLOB(65535) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

--
-- Tabellenstruktur für Tabelle `EMBL_Transcripts`
--

CREATE TABLE `EMBL_Transcripts` (
  `Transcript_ID` varchar(32) COLLATE utf8_unicode_ci NOT NULL,
  `Gene_ID` varchar(32) COLLATE utf8_unicode_ci NOT NULL,  
  `Sequence` blob
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

-- --------------------------------------------------------
--
-- Indizes der exportierten Tabellen
--

--
-- Indizes für die Tabelle `UNIPROT`
--
ALTER TABLE `UNIPROT`
  ADD PRIMARY KEY (`Uniprot_Ac`),
  ADD KEY `Uniprot_Id` (`Uniprot_Id`),
  ADD KEY `RefSeq` (`RefSeq`),
  ADD KEY `RefSeq_NT` (`RefSeq_NT`);


--
-- Indizes für die Tabelle `EMBL_Transcripts`
--
ALTER TABLE `EMBL_Transcripts`
  ADD PRIMARY KEY (`Transcript_ID`);


COMMIT;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
