-- phpMyAdmin SQL Dump
-- version 4.9.1
-- https://www.phpmyadmin.net/
--
-- Host: wibi-guardian.helmholtz-hzi.de
-- Erstellungszeit: 25. Jan 2021 um 12:39
-- Server-Version: 10.1.48-MariaDB
-- PHP-Version: 7.0.33

SET SQL_MODE = "NO_AUTO_VALUE_ON_ZERO";
SET AUTOCOMMIT = 0;
START TRANSACTION;
SET time_zone = "+00:00";


/*!40101 SET @OLD_CHARACTER_SET_CLIENT=@@CHARACTER_SET_CLIENT */;
/*!40101 SET @OLD_CHARACTER_SET_RESULTS=@@CHARACTER_SET_RESULTS */;
/*!40101 SET @OLD_COLLATION_CONNECTION=@@COLLATION_CONNECTION */;
/*!40101 SET NAMES utf8mb4 */;

--
-- Datenbank: `agress_structman_test_3`
--
USE ; -- has to be replaced
-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Alignment`
--

CREATE TABLE `Alignment` (
  `Alignment_Id` int(11) UNSIGNED NOT NULL,
  `Protein` int(11) UNSIGNED NOT NULL,
  `Structure` int(11) UNSIGNED NOT NULL,
  `Sequence_Identity` float DEFAULT NULL,
  `Coverage` float DEFAULT NULL,
  `Alignment` BLOB(65535) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Complex`
--

CREATE TABLE `Complex` (
  `Complex_Id` int(11) UNSIGNED NOT NULL,   
  `PDB` varchar(16) NOT NULL,
  `Resolution` float DEFAULT NULL,
  `Chains` text,
  `Homooligomers` text,
  `Ligand_Profile` text,
  `Metal_Profile` text,
  `Ion_Profile` text,
  `Chain_Chain_Profile` text
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Interface`
--

CREATE TABLE `Interface` (
  `Interface_Id` int(11) UNSIGNED NOT NULL,
  `Protein` int(10) UNSIGNED NOT NULL,
  `Structure_Recommendation` varchar(32) NOT NULL,
  `First_Position` int(8) UNSIGNED NOT NULL,
  `Last_Position` int(8) UNSIGNED NOT NULL,
  `Mean_Position` int(8) UNSIGNED NOT NULL,
  `Interface_Size` int(8) UNSIGNED NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Protein_Protein_Interaction`
--

CREATE TABLE `Protein_Protein_Interaction` (
  `Protein_Protein_Interaction_Id` int(11) UNSIGNED NOT NULL,
  `Interface_A` int(10) UNSIGNED NOT NULL,
  `Interface_B` int(10) UNSIGNED NOT NULL,
  `Complex` int(11) UNSIGNED NOT NULL,
  `Chain_A` char(1) DEFAULT NULL,
  `Chain_B` char(1) DEFAULT NULL,
  `Interaction_Score` float DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Position_Position_Interaction`
--

CREATE TABLE `Position_Position_Interaction` (
  `Position_Position_Interaction_Id` int(11) NOT NULL,
  `Position_A` int(11) NOT NULL,
  `Position_B` int(11) NOT NULL,
  `Residue_A` int(11) UNSIGNED NOT NULL,
  `Residue_B` int(11) UNSIGNED NOT NULL,
  `Interaction_Score` float DEFAULT NULL  
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Protein`
--

CREATE TABLE `Protein` (
  `Protein_Id` int(10) UNSIGNED NOT NULL,
  `Gene` int(10) UNSIGNED DEFAULT NULL,
  `Primary_Protein_Id` varchar(255) DEFAULT NULL,
  `Uniprot_Ac` varchar(16) DEFAULT NULL,
  `RefSeq_Ids` text,
  `Uniprot_Id` varchar(32) DEFAULT NULL,
  `Mutant_Type` varchar(32) DEFAULT NULL,
  `Mutant_Positions` text,
  `Original_Session` int(11) DEFAULT NULL,
  `Error_Code` int(11) DEFAULT NULL,
  `Error` varchar(255) DEFAULT NULL,
  `Sequence` text,
  `Species` varchar(255) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Gene`
--

CREATE TABLE `Gene` (
  `Gene_Id` int(10) UNSIGNED NOT NULL,
  `Primary_Gene_Id` varchar(255) DEFAULT NULL,
  `Gene_Name` varchar(255) DEFAULT NULL,  
  `Species` varchar(255) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `GO_Term`
--

CREATE TABLE `GO_Term` (
  `GO_Term_Id` int(11) NOT NULL,
  `Name` tinytext NOT NULL,
  `Id` varchar(64) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Indel`
--

CREATE TABLE `Indel` (
  `Indel_Id` int(10) UNSIGNED NOT NULL,
  `Wildtype_Protein` int(10) UNSIGNED NOT NULL,
  `Mutant_Protein` int(10) UNSIGNED,
  `Indel_Notation` text COLLATE utf8_unicode_ci NOT NULL,
  `Analysis_Results` varbinary(8000) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Ligand`
--

CREATE TABLE `Ligand` (
  `Ligand_Id` int(11) NOT NULL,
  `Name` varchar(12) NOT NULL,
  `Smiles` text NOT NULL,
  `Inchi` text NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Position`
--

CREATE TABLE `Position` (
  `Position_Id` int(11) NOT NULL,
  `Position_Number` int(7) DEFAULT NULL,
  `Residue_Id` varchar(8) DEFAULT NULL,
  `Wildtype_Residue` char(1) DEFAULT NULL,
  `Protein` int(10) UNSIGNED DEFAULT NULL,
  `IUPRED` float DEFAULT NULL,
  `IUPRED_Glob` varchar(16) DEFAULT NULL,
  `Recommended_Structure_Data` varbinary(512) DEFAULT NULL,
  `Position_Data` varbinary(2048) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `SNV`
--

CREATE TABLE `SNV` (
  `SNV_Id` int(11) NOT NULL,
  `Position` int(11) NOT NULL,
  `New_AA` char(1) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Multi_Mutation`
--

CREATE TABLE `Multi_Mutation` (
  `Multi_Mutation_Id` int(11) NOT NULL,
  `SNVs` varchar(256) DEFAULT NULL,
  `Indels` varchar(256) DEFAULT NULL,
  `Wildtype_Protein` int(10) UNSIGNED NOT NULL,
  `Mutant_Protein` int(10) UNSIGNED DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Pathway`
--

CREATE TABLE `Pathway` (
  `Pathway_Id` int(11) NOT NULL,
  `Reactome_Id` varchar(32) NOT NULL,
  `Name` varchar(255) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Residue`
--

CREATE TABLE `Residue` (
  `Residue_Id` int(11) UNSIGNED NOT NULL,
  `Structure` int(11) UNSIGNED NOT NULL,
  `Number` varchar(32) NOT NULL,
  `Residue_Data` BLOB(65535) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Protein_GO_Term`
--

CREATE TABLE `RS_Protein_GO_Term` (
  `Protein` int(10) UNSIGNED NOT NULL,
  `GO_Term` int(11) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Isoform`
--

CREATE TABLE `RS_Isoform` (
  `Protein_A` int(10) UNSIGNED NOT NULL,
  `Protein_B` int(10) UNSIGNED NOT NULL,
  `Gene` int(10) UNSIGNED NOT NULL,
  `Session` int(11) NOT NULL,
  `Multi_Mutation` int(11) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Protein_Pathway`
--

CREATE TABLE `RS_Protein_Pathway` (
  `Protein` int(10) UNSIGNED NOT NULL,
  `Pathway` int(11) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Protein_Session`
--

CREATE TABLE `RS_Protein_Session` (
  `Protein` int(10) UNSIGNED NOT NULL,
  `Session` int(11) NOT NULL,
  `Input_Id` varchar(255) DEFAULT NULL,
  `Tags` mediumtext COLLATE utf8_unicode_ci
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Tabellenstruktur für Tabelle `RS_Position_Interface`
--

CREATE TABLE `RS_Position_Interface` (
  `Interface` int(10) UNSIGNED NOT NULL,
  `Position` int(11) NOT NULL,
  `Recommended_Residue` varchar(32) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Residue_Interface`
--

CREATE TABLE `RS_Residue_Interface` (
  `Residue` int(11) UNSIGNED NOT NULL,
  `Structure` int(11) UNSIGNED NOT NULL,
  `Interacting_Residue` BLOB(65535) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Indel_Session`
--

CREATE TABLE `RS_Indel_Session` (
  `Session` int(11) NOT NULL,
  `Indel` int(10) UNSIGNED NOT NULL,
  `Tags` mediumtext COLLATE utf8_unicode_ci
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

-- --------------------------------------------------------

--
-- Table structure for table `Database_Metadata`
--

CREATE TABLE `Database_Metadata` (
  `StructMAn_Version` varchar(32) COLLATE utf8_unicode_ci NOT NULL,
  `PPI_Feature` tinyint(1) NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

-- --------------------------------------------------------
--
-- Tabellenstruktur für Tabelle `RS_Multi_Mutation_Session`
--

CREATE TABLE `RS_Multi_Mutation_Session` (
  `Session` int(11) NOT NULL,
  `Multi_Mutation` int(11) NOT NULL,
  `Tags` mediumtext COLLATE utf8_unicode_ci
) ENGINE=InnoDB DEFAULT CHARSET=utf8 COLLATE=utf8_unicode_ci;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Ligand_Structure`
--

CREATE TABLE `RS_Ligand_Structure` (
  `Ligand` int(11) NOT NULL,
  `Complex` int(11) UNSIGNED NOT NULL,
  `Chain` char(1) DEFAULT NULL,
  `Residue` varchar(16) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Position_Session`
--

CREATE TABLE `RS_Position_Session` (
  `Position` int(11) NOT NULL,
  `Session` int(11) NOT NULL,
  `Tag` mediumtext
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_SNV_Session`
--

CREATE TABLE `RS_SNV_Session` (
  `SNV` int(11) NOT NULL,
  `Session` int(11) NOT NULL,
  `Tag` mediumtext
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `RS_Residue_Residue`
--

CREATE TABLE `RS_Residue_Residue` (
  `Residue_1` int(11) UNSIGNED NOT NULL,
  `Structure_1` int(11) UNSIGNED NOT NULL,
  `Residue_2` int(11) UNSIGNED NOT NULL,
  `Structure_2` int(11) UNSIGNED NOT NULL,
  `Distance` float NOT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Session`
--

CREATE TABLE `Session` (
  `Session_Id` int(11) NOT NULL,
  `Input_File` varchar(255) NOT NULL,
  `Checksum` INT(4) UNSIGNED DEFAULT NULL,
  `Start` datetime NOT NULL,
  `End` datetime DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Structure`
--

CREATE TABLE `Structure` (
  `Structure_Id` int(11) UNSIGNED NOT NULL,
  `PDB` varchar(16) NOT NULL,
  `Chain` char(1) NOT NULL,
  `Homooligomer` varchar(256) DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

-- --------------------------------------------------------

--
-- Tabellenstruktur für Tabelle `Suggestion`
--

CREATE TABLE `Suggestion` (
  `Suggestion_Id` int(11) UNSIGNED NOT NULL,
  `Protein` int(11) UNSIGNED NOT NULL,
  `Structure` int(11) UNSIGNED NOT NULL,
  `Method` varchar(32) DEFAULT NULL,
  `Score` int(4) UNSIGNED DEFAULT NULL,
  `Rank` int(2) UNSIGNED DEFAULT NULL
) ENGINE=InnoDB DEFAULT CHARSET=latin1;

--
-- Indizes der exportierten Tabellen
--

--
-- Indizes für die Tabelle `Alignment`
--
ALTER TABLE `Alignment`
  ADD PRIMARY KEY (`Alignment_Id`,`Protein`,`Structure`),
  ADD KEY `Protein` (`Protein`),
  ADD KEY `Structure` (`Structure`);

--
-- Indizes für die Tabelle `Complex`
--
ALTER TABLE `Complex`
  ADD PRIMARY KEY (`Complex_Id`);

--
-- Indizes für die Tabelle `Protein`
--
ALTER TABLE `Protein`
  ADD PRIMARY KEY (`Protein_Id`),
  ADD UNIQUE KEY `Name` (`Uniprot_Ac`,`Uniprot_Id`),
  ADD KEY `Gene` (`Gene`),
  ADD KEY `Uniprot_Ac` (`Uniprot_Ac`);

--
-- Indizes für die Tabelle `Gene`
--
ALTER TABLE `Gene`
  ADD PRIMARY KEY (`Gene_Id`);

--
-- Indizes für die Tabelle `Interface`
--
ALTER TABLE `Interface`
  ADD PRIMARY KEY (`Interface_Id`),
  ADD KEY `Protein` (`Protein`);

--
-- Indizes für die Tabelle `Protein_Protein_Interaction`
--
ALTER TABLE `Protein_Protein_Interaction`
  ADD PRIMARY KEY (`Protein_Protein_Interaction_Id`),
  ADD KEY `Interface_A` (`Interface_A`),
  ADD KEY `Interface_B` (`Interface_B`),
  ADD KEY `Complex` (`Complex`);


--
-- Indizes für die Tabelle `GO_Term`
--
ALTER TABLE `GO_Term`
  ADD PRIMARY KEY (`GO_Term_Id`),
  ADD UNIQUE KEY `Id` (`Id`) USING BTREE;

--
-- Indizes für die Tabelle `Indel`
--
ALTER TABLE `Indel`
  ADD PRIMARY KEY (`Indel_Id`),
  ADD KEY `Mutant_Protein` (`Mutant_Protein`),
  ADD KEY `Wildtype_Protein` (`Wildtype_Protein`);

--
-- Indizes für die Tabelle `Multi_Mutation`
--
ALTER TABLE `Multi_Mutation`
  ADD PRIMARY KEY (`Multi_Mutation_Id`),
  ADD KEY `Mutant_Protein` (`Mutant_Protein`),
  ADD KEY `Wildtype_Protein` (`Wildtype_Protein`);

--
-- Indizes für die Tabelle `SNV`
--
ALTER TABLE `SNV`
  ADD PRIMARY KEY (`SNV_Id`),
  ADD KEY `Position` (`Position`);

--
-- Indizes für die Tabelle `Ligand`
--
ALTER TABLE `Ligand`
  ADD PRIMARY KEY (`Ligand_Id`);

--
-- Indizes für die Tabelle `Position`
--
ALTER TABLE `Position`
  ADD PRIMARY KEY (`Position_Id`),
  ADD KEY `Protein` (`Protein`);

--
-- Indizes für die Tabelle `Position_Position_Interaction`
--
ALTER TABLE `Position_Position_Interaction`
  ADD PRIMARY KEY (`Position_Position_Interaction_Id`),
  ADD KEY `Position_A` (`Position_A`),
  ADD KEY `Position_B` (`Position_B`),
  ADD KEY `Residue_A` (`Residue_A`),
  ADD KEY `Residue_B` (`Residue_B`);

--
-- Indizes für die Tabelle `Pathway`
--
ALTER TABLE `Pathway`
  ADD PRIMARY KEY (`Pathway_Id`);

--
-- Indizes für die Tabelle `Residue`
--
ALTER TABLE `Residue`
  ADD PRIMARY KEY (`Residue_Id`,`Structure`),
  ADD KEY `Structure` (`Structure`),
  ADD KEY `Residue_Id` (`Residue_Id`);

--
-- Indizes für die Tabelle `RS_Protein_GO_Term`
--
ALTER TABLE `RS_Protein_GO_Term`
  ADD KEY `Protein` (`Protein`),
  ADD KEY `GO_Term` (`GO_Term`);

--
-- Indizes für die Tabelle `RS_Isoform`
--
ALTER TABLE `RS_Isoform`
  ADD KEY `Protein_A` (`Protein_A`),
  ADD KEY `Protein_B` (`Protein_B`),
  ADD KEY `Gene` (`Gene`),
  ADD KEY `Session` (`Session`),
  ADD KEY `Multi_Mutation` (`Multi_Mutation`);

--
-- Indizes für die Tabelle `RS_Protein_Pathway`
--
ALTER TABLE `RS_Protein_Pathway`
  ADD KEY `Protein` (`Protein`),
  ADD KEY `Pathway` (`Pathway`);

--
-- Indizes für die Tabelle `RS_Protein_Session`
--
ALTER TABLE `RS_Protein_Session`
  ADD KEY `Protein` (`Protein`,`Session`),
  ADD KEY `Session` (`Session`);

--
-- Indizes für die Tabelle `RS_Indel_Session`
--
ALTER TABLE `RS_Indel_Session`
  ADD KEY `Indel` (`Indel`),
  ADD KEY `Session` (`Session`);

--
-- Indizes für die Tabelle `RS_Multi_Mutation_Session`
--
ALTER TABLE `RS_Multi_Mutation_Session`
  ADD KEY `Multi_Mutation` (`Multi_Mutation`),
  ADD KEY `Session` (`Session`);

--
-- Indizes für die Tabelle `RS_Ligand_Structure`
--
ALTER TABLE `RS_Ligand_Structure`
  ADD KEY `Ligand` (`Ligand`,`Complex`),
  ADD KEY `RS_Ligand_Structure_ibfk_2` (`Complex`);

--
-- Indizes für die Tabelle `RS_Position_Session`
--
ALTER TABLE `RS_Position_Session`
  ADD KEY `Position` (`Position`),
  ADD KEY `Session` (`Session`);

--
-- Indizes für die Tabelle `RS_Position_Interface`
--
ALTER TABLE `RS_Position_Interface`
  ADD KEY `Position` (`Position`),
  ADD KEY `Interface` (`Interface`);

--
-- Indizes für die Tabelle `RS_Residue_Interface`
--
ALTER TABLE `RS_Residue_Interface`
  ADD KEY `Residue` (`Residue`),
  ADD KEY `Structure` (`Structure`);


--
-- Indizes für die Tabelle `RS_SNV_Session`
--
ALTER TABLE `RS_SNV_Session`
  ADD KEY `SNV` (`SNV`),
  ADD KEY `Session` (`Session`);

--
-- Indizes für die Tabelle `RS_Residue_Residue`
--
ALTER TABLE `RS_Residue_Residue`
  ADD KEY `Residue_1` (`Residue_1`),
  ADD KEY `Structure_1` (`Structure_1`),
  ADD KEY `Residue_2` (`Residue_2`),
  ADD KEY `Structure_2` (`Structure_2`);

--
-- Indizes für die Tabelle `Session`
--
ALTER TABLE `Session`
  ADD PRIMARY KEY (`Session_Id`);

--
-- Indizes für die Tabelle `Structure`
--
ALTER TABLE `Structure`
  ADD PRIMARY KEY (`Structure_Id`),
  ADD KEY `PDB` (`PDB`);

--
-- Indizes für die Tabelle `Suggestion`
--
ALTER TABLE `Suggestion`
  ADD PRIMARY KEY (`Suggestion_Id`),
  ADD KEY `Protein` (`Protein`),
  ADD KEY `Structure` (`Structure`);

--
-- AUTO_INCREMENT für exportierte Tabellen
--

--
-- AUTO_INCREMENT für Tabelle `Alignment`
--
ALTER TABLE `Alignment`
  MODIFY `Alignment_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Complex`
--
ALTER TABLE `Complex`
  MODIFY `Complex_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Protein`
--
ALTER TABLE `Protein`
  MODIFY `Protein_Id` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Gene`
--
ALTER TABLE `Gene`
  MODIFY `Gene_Id` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Interface`
--
ALTER TABLE `Interface`
  MODIFY `Interface_Id` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Protein_Protein_Interaction`
--
ALTER TABLE `Protein_Protein_Interaction`
  MODIFY `Protein_Protein_Interaction_Id` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `GO_Term`
--
ALTER TABLE `GO_Term`
  MODIFY `GO_Term_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Indel`
--
ALTER TABLE `Indel`
  MODIFY `Indel_Id` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Multi_Mutation`
--
ALTER TABLE `Multi_Mutation`
  MODIFY `Multi_Mutation_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Ligand`
--
ALTER TABLE `Ligand`
  MODIFY `Ligand_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Position`
--
ALTER TABLE `Position`
  MODIFY `Position_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Position_Position_Interaction`
--
ALTER TABLE `Position_Position_Interaction`
  MODIFY `Position_Position_Interaction_Id` int(10) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `SNV`
--
ALTER TABLE `SNV`
  MODIFY `SNV_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Pathway`
--
ALTER TABLE `Pathway`
  MODIFY `Pathway_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Residue`
--
ALTER TABLE `Residue`
  MODIFY `Residue_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Session`
--
ALTER TABLE `Session`
  MODIFY `Session_Id` int(11) NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Structure`
--
ALTER TABLE `Structure`
  MODIFY `Structure_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- AUTO_INCREMENT für Tabelle `Suggestion`
--
ALTER TABLE `Suggestion`
  MODIFY `Suggestion_Id` int(11) UNSIGNED NOT NULL AUTO_INCREMENT;

--
-- Constraints der exportierten Tabellen
--

--
-- Constraints der Tabelle `Alignment`
--
ALTER TABLE `Alignment`
  ADD CONSTRAINT `Alignment_ibfk_1` FOREIGN KEY (`Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Alignment_ibfk_2` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Protein`
--
ALTER TABLE `Protein`
  ADD CONSTRAINT `Protein_ibfk_1` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE;


--
-- Constraints der Tabelle `Protein_Protein_Interaction`
--
ALTER TABLE `Protein_Protein_Interaction`
  ADD CONSTRAINT `Interface_A_ibfk_1` FOREIGN KEY (`Interface_A`) REFERENCES `Interface` (`Interface_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Interface_B_ibfk_2` FOREIGN KEY (`Interface_B`) REFERENCES `Interface` (`Interface_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Complex_ibfk_3` FOREIGN KEY (`Complex`) REFERENCES `Complex` (`Complex_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Indel`
--
ALTER TABLE `Indel`
  ADD CONSTRAINT `Indel_ibfk_1` FOREIGN KEY (`Mutant_Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Indel_ibfk_2` FOREIGN KEY (`Wildtype_Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Multi_Mutation`
--
ALTER TABLE `Multi_Mutation`
  ADD CONSTRAINT `Multi_Mutation_ibfk_1` FOREIGN KEY (`Mutant_Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Multi_Mutation_ibfk_2` FOREIGN KEY (`Wildtype_Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Position`
--
ALTER TABLE `Position`
  ADD CONSTRAINT `Position_ibfk_1` FOREIGN KEY (`Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Interface`
--
ALTER TABLE `Interface`
  ADD CONSTRAINT `Interface_ibfk_1` FOREIGN KEY (`Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Position_Position_Interaction`
--
ALTER TABLE `Position_Position_Interaction`
  ADD CONSTRAINT `Position_A_ibfk_1` FOREIGN KEY (`Position_A`) REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Position_B_ibfk_2` FOREIGN KEY (`Position_B`) REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Residue_A_ibfk_1` FOREIGN KEY (`Residue_A`) REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Residue_B_ibfk_2` FOREIGN KEY (`Residue_B`) REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `SNV`
--
ALTER TABLE `SNV`
  ADD CONSTRAINT `SNV_ibfk_1` FOREIGN KEY (`Position`) REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Suggestion`
--
ALTER TABLE `Suggestion`
  ADD CONSTRAINT `Suggestion_ibfk_1` FOREIGN KEY (`Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `Suggestion_ibfk_2` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `Residue`
--
ALTER TABLE `Residue`
  ADD CONSTRAINT `Residue_ibfk_1` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Protein_GO_Term`
--
ALTER TABLE `RS_Protein_GO_Term`
  ADD CONSTRAINT `RS_Protein_GO_Term_ibfk_1` FOREIGN KEY (`Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Protein_GO_Term_ibfk_2` FOREIGN KEY (`GO_Term`) REFERENCES `GO_Term` (`GO_Term_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Isoform`
--
ALTER TABLE `RS_Isoform`
  ADD CONSTRAINT `RS_Isoform_ibfk_1` FOREIGN KEY (`Protein_A`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Isoform_ibfk_2` FOREIGN KEY (`Protein_B`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Isoform_ibfk_3` FOREIGN KEY (`Gene`) REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Isoform_ibfk_4` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Isoform_ibfk_5` FOREIGN KEY (`Multi_Mutation`) REFERENCES `Multi_Mutation` (`Multi_Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE;


--
-- Constraints der Tabelle `RS_Protein_Pathway`
--
ALTER TABLE `RS_Protein_Pathway`
  ADD CONSTRAINT `RS_Protein_Pathway_ibfk_1` FOREIGN KEY (`Pathway`) REFERENCES `Pathway` (`Pathway_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Protein_Pathway_ibfk_2` FOREIGN KEY (`Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Protein_Session`
--
ALTER TABLE `RS_Protein_Session`
  ADD CONSTRAINT `RS_Protein_Session_ibfk_1` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Protein_Session_ibfk_2` FOREIGN KEY (`Protein`) REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Indel_Session`
--
ALTER TABLE `RS_Indel_Session`
  ADD CONSTRAINT `RS_Indel_Session_ibfk_1` FOREIGN KEY (`Indel`) REFERENCES `Indel` (`Indel_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Indel_Session_ibfk_2` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Indel_Session`
--
ALTER TABLE `RS_Multi_Mutation_Session`
  ADD CONSTRAINT `RS_Multi_Mutation_Session_ibfk_1` FOREIGN KEY (`Multi_Mutation`) REFERENCES `Multi_Mutation` (`Multi_Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Multi_Mutation_Session_ibfk_2` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Ligand_Structure`
--
ALTER TABLE `RS_Ligand_Structure`
  ADD CONSTRAINT `RS_Ligand_Structure_ibfk_1` FOREIGN KEY (`Ligand`) REFERENCES `Ligand` (`Ligand_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Ligand_Structure_ibfk_2` FOREIGN KEY (`Complex`) REFERENCES `Complex` (`Complex_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Position_Session`
--
ALTER TABLE `RS_Position_Session`
  ADD CONSTRAINT `RS_Position_Session_ibfk_1` FOREIGN KEY (`Position`) REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Position_Session_ibfk_2` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Position_Interface`
--
ALTER TABLE `RS_Position_Interface`
  ADD CONSTRAINT `RS_Position_Interface_ibfk_1` FOREIGN KEY (`Position`) REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Position_Interface_ibfk_2` FOREIGN KEY (`Interface`) REFERENCES `Interface` (`Interface_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Residue_Interface`
--
ALTER TABLE `RS_Residue_Interface`
  ADD CONSTRAINT `RS_Residue_Interface_ibfk_1` FOREIGN KEY (`Residue`) REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Residue_Interface_ibfk_2` FOREIGN KEY (`Structure`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_SNV_Session`
--
ALTER TABLE `RS_SNV_Session`
  ADD CONSTRAINT `RS_SNV_Session_ibfk_1` FOREIGN KEY (`SNV`) REFERENCES `SNV` (`SNV_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_SNV_Session_ibfk_2` FOREIGN KEY (`Session`) REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE;

--
-- Constraints der Tabelle `RS_Residue_Residue`
--
ALTER TABLE `RS_Residue_Residue`
  ADD CONSTRAINT `RS_Residue_Residue_ibfk_1` FOREIGN KEY (`Residue_1`) REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Residue_Residue_ibfk_2` FOREIGN KEY (`Structure_1`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Residue_Residue_ibfk_3` FOREIGN KEY (`Residue_2`) REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  ADD CONSTRAINT `RS_Residue_Residue_ibfk_4` FOREIGN KEY (`Structure_2`) REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE;
COMMIT;

/*!40101 SET CHARACTER_SET_CLIENT=@OLD_CHARACTER_SET_CLIENT */;
/*!40101 SET CHARACTER_SET_RESULTS=@OLD_CHARACTER_SET_RESULTS */;
/*!40101 SET COLLATION_CONNECTION=@OLD_COLLATION_CONNECTION */;
