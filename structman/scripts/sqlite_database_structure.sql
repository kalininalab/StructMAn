PRAGMA synchronous = OFF;
PRAGMA journal_mode = MEMORY;
BEGIN TRANSACTION;

CREATE TABLE `Alignment` (
  `Alignment_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Protein` integer  NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Structure` integer  NOT NULL REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Sequence_Identity` float DEFAULT NULL,
  `Coverage` float DEFAULT NULL,
  `Alignment` BLOB(65535) DEFAULT NULL,
  UNIQUE (`Alignment_Id`,`Protein`,`Structure`)
);

CREATE TABLE `Complex` (
  `Complex_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `PDB` varchar(16) NOT NULL,
  `Resolution` float DEFAULT NULL,
  `Chains` text,
  `Homooligomers` text,
  `Ligand_Profile` text,
  `Metal_Profile` text,
  `Ion_Profile` text,
  `Chain_Chain_Profile` text
);

CREATE TABLE `Interface` (
  `Interface_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Protein` integer  NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Structure_Recommendation` varchar(32) NOT NULL,
  `First_Position` integer  NOT NULL,
  `Last_Position` integer  NOT NULL,
  `Mean_Position` integer  NOT NULL,
  `Interface_Size` integer  NOT NULL
);

CREATE TABLE `Protein_Protein_Interaction` (
  `Protein_Protein_Interaction_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Interface_A` integer  NOT NULL REFERENCES `Interface` (`Interface_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Interface_B` integer  NOT NULL REFERENCES `Interface` (`Interface_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Complex` integer  NOT NULL REFERENCES `Complex` (`Complex_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Chain_A` char(1) DEFAULT NULL,
  `Chain_B` char(1) DEFAULT NULL,
  `Interaction_Score` float DEFAULT NULL
);

CREATE TABLE `Position_Position_Interaction` (
  `Position_Position_Interaction_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Position_A` integer NOT NULL REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Position_B` integer NOT NULL REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Residue_A` integer  NOT NULL REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Residue_B` integer  NOT NULL REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Interaction_Score` float DEFAULT NULL
);

CREATE TABLE `Protein` (
  `Protein_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Gene` integer  DEFAULT NULL REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Primary_Protein_Id` varchar(255) DEFAULT NULL,
  `Uniprot_Ac` varchar(16) DEFAULT NULL UNIQUE,
  `RefSeq_Ids` text,
  `Uniprot_Id` varchar(32) DEFAULT NULL,
  `Mutant_Type` varchar(32) DEFAULT NULL,
  `Mutant_Positions` text,
  `Original_Session` integer DEFAULT NULL,
  `Error_Code` integer DEFAULT NULL,
  `Error` varchar(255) DEFAULT NULL,
  `Sequence` text,
  `Species` varchar(255) DEFAULT NULL
);

CREATE TABLE `Gene` (
  `Gene_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Primary_Gene_Id` varchar(255) DEFAULT NULL,
  `Gene_Name` varchar(255) DEFAULT NULL,
  `Species` varchar(255) DEFAULT NULL
);

CREATE TABLE `GO_Term` (
  `GO_Term_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Name` tinytext NOT NULL,
  `Id` varchar(64) NOT NULL UNIQUE
);

CREATE TABLE `Indel` (
  `Indel_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Wildtype_Protein` integer  NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Mutant_Protein` integer REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Indel_Notation` text NOT NULL,
  `Analysis_Results` BLOB DEFAULT NULL
);

CREATE TABLE `Ligand` (
  `Ligand_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Name` varchar(12) NOT NULL,
  `Smiles` text NOT NULL,
  `Inchi` text NOT NULL
);

CREATE TABLE `Position` (
  `Position_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Position_Number` integer DEFAULT NULL,
  `Residue_Id` varchar(8) DEFAULT NULL,
  `Wildtype_Residue` char(1) DEFAULT NULL,
  `Protein` integer  DEFAULT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `IUPRED` float DEFAULT NULL,
  `IUPRED_Glob` varchar(16) DEFAULT NULL,
  `Recommended_Structure_Data` BLOB DEFAULT NULL,
  `Position_Data` BLOB DEFAULT NULL
);

CREATE TABLE `SNV` (
  `SNV_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Position` integer NOT NULL REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `New_AA` char(1) NOT NULL
);

CREATE TABLE `Multi_Mutation` (
  `Multi_Mutation_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `SNVs` varchar(256) DEFAULT NULL,
  `Indels` varchar(256) DEFAULT NULL,
  `Wildtype_Protein` integer  NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Mutant_Protein` integer  DEFAULT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE
);

CREATE TABLE `Pathway` (
  `Pathway_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Reactome_Id` varchar(32) NOT NULL,
  `Name` varchar(255) NOT NULL
);

CREATE TABLE `Residue` (
  `Residue_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Structure` integer  NOT NULL REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Number` varchar(32) NOT NULL,
  `Residue_Data` BLOB(65535) DEFAULT NULL
);

CREATE TABLE `RS_Protein_GO_Term` (
  `Protein` integer NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `GO_Term` integer NOT NULL REFERENCES `GO_Term` (`GO_Term_Id`) ON DELETE CASCADE ON UPDATE CASCADE
);

CREATE TABLE `RS_Isoform` (
  `Protein_A` integer  NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Protein_B` integer  NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Gene` integer  NOT NULL REFERENCES `Gene` (`Gene_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Session` integer NOT NULL REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Multi_Mutation` integer NOT NULL REFERENCES `Multi_Mutation` (`Multi_Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE
);

CREATE TABLE `RS_Protein_Pathway` (
  `Protein` integer  NOT NULL REFERENCES `Pathway` (`Pathway_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Pathway` integer NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE
);

CREATE TABLE `RS_Protein_Session` (
  `Protein` integer  NOT NULL REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Session` integer NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Input_Id` varchar(255) DEFAULT NULL,
  `Tags` mediumtext COLLATE BINARY
);

CREATE TABLE `RS_Position_Interface` (
  `Interface` integer  NOT NULL REFERENCES `Interface` (`Interface_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Position` integer NOT NULL REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Recommended_Residue` varchar(32) DEFAULT NULL
);

CREATE TABLE `RS_Residue_Interface` (
  `Residue` integer  NOT NULL REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Structure` integer  NOT NULL REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Interacting_Residue` BLOB(65535) DEFAULT NULL
);

CREATE TABLE `RS_Indel_Session` (
  `Session` integer NOT NULL REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Indel` integer  NOT NULL REFERENCES `Indel` (`Indel_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Tags` mediumtext COLLATE BINARY
);

CREATE TABLE `Database_Metadata` (
  `StructMAn_Version` varchar(32) NOT NULL,
  `PPI_Feature` integer NOT NULL
);

CREATE TABLE `RS_Multi_Mutation_Session` (
  `Session` integer NOT NULL REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Multi_Mutation` integer NOT NULL REFERENCES `Multi_Mutation` (`Multi_Mutation_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Tags` mediumtext COLLATE BINARY
);

CREATE TABLE `RS_Ligand_Structure` (
  `Ligand` integer NOT NULL REFERENCES `Ligand` (`Ligand_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Complex` integer  NOT NULL REFERENCES `Complex` (`Complex_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Chain` char(1) DEFAULT NULL,
  `Residue` varchar(16) DEFAULT NULL
);

CREATE TABLE `RS_Position_Session` (
  `Position` integer NOT NULL REFERENCES `Position` (`Position_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Session` integer NOT NULL REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Tag` mediumtext
);

CREATE TABLE `RS_SNV_Session` (
  `SNV` integer NOT NULL REFERENCES `SNV` (`SNV_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Session` integer NOT NULL REFERENCES `Session` (`Session_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Tag` mediumtext
);

CREATE TABLE `RS_Residue_Residue` (
  `Residue_1` integer  NOT NULL REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Structure_1` integer  NOT NULL REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Residue_2` integer  NOT NULL REFERENCES `Residue` (`Residue_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Structure_2` integer  NOT NULL REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Distance` float NOT NULL
);

CREATE TABLE `Session` (
  `Session_Id` integer NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Input_File` varchar(255) NOT NULL,
  `Checksum` integer  DEFAULT NULL,
  `Start` datetime NOT NULL,
  `End` datetime DEFAULT NULL
);

CREATE TABLE `Structure` (
  `Structure_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `PDB` varchar(16) NOT NULL,
  `Chain` char(1) NOT NULL,
  `Homooligomer` varchar(256) DEFAULT NULL
);

CREATE TABLE `Suggestion` (
  `Suggestion_Id` integer  NOT NULL PRIMARY KEY AUTOINCREMENT,
  `Protein` integer  NOT NULL REFERENCES `Protein` (`Protein_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Structure` integer  NOT NULL REFERENCES `Structure` (`Structure_Id`) ON DELETE CASCADE ON UPDATE CASCADE,
  `Method` varchar(32) DEFAULT NULL,
  `Score` integer  DEFAULT NULL,
  `Rank` integer  DEFAULT NULL
);


CREATE INDEX idx_uniprot_ac ON `Protein` (`Uniprot_Ac`);

CREATE INDEX idx_protein ON `Interface` (`Protein`);

CREATE INDEX idx_residue ON `Residue` (`Residue_Id`);

CREATE INDEX idx_pdb ON `Structure` (`PDB`);

END TRANSACTION;
