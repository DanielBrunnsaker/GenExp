#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 10:15:17 2025

@author: danbru
"""

import sqlite3
import pandas as pd
import argparse

def create_ion_mobility_db(csv_file, imsdb_file, sep = ','):
    """
    Create a Skyline 22.2+ Ion Mobility Database (.imsdb) with the 
    required tables: molecules, ions, and ionMobilityValues, 
    storing CCS values (IonMobilityUnits=2).
    
    The CSV must contain columns: Molecule, Formula, Adduct, PrecursorMz, CCS
    """
    df = pd.read_csv(csv_file,sep = sep)

    required_cols = {"Molecule", "Formula", "Adduct", "PrecursorMz", "CCS"}
    if not required_cols.issubset(df.columns):
        raise ValueError(f"CSV must contain columns: {required_cols}")

    # Connect (or create) the SQLite DB
    conn = sqlite3.connect(imsdb_file)
    c = conn.cursor()

    # Drop old tables if they exist, so we start clean
    c.execute("DROP TABLE IF EXISTS ionMobilityValues")
    c.execute("DROP TABLE IF EXISTS ions")
    c.execute("DROP TABLE IF EXISTS molecules")

    # Create 'molecules' table
    #
    # Skyline looks for columns:
    #   id (PK)
    #   PeptideModifiedSequence (TEXT)
    #   moleculeName (TEXT)
    #   chemicalFormula (TEXT)
    #   inchiKey (TEXT)
    #   otherKeys (TEXT)
    c.execute("""
        CREATE TABLE molecules (
            id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
            PeptideModifiedSequence TEXT,
            moleculeName TEXT,
            chemicalFormula TEXT,
            inchiKey TEXT,
            otherKeys TEXT
        )
    """)

    # Create 'ions' table
    # Skyline expects columns:
    #   id (PK)
    #   moleculeId (FK into molecules.id)
    #   precursorAdduct (TEXT)
    c.execute("""
        CREATE TABLE ions (
            id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
            moleculeId INTEGER NOT NULL,
            precursorAdduct TEXT,
            FOREIGN KEY (moleculeId) REFERENCES molecules(id)
        )
    """)

    # Create 'ionMobilityValues' table
    # The critical part: IonMobilityUnits is an integer enum, not TEXT
    #   0 => None/Unknown
    #   1 => Drift time (msec)
    #   2 => CCS (square Angstroms)
    #   3 => Inverse K0
    #
    # Skyline reads:
    #   CollisionalCrossSectionSqA => your CCS
    #   IonMobility => your drift time if IonMobilityUnits=1 (or None if CCS-only)
    #   IonMobilityUnits => integer code
    #   HighEnergyIonMobilityOffset => additional offset for fragments (use 0 if unknown)
    c.execute("""
        CREATE TABLE ionMobilityValues (
            id INTEGER NOT NULL PRIMARY KEY AUTOINCREMENT,
            ionId INTEGER NOT NULL,
            CollisionalCrossSectionSqA REAL,
            IonMobility REAL,
            IonMobilityUnits INTEGER,
            HighEnergyIonMobilityOffset REAL,
            FOREIGN KEY (ionId) REFERENCES ions(id)
        )
    """)

    # Insert rows from your CSV
    for _, row in df.iterrows():
        # Insert into molecules
        c.execute("""
            INSERT INTO molecules 
                (PeptideModifiedSequence, moleculeName, chemicalFormula, inchiKey, otherKeys)
            VALUES (?, ?, ?, ?, ?)
        """, (
            "",  # PeptideModifiedSequence => empty for small molecules
            row["Molecule"],
            row["Formula"],
            "",               # InChIKey => unknown
            str(row["PrecursorMz"])  # store precursor m/z in 'otherKeys' for reference
        ))
        molecule_id = c.lastrowid

        # Insert into ions
        c.execute("""
            INSERT INTO ions (moleculeId, precursorAdduct)
            VALUES (?, ?)
        """, (molecule_id, row["Adduct"]))
        ion_id = c.lastrowid

        # Insert into ionMobilityValues
        # For CCS data: IonMobility=NULL, IonMobilityUnits=2
        c.execute("""
            INSERT INTO ionMobilityValues
                (ionId, CollisionalCrossSectionSqA, IonMobility, IonMobilityUnits, HighEnergyIonMobilityOffset)
            VALUES (?, ?, ?, ?, ?)
        """, (
            ion_id,
            row["CCS"],    # CollisionalCrossSectionSqA
            None,          # IonMobility => None (we only have CCS)
            2,             # IonMobilityUnits => 2 => CCS
            0.0            # HighEnergyIonMobilityOffset => 0 if unknown
        ))

    conn.commit()
    conn.close()
    print(f"\nCreated Skyline Ion Mobility DB: {imsdb_file}\n")



if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert CSV to Skyline-compatible .imsdb file.")
    parser.add_argument("csv_file", help="Input CSV file with Molecule, Formula, Adduct, PrecursorMz, and CCS columns.")
    parser.add_argument("output_db", help="Output .imsdb file name.")
    
    #csv_file = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/y9_posTransitionList copy 2.csv'
    csv_file = '/Users/danbru/Downloads/20250312_YeastMetabolites_CCSlibrary 1.csv'
    #output_db = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/ion_mobility_library/yeast9v2.imsdb'
    output_db = '/Users/danbru/Downloads/20250312_YeastMetabolites_CCSlibrary 1.imsdb'
    #args = parser.parse_args()
    #create_imsdb(csv_file, output_db)
    
    #csv_to_imsdb(csv_file, output_db)
    create_ion_mobility_db(csv_file, output_db, sep =';')
    




