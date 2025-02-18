#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb  4 10:04:58 2025

@author: danbru
"""

import pandas as pd
import random
from pathlib import Path
import config

def generate_runlist(wells_df, samples_in_sequence=3, tune_interval=4, randomize=True):
    """
    Generate a runlist with Wash, Blank, Random Samples, QC samples, and Tune injections.

    Parameters:
        wells_df (pd.DataFrame): DataFrame containing well information (columns: "Well", "Notes").
        output_file (str): Path to save the generated runlist (CSV format).
        samples_in_sequence (int): Number of random samples between Wash/Blank and QC.
        tune_interval (int): Number of sequences after which a Tune injection is added.
        randomize (bool): Whether to randomize the sample order or maintain the original order.
    """
    
    def structured_randomization(wells_df, samples_in_sequence):
        """Ensures each sequence contains at least one of each experimental group."""
        if not randomize:
            return wells_df.reset_index(drop=True)
        
        grouped = [group_df.sample(frac=1, random_state=42).reset_index(drop=True)
                   for _, group_df in wells_df.groupby("Summary")]
        samples = []
        while any(len(g) > 0 for g in grouped):
            batch = [g.iloc[0] for g in grouped if not g.empty]
            samples.extend(batch)
            for i, g in enumerate(grouped):
                if not g.empty:
                    grouped[i] = g.iloc[1:].reset_index(drop=True)
            random.shuffle(samples)
        
        return pd.DataFrame(samples)
    
    # Apply structured randomization or maintain original order
    randomized_samples = structured_randomization(wells_df[['well','Summary']], samples_in_sequence)
    
    # Initialize the runlist
    runlist = []
    
    # Pre run wash
    #for i in range(2):
    #    runlist.append({"Description": 'Wash', "Well": 'WASH1', "Notes": 'Acqeous Wash', "Sample_Type": 'WASH'})
    #    runlist.append({"Description": 'Wash', "Well": 'WASH2', "Notes": 'Organic Wash', "Sample_Type": 'WASH'})
    
    # Add a tune injection at the start
    runlist.append({"Description": "Tune", "Well": "MAT3", "Notes": 'Tuning mix injection', "Sample_Type": 'TUNE'})
    sequence = ["Blank", "Blank", "QC"] + ["Random Sample"] * samples_in_sequence + ["QC","Blank"]
    
    # Iterate over randomized samples and construct the runlist
    sample_idx = 0  # Track position in the randomized sample list
    sequence_count = 0  # Track number of completed sequences

    while sample_idx < len(randomized_samples):
        for step in sequence:
            if step == "Wash":
                runlist.append({"Description": step, "Well": 'WASH1', "Notes": 'Acqeous Wash', "Sample_Type": 'WASH'})
                runlist.append({"Description": step, "Well": 'WASH2', "Notes": 'Organic Wash', "Sample_Type": 'WASH'})
            elif step == "Blank":
                runlist.append({"Description": step, "Well": "MAT1", "Notes": 'Blank', "Sample_Type": 'BLANK'})
            elif step == "QC":
                runlist.append({"Description": step, "Well": "MAT2", "Notes": 'QC', "Sample_Type": 'SAMPLE'})
            elif step == "Random Sample":
                if sample_idx < len(randomized_samples):
                    runlist.append({
                        "Description": "Sample",
                        "Well": randomized_samples.iloc[sample_idx, 0],
                        "Notes": randomized_samples.iloc[sample_idx, 1],
                        'Sample_Type': 'SAMPLE'
                    })
                    sample_idx += 1

        sequence_count += 1

        # Add a Tune injection every tune_interval sequences
        if sequence_count % tune_interval == 0:
            #runlist.append({"Description": step, "Well": 'WASH1', "Notes": 'Acqeous Wash', "Sample_Type": 'WASH'})
            #runlist.append({"Description": step, "Well": 'WASH2', "Notes": 'Organic Wash', "Sample_Type": 'WASH'})
            runlist.append({"Description": step, "Well": "MAT1", "Notes": 'Blank', "Sample_Type": 'BLANK'})
            runlist.append({"Description": "Tune", "Well": "MAT3", "Notes": 'Tuning mix injection', "Sample_Type": 'TUNE'})
            runlist.append({"Description": step, "Well": "MAT1", "Notes": 'Blank', "Sample_Type": 'BLANK'})
    
    # Post run wash
    for i in range(2):
        runlist.append({"Description": step, "Well": 'WASH1', "Notes": 'Acqeous Wash', "Sample_Type": 'WASH'})
        runlist.append({"Description": step, "Well": 'WASH2', "Notes": 'Organic Wash', "Sample_Type": 'WASH'})
        
    return pd.DataFrame(runlist)

def generate_mass_spec_metadata(output_folder):
    
    EXPERIMENT_DIR = Path(output_folder)
    layout = pd.read_csv(EXPERIMENT_DIR / 'protocol/plate_layout/layout.tsv', sep = '\t', index_col = 0)
    
    # So, likely randomize the injection order?
    samples_in_sequence = config.SAMPLES_IN_SEQUENCE
    tune_interval = config.TUNE_INTERVAL
    
    # Generate the runlist
    runlist = generate_runlist(layout, samples_in_sequence, tune_interval, randomize = config.RUNLIST_RANDOMIZATION)
    
    # Reformat runlist into AutonoMS formatting
    runlist['Sequence'] = config.MS_POLARITY
    runlist['6560_Method'] = config.METHOD_NAME
    runlist['Column_Type'] = config.COLUMN_TYPE
    runlist['Plate_Type'] = config.PLATE_TYPE
    runlist['Sample_Number'] = 1
    runlist['Replicate_Number'] = 1
    
    runlist = runlist[['Well','Description','Sequence','Sample_Number','Replicate_Number','Sample_Type','6560_Method','Plate_Type','Column_Type','Notes']]
    runlist['Well'] = runlist['Well'].str.replace(r'(\D)0+(\d+)', r'\1\2', regex=True)
    
    # Now, for the rf_params sheet? Maybe have this fixed as is (e.g. load from excel sheet)
    # Same with sample sheet
    
    with pd.ExcelWriter(EXPERIMENT_DIR / 'protocol/mass_spectrometry/ExperimentTemplate.xlsx', engine='openpyxl') as writer:
        runlist.to_excel(writer, sheet_name='samples', index=False)
    
    rf_settings = pd.read_excel(config.RF_SETTINGS_FILE , sheet_name='rf_params')
    headers = list(rf_settings.columns)  # Get the first row as headers
    headers = [headers[i] if i < 2 else "" for i in range(len(headers))]  # Keep first two, empty the rest
    rf_settings.columns = headers  # Assign the modified headers
    rf_settings = rf_settings.iloc[1:]  # Remove the first row (since it's now set as headers)
    
    analysis_settings = pd.read_excel(config.RF_SETTINGS_FILE, sheet_name='data_analysis')
    
    # Append the existing sheet to the newly created Excel file
    with pd.ExcelWriter(EXPERIMENT_DIR / 'protocol/mass_spectrometry/ExperimentTemplate.xlsx', engine='openpyxl', mode='a') as writer:
        rf_settings.to_excel(writer, sheet_name='rf_params', index=False)
    
    with pd.ExcelWriter(EXPERIMENT_DIR / 'protocol/mass_spectrometry/ExperimentTemplate.xlsx', engine='openpyxl', mode='a') as writer:
        analysis_settings.to_excel(writer, sheet_name='data_analysis', index=False)
    
    

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Experiment folder")
    parser.add_argument("--output_folder", required=True, type=str, help="Experiment folder")

    args = parser.parse_args()
    generate_mass_spec_metadata(args.output_folder)
    
    
    # Example
    # python ms_runlist_generation.py --output_folder "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_202501281618"
    
    
    

