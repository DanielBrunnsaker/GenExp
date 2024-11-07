#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  7 13:14:52 2024

@author: danbru
"""

import json
import pubchempy as pcp
from pint import UnitRegistry
import pandas as pd
from tabulate import tabulate

def load_json(file_path):
    """
    Load a JSON file and return its content.
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return None

def convert_to_molar(concentration, molecular_weight):
    if concentration.check('[mass] / [volume]'):  # Check if concentration is mass/volume
        return (concentration / molecular_weight).to('mM')  # Convert mass-based to molar
    return concentration.to('mM')  # Already in molar units

def calculate_volume(C_stock, C_final, V_total, ureg):
    """
    Calculate the volume to dispense from stock to achieve desired final concentration.
    """
    V1 = (C_final * V_total / C_stock).to('microliter')
    return V1.magnitude  # Return as float

def percentage_to_volume(final_concentration_percent, total_volume, ureg):
    """
    Convert percentage concentration to volume of solute.
    """
    return (final_concentration_percent / 100) * total_volume.to('microliter').magnitude  # Return as float

def return_hamilton_concentrations(protocol_path, plate_layout_path):
    """
    Calculate and return the concentrations and volumes to dispense into each well.
    """
    # Initialize Unit Registry
    ureg = UnitRegistry()
    Q_ = ureg.Quantity


    #protocol_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/generated_outputs/20241107_1440_protocol_alanine_10_0.75_1.0_0.1.json'
    #plate_layout_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/generated_outputs/20241107_1440_layoutTable_alanine_10_0.75_1.0_0.1.tsv'
    
    # Load protocol and plate layout
    protocol = load_json(protocol_path)
    if protocol is None:
        return

    plate_layout = pd.read_csv(plate_layout_path, sep='\t', index_col=0)
    plate_layout = plate_layout.fillna('Media control')

    # Initialize dosing table with necessary columns
    dosing_table = pd.DataFrame(columns=['Summary', 'Media (uL)', 'Inoculated Media (uL)', 
                                         'Supplement (uL)', 'Treatment (uL)', 'MilliQ (uL)'])

    # Define fixed parameters
    V_total = Q_(225, 'microliter')  # Total final volume per well
    D_pre_culture = 20  # Pre-culture dilution factor
    V_pre_culture = (V_total / D_pre_culture).to('microliter').magnitude  # Calculated inoculation volume
    
    # Extract supplement and treatment names from the first experiment
    first_experiment = next((exp for exp in protocol['experiments'] if exp.get('type', '').lower() == 'experiment'), None)
    if first_experiment:
        s_name = first_experiment.get('media_supplementation', 'Supplement')
        t_name = first_experiment.get('treatment', 'Treatment')
    else:
        s_name = 'Supplement'
        t_name = 'Treatment'
    
    # Prompt user for stock concentrations
    print(f"Enter stock concentration for {s_name} (e.g., '10 mM' or '5 mg/mL'):")
    s_stock_input = input().strip()
    s_C_stock = Q_(s_stock_input)

    print(f"Enter stock concentration for {t_name} (e.g., '5 mg/mL' or '10%'):")
    t_stock_input = input().strip()
    t_C_stock = Q_(t_stock_input)

    for experiment in protocol['experiments']:
        experiment_type = experiment.get('type', '').lower()

        # Initialize concentrations and volumes
        s_V1 = 0.0  # Supplement volume
        t_V1 = 0.0  # Treatment volume

        # Handle different experiment types
        if experiment_type == 'control':
            # No supplement or treatment
            s_C_final = Q_(0, 'mM')
            t_C_final = Q_(0, 'mM')
        elif experiment_type == 'supplementation control':
            
            
            if '%' in experiment['media_supplementation_doses']:
                s_dose = experiment['media_supplementation_doses'][:-1]
                s_dose_unit = experiment['media_supplementation_doses'][-1]
            else:
                s_dose = experiment['media_supplementation_doses'].split(' ')[0]
                s_dose_unit = experiment['media_supplementation_doses'].split(' ')[1]
                
            # Only supplement is added
            s_C_final = Q_(float(s_dose), s_dose_unit)
            t_C_final = Q_(0, 'mM')
        elif experiment_type == 'treatment control':
            
            #t_dose = experiment['treatment_parameters'].split(' ')[0]
            #t_dose_unit = experiment['treatment_parameters'].split(' ')[1]
            
            if '%' in experiment['treatment_parameters']:
                t_dose = experiment['treatment_parameters'][:-1]
                t_dose_unit = experiment['treatment_parameters'][-1]
            else:
                t_dose = experiment['treatment_parameters'].split(' ')[0]
                t_dose_unit = experiment['treatment_parameters'].split(' ')[1]
            
            # Only treatment is added
            s_C_final = Q_(0, 'mM')
            t_C_final = Q_(float(t_dose), t_dose_unit)
        elif experiment_type == 'experiment':
            
            
            # S
            if '%' in experiment['media_supplementation_doses']:
                s_dose = experiment['media_supplementation_doses'][:-1]
                s_dose_unit = experiment['media_supplementation_doses'][-1]
            else:
                s_dose = experiment['media_supplementation_doses'].split(' ')[0]
                s_dose_unit = experiment['media_supplementation_doses'].split(' ')[1]
            
            
            # T
            if '%' in experiment['treatment_parameters']:
                t_dose = experiment['treatment_parameters'][:-1]
                t_dose_unit = experiment['treatment_parameters'][-1]
            else:
                t_dose = experiment['treatment_parameters'].split(' ')[0]
                t_dose_unit = experiment['treatment_parameters'].split(' ')[1]
            
            
            # Both supplement and treatment are added
            s_C_final = Q_(float(s_dose), s_dose_unit)
            t_C_final = Q_(float(t_dose), t_dose_unit)
        else:
            # Undefined experiment type; skip or handle accordingly
            continue

        # Convert supplement final concentration to molar if necessary
        if s_C_final.magnitude > 0:
            # Get molecular weight from PubChem
            compounds = pcp.get_compounds(s_name, 'name')
            if compounds:
                molecular_weight = float(compounds[0].molecular_weight) * ureg('g/mol')
                s_C_final = convert_to_molar(s_C_final, molecular_weight)
                s_C_stock_converted = convert_to_molar(s_C_stock, molecular_weight)
            else:
                print(f"Error: Compound '{s_name}' not found in PubChem.")
                continue

            # Calculate Supplement Volume (s_V1)
            s_V1 = calculate_volume(s_C_stock_converted, s_C_final, V_total, ureg)
        else:
            s_V1 = 0.0

        # Convert treatment final concentration to molar if necessary
        if t_C_final.magnitude > 0:
            if '%' in t_dose_unit:
                # Handle percentage-based treatment
                final_conc_percent = float(t_dose)
                t_V1 = percentage_to_volume(final_conc_percent, V_total, ureg)
            else:
                # Get molecular weight from PubChem
                compounds = pcp.get_compounds(t_name, 'name')
                if compounds:
                    molecular_weight = float(compounds[0].molecular_weight) * ureg('g/mol')
                    t_C_final = convert_to_molar(t_C_final, molecular_weight)
                    t_C_stock_converted = convert_to_molar(t_C_stock, molecular_weight)
                else:
                    print(f"Error: Compound '{t_name}' not found in PubChem.")
                    continue

                # Calculate Treatment Volume (t_V1)
                t_V1 = calculate_volume(t_C_stock_converted, t_C_final, V_total, ureg)
        else:
            t_V1 = 0.0

        V_rest = V_total - (V_pre_culture * ureg('uL'))-(t_V1*ureg('uL'))-(s_V1*ureg('uL'))
        V_additions = (t_V1*ureg('uL'))+(s_V1*ureg('uL'))

        # How much 2x media to add 
        V_media = V_additions+(V_rest-V_additions)/2
        V_milliQ = (V_rest-V_additions)/2

        # Check so that volumes match
        V_pre_culture*ureg('uL') + V_additions + V_media + V_milliQ
        
        # Simple asserts to ensure valid volumes
        # E.g. Total volume adds up and that no volumes are negative (due to too low stock concentrations)
        assert V_milliQ > 0, 'Too low stock concentrations'
        assert V_pre_culture*ureg('uL') + V_additions + V_media + V_milliQ == V_total, 'Incorrect total volume.'

        # Create a new row with calculated volumes
        new_row = {
            'Summary': experiment.get('summary', 'Unnamed Experiment'),
            'Media (uL)': V_media.m,
            'Inoculated Media (uL)': V_pre_culture,
            'Supplement (uL)': s_V1,
            'Treatment (uL)': t_V1,
            'MilliQ (uL)': V_milliQ.m
        }

        # Append the new row to the dosing table
        dosing_table = pd.concat([dosing_table, pd.DataFrame([new_row])], ignore_index=True)

    
    # Add Media Control (no pre-culture, supplement, or treatment)
    media_control_row = {
        'Summary': 'Media control',
        'Media (uL)': (V_total/2).m, # Complete to total volume (µL)
        'Inoculated Media (uL)': 0.0,
        'Supplement (uL)': 0.0,
        'Treatment (uL)': 0.0,
        'MilliQ (uL)': (V_total/2).m  # Complete to total volume (µL)
    }
    dosing_table = pd.concat([dosing_table, pd.DataFrame([media_control_row])], ignore_index=True)


    # Merge dosing table with plate layout
    plate_layout = plate_layout.reset_index()  # Ensure 'Summary' is a column
    merged_df = pd.merge(dosing_table, plate_layout[['Summary', 'well']], on='Summary', how='left')
    
    # Group by "Summary" in df1 and aggregate the "well" values into a list
    result = merged_df.groupby('Summary')['well'].apply(list).reset_index()
    
    # Merge the result back into df1 to keep the original structure
    dosing_with_wells = pd.merge(dosing_table, result, on='Summary', how='left')
    dosing_with_wells.iloc[:, [1, 2, 3]] = dosing_with_wells.iloc[:, [1, 2, 3]].apply(pd.to_numeric, errors='coerce')
    
    
    #['Summary', 'Media (uL)', 'Inoculated Media (uL)', 'Supplement (uL)', 'Treatment (uL)', 'MilliQ (uL)', 'well']
    print(tabulate(dosing_with_wells[['Media (uL)', 'Inoculated Media (uL)', 'Supplement (uL)', 'Treatment (uL)', 'MilliQ (uL)', 'well']], headers='keys', tablefmt='psql'))
