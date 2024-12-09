#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:14:20 2024

@author: danbru
"""

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Hypothesis Generator with adjustable parameters.")
    
    # Define arguments with default values
    parser.add_argument('--folder', type=str, help='Working folder for experiment.')
    parser.add_argument('--volume', type=int, help='Working well volume in microliters.')
    parser.add_argument('--S', type=str, help='Stock concentration for supplement.')
    parser.add_argument('--Treatment', type=str, help='Stock concentration for treatment.')
    parser.add_argument('--environmental', type=int, default=0, help='1 of if environmental treatment, 0 otherwise.')

    # Parse arguments
    args = parser.parse_args()
    
    # Return as a dictionary for easy access and logging
    return vars(args)  

def main():
    '''
    
    Description of what this does? Not the prettiest thing ever...
    
    '''
    
    # Parse arguments
    args = parse_arguments()

    # Log each variable setting
    for arg, value in args.items():
        print(f'Variable "{arg}" set to: {value}')

    # Example usage of parsed values
    folder = args['folder']
    vol = int(args['volume'])
    s_stock_input = args['S']
    t_stock_input = args['Treatment']
    environmental = args['environmental']
    os.chdir(folder)
    
    """
    Calculate and return the concentrations and volumes to dispense into each well.
    """
    # Initialize Unit Registry
    ureg = UnitRegistry()
    Q_ = ureg.Quantity

    protocol_path = 'protocol/protocol.json'
    plate_layout_path = 'protocol/plate_layout/layout.tsv'
    
    # Load protocol and plate layout
    protocol = load_json(protocol_path)
    if protocol is None:
        raise ValueError("Protocol-JSON could not be loaded.")

    plate_layout = pd.read_csv(plate_layout_path, sep='\t', index_col=0)
    plate_layout = plate_layout.fillna('Media control')

    # Initialize dosing table with necessary columns
    dosing_table = pd.DataFrame(columns=['Summary', 'Media (uL)', 'Inoculated Media (uL)', 
                                         'Supplement (uL)', 'Treatment (uL)', 'MilliQ (uL)'])
    
    # Define fixed parameters
    V_total = Q_(vol, 'microliter')  # Total final volume per well
    D_pre_culture = 10  # Pre-culture dilution factor
    V_pre_culture = (V_total / D_pre_culture).to('microliter').magnitude  # Calculated inoculation volume
    
    # Extract supplement and treatment names from the first experiment
    first_experiment = next((exp for exp in protocol['experiments'] if exp.get('type', '').lower() == 'experiment'), None)
    if first_experiment:
        s_name = first_experiment.get('media_supplementation', 'Supplement')
        t_name = first_experiment.get('treatment', 'Treatment')
    else:
        s_name = 'Supplement'
        t_name = 'Treatment'
    
    # Parse stock concentrations
    s_stock_type, s_C_stock, s_stock_percentage_type = parse_concentration(s_stock_input, Q_)
    t_stock_type, t_C_stock, t_stock_percentage_type = parse_concentration(t_stock_input, Q_)
    
    # Adjust stock_percentage handling
    if s_stock_type == 'percentage':
        s_stock_percentage = s_C_stock.magnitude
    else:
        s_stock_percentage = None
    
    if t_stock_type == 'percentage':
        t_stock_percentage = t_C_stock.magnitude
    else:
        t_stock_percentage = None
    
    for experiment in protocol['experiments']:
        experiment_type = experiment.get('type', '').lower()
    
        # Initialize volumes
        s_V1 = 0.0  # Supplement volume
        t_V1 = 0.0  # Treatment volume
    
        # Extract and parse concentrations
        s_conc_str = experiment.get('media_supplementation_doses', '0')
        t_conc_str = experiment.get('treatment_parameters', '0')
    
        s_C_final_type, s_C_final, s_final_percentage_type = parse_concentration(s_conc_str, Q_)
        t_C_final_type, t_C_final, t_final_percentage_type = parse_concentration(t_conc_str, Q_)
    
        # Handle Supplement
        if s_C_final.magnitude > 0:
            if s_C_final_type == 'percentage':
                if s_stock_type != 'percentage':
                    raise ValueError("Stock concentration and final concentration units do not match for supplement.")
                s_V1 = calculate_percentage_volume(
                    s_C_final.magnitude,
                    V_total,
                    s_final_percentage_type,
                    s_stock_percentage,
                    Q_
                )
            else:
                # Convert to molar if needed
                compounds = pcp.get_compounds(s_name, 'name')
                if compounds:
                    molecular_weight = float(compounds[0].molecular_weight) * ureg('g/mol')
                    s_C_final_molar = convert_to_molar(s_C_final, molecular_weight)
                    s_C_stock_converted = convert_to_molar(s_C_stock, molecular_weight)
                    s_V1 = calculate_volume(s_C_stock_converted, s_C_final_molar, V_total)
                else:
                    print(f"Error: Compound '{s_name}' not found in PubChem.")
                    continue
        else:
            s_V1 = 0.0
    
        # Handle Treatment
        if t_C_final.magnitude > 0:
            if t_C_final_type == 'percentage':
                if t_stock_type != 'percentage':
                    raise ValueError("Stock concentration and final concentration units do not match for treatment.")
                t_V1 = calculate_percentage_volume(
                    t_C_final.magnitude,
                    V_total,
                    t_final_percentage_type,
                    t_stock_percentage,
                    Q_
                )
            else:
                # Convert to molar if needed
                compounds = pcp.get_compounds(t_name, 'name')
                if compounds:
                    molecular_weight = float(compounds[0].molecular_weight) * ureg('g/mol')
                    t_C_final_molar = convert_to_molar(t_C_final, molecular_weight)
                    t_C_stock_converted = convert_to_molar(t_C_stock, molecular_weight)
                    t_V1 = calculate_volume(t_C_stock_converted, t_C_final_molar, V_total)
                else:
                    print(f"Error: Compound '{t_name}' not found in PubChem.")
                    continue
        else:
            t_V1 = 0.0
    
        # Calculate remaining volumes
        V_rest = V_total - (V_pre_culture * ureg('uL')) - (t_V1 * ureg('uL')) - (s_V1 * ureg('uL'))
        V_additions = (t_V1 * ureg('uL')) + (s_V1 * ureg('uL'))
    
        # Calculate volumes for media and MilliQ
        V_media = V_additions + (V_rest - V_additions) / 2
        V_milliQ = (V_rest - V_additions) / 2
    
        # Ensure volumes are valid
        assert V_milliQ > 0, 'Too low stock concentrations resulting in negative MilliQ volume.'
        assert round(V_pre_culture + V_additions.magnitude + V_media.magnitude + V_milliQ.magnitude) == V_total.magnitude, 'Incorrect total volume.'
    
        # Create a new row with calculated volumes
        new_row = {
            'Summary': experiment.get('summary', 'Unnamed Experiment'),
            'Media (uL)': V_media.magnitude,
            'Inoculated Media (uL)': V_pre_culture,
            'Supplement (uL)': s_V1,
            'Treatment (uL)': t_V1,
            'MilliQ (uL)': V_milliQ.magnitude
        }
    
        # Append the new row to the dosing table
        dosing_table = pd.concat([dosing_table, pd.DataFrame([new_row])], ignore_index=True)

    # Add Media Control (no pre-culture, supplement, or treatment)
    media_control_row = {
        'Summary': 'Media control',
        'Media (uL)': (V_total / 2).magnitude,  # Complete to total volume (µL)
        'Inoculated Media (uL)': 0.0,
        'Supplement (uL)': 0.0,
        'Treatment (uL)': 0.0,
        'MilliQ (uL)': (V_total / 2).magnitude  # Complete to total volume (µL)
    }
    dosing_table = pd.concat([dosing_table, pd.DataFrame([media_control_row])], ignore_index=True)

    # Merge dosing table with plate layout
    plate_layout = plate_layout.reset_index()  # Ensure 'Summary' is a column
    merged_df = pd.merge(dosing_table, plate_layout[['Summary', 'well']], on='Summary', how='left')
    merged_df['well'] = merged_df['well'].str.replace(r'(\D)0*(\d+)', r'\1\2', regex=True)
    
    # Step 2: Sort by letter and then by number
    merged_df['letter'] = merged_df['well'].str[0]  # Extract letter part
    merged_df['number'] = merged_df['well'].str[1:].astype(int)  # Extract number part as integer
    merged_df = merged_df.sort_values(by=['number', 'letter']).reset_index(drop=True)
    
    # Iterate over each column and create a new column
    for col in merged_df.columns[1:-3]:
        # Add a new column with the logic: 1 if value != 0, else 0
        merged_df[f'{col}_ch'] = merged_df[col].apply(lambda x: 1 if x != 0 else 0)
    
    # Create a new folder for channels
    base_path = 'protocol/hamilton' 
    main_folder = "channels"
    subfolders_structure = {
        'milliq': [], 
        'media': [],
        'supplement': [],
        'treatment': [],
        'inoculant': [],
    }

    create_folder_structure(base_path, main_folder, subfolders_structure)    
    # Step 3: Save each number group to a separate sheet in Excel
    with pd.ExcelWriter('protocol/hamilton/runlist.xlsx') as writer:
        for num in range(1, 13):
            # Filter rows for the current number
            group_df = merged_df[merged_df['number'] == num].drop(columns=['letter', 'number'])
            
            # Write to a sheet named after the number (e.g., "Group 1", "Group 2", ...)
            group_df.to_excel(writer, sheet_name=f'{num}', index=False)
            
            # Write channel-valie
            
            for i, ch in enumerate(['media','inoculant','supplement','treatment','milliq']):
                channels = group_df.iloc[:,7+i]
                save_series_to_file(channels, f"protocol/hamilton/channels/{ch}/{num}.txt")
                
            
    merged_df.to_excel('protocol/hamilton/pipetting_layout.xlsx', index = False)
    # Group by "Summary" in df1 and aggregate the "well" values into a list
    result = merged_df.groupby('Summary')['well'].apply(list).reset_index()
    
    # Merge the result back into df1 to keep the original structure
    dosing_with_wells = pd.merge(dosing_table, result, on='Summary', how='left')
    dosing_with_wells.iloc[:, [1, 2, 3]] = dosing_with_wells.iloc[:, [1, 2, 3]].apply(pd.to_numeric, errors='coerce')

    # Extract "treatment_parameters" values
    if environmental > 0:
        copy_excel_file('protocol/hamilton/runlist.xlsx', len({experiment["treatment_parameters"] for experiment in protocol["experiments"]}))
    

import os
os.chdir(os.path.dirname(os.path.abspath(__file__)))
import pubchempy as pcp
from pint import UnitRegistry
import pandas as pd
from utils.plate_utils import *

if __name__ == "__main__":
    main()

