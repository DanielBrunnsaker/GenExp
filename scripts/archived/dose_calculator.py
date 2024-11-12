#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 21 14:27:41 2024

@author: danbru
"""

### Now i need to estimate concentrations, given a protocol?
def load_json(file_path):
    
    import json
    
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

def calculate_V1(C1, C2, V2):
    V1 = ((C2 * V2) / C1).to('uL')
    return V1

def percentage_to_volume(final_concentration_percent, total_volume):
    return (final_concentration_percent / 100) * total_volume  # Volume of solute in liters


def return_hamilton_concentrations(protocol_path, plate_layout):
    
    import pubchempy as pcp
    from pint import UnitRegistry
    import pandas as pd
    from tabulate import tabulate
    
    #plate_layout = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/generated_outputs/20241021_1833_layoutTable_glutamine_25_0.5_1.0_0.2.tsv', sep = '\t', index_col = 0)
    
   # protocol = load_json('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/generated_outputs/20241021_1833_protocol_glutamine_25_0.5_1.0_0.2.json')
    #experiment = protocol['experiments'][5]
    protocol = load_json(protocol_path)
    
    dosing_table = pd.DataFrame(columns = ['Summary','Media (uL)','Treatment (uL)','Supplement (uL)'])
        
    
    experiment_entry = [entry for entry in  protocol['experiments'] if 'type' in entry and 'experiment' in entry['type'].lower()][0]
    ureg = UnitRegistry()

    V2 = 0.25 * ureg('mL')  # Final volume

    s_name = experiment_entry['media_supplementation']
    t_name = experiment_entry['treatment']

    print(f'Enter stock concentration of {s_name} in mM')
    s_C1 = input()*ureg('mM')
    
    if'%' in experiment_entry['treatment_parameters']: 
        print(f'Enter stock concentration of {t_name} in %')
        t_C1 = input()*ureg('percent')
    else:
        print(f'Enter stock concentration of {t_name} in mM')
        t_C1 = input()*ureg('mM')

    for experiment in protocol['experiments']:
        
        if experiment['type'] == 'control':
            t_C2 = 0*ureg('uM')
            t_V1 = 0*ureg('uL')
            
            s_C2 = 0*ureg('uM')
            s_V1 = 0*ureg('uL')
            
        if experiment['type'] == 'supplementation control':
            t_C2 = 0*ureg('uM')
            t_V1 = 0*ureg('uL')
            
            s_dose = experiment['media_supplementation_doses'].split(' ')[0]
            s_dose_unit = experiment['media_supplementation_doses'].split(' ')[1]
            
            s_C2 = float(s_dose) * ureg(s_dose_unit)  # Final concentration
            s_C2_converted = convert_to_molar(s_C2, float(pcp.get_compounds(s_name, 'name')[0].molecular_weight)* ureg('g/mol'))
            #s_V1 = calculate_V1(s_C1, s_C2_converted, V2)
            temp = (s_C2_converted.to('mM') * V2.to('mL'))
            s_V1 = float(temp.m) / float(s_C1.m)*ureg('mL').to('uL')
            
            
        if experiment['type'] == 'treatment control':
            s_C2 = 0*ureg('uM')
            s_V1 = 0*ureg('uL')
            
            t_dose = experiment['treatment_parameters'].split(' ')[0]
            t_dose_unit = experiment['treatment_parameters'].split(' ')[1]
            
            if '%' in t_dose_unit:
                
                final_concentration_percent = float(t_dose)  # Desired final concentration as a percentage (e.g., 10%)
                solute_volume = percentage_to_volume(final_concentration_percent, V2)
            
                # Calculate volume of stock solution needed
                t_V1 = (solute_volume / (t_C1 / 100).m)*ureg('uL')
            
            else:
              
                t_C2 = float(t_dose) * ureg(t_dose_unit)  # Final concentration
                t_C2_converted = convert_to_molar(t_C2, float(pcp.get_compounds(t_name, 'name')[0].molecular_weight)* ureg('g/mol'))
                #t_V1 = calculate_V1(t_C1, t_C2_converted, V2)
                temp = (t_C2_converted * V2.to('mL'))
                t_V1 = float(temp.m) / float(t_C1.m)*ureg('mL').to('uL')
            
        if experiment['type'] == 'experiment':
            
            s_dose = experiment['media_supplementation_doses'].split(' ')[0]
            s_dose_unit = experiment['media_supplementation_doses'].split(' ')[1]
            
            s_C2 = float(s_dose) * ureg(s_dose_unit)  # Final concentration
            s_C2_converted = convert_to_molar(s_C2, float(pcp.get_compounds(s_name, 'name')[0].molecular_weight)* ureg('g/mol'))
            
            temp = (s_C2_converted.to('mM') * V2.to('mL'))
            s_V1 = float(temp.m) / float(s_C1.m)*ureg('mL').to('uL')
            #s_V1 = calculate_V1(s_C1, s_C2_converted, V2)
            
            t_dose = experiment['treatment_parameters'].split(' ')[0]
            t_dose_unit = experiment['treatment_parameters'].split(' ')[1]
            
            if '%' in t_dose_unit:
                
                final_concentration_percent = float(t_dose)  # Desired final concentration as a percentage (e.g., 10%)
                solute_volume = percentage_to_volume(final_concentration_percent, V2)
            
                # Calculate volume of stock solution needed
                t_V1 = (solute_volume / (t_C1 / 100).m)*ureg('uL')
            
            else:
              
                t_C2 = float(t_dose) * ureg(t_dose_unit)  # Final concentration
                t_C2_converted = convert_to_molar(t_C2, float(pcp.get_compounds(t_name, 'name')[0].molecular_weight)* ureg('g/mol'))
                #t_V1 = calculate_V1(t_C1, t_C2_converted, V2)
                temp = (t_C2_converted * V2.to('mL'))
                t_V1 = float(temp.m) / float(t_C1.m)*ureg('mL').to('uL')
            
            
            
            
        # Unit handling
        new_row = {'Summary': experiment['summary'],
                   'Media (uL)': (250*0.5-t_V1.m-s_V1.m),
                   'Inoculated Media (uL)': 250*0.5,
                   'Supplement (uL)': s_V1.m,
                   'Treatment (uL)': t_V1.m}
        # Add row in table?
        #dosing_table = dosing_table.append(new_row, ignore_index=True)
        dosing_table = pd.concat([dosing_table, pd.DataFrame([new_row])], ignore_index=True)

    dosing_table = pd.concat([dosing_table, 
                              pd.DataFrame([{'Summary': 'Media control',
                                           'Media (uL)': float(V2.to('uL').m),
                                           'Inoculated Media (uL)': 0,
                                           'Supplement (uL)': 0,
                                           'Treatment (uL)': 0}])], ignore_index=True)
    
    plate_layout['Summary'] = plate_layout['Summary'].fillna('Media control')
    #dosing_table.merge(plate_layout, left_on = 'Summary', right_on = 'Summary')        
    
    
    # Perform a left merge on the "Summary" column
    merged_df = pd.merge(dosing_table, plate_layout[['Summary', 'well']], on='Summary', how='left')
    
    # Group by "Summary" in df1 and aggregate the "well" values into a list
    result = merged_df.groupby('Summary')['well'].apply(list).reset_index()
    
    # Merge the result back into df1 to keep the original structure
    dosing_with_wells = pd.merge(dosing_table, result, on='Summary', how='left')
    dosing_with_wells.iloc[:, [1, 2, 3]] = dosing_with_wells.iloc[:, [1, 2, 3]].apply(pd.to_numeric, errors='coerce')
    
    print(tabulate(dosing_with_wells[['Media (uL)','Inoculated Media (uL)','Treatment (uL)','Supplement (uL)','well']], headers='keys', tablefmt='psql'))
        
