#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 11:14:20 2024

@author: danbru
"""
from pathlib import Path
import pandas as pd

from pint import UnitRegistry
import pubchempy as pcp

# util imports
from utils.hamilton_utils import *
from utils.plate_utils import *
from utils.hgen_support import *
import config

def hamilton_protocol(output_folder):

    """
    Calculate and return the concentrations and volumes to dispense into each well.
    Logic seems robust, but could be prettier. Someone please fix?
    """

    EXPERIMENT_DIR = Path(output_folder)
    protocol_path = EXPERIMENT_DIR / 'protocol/protocol.json'
    plate_layout_path = EXPERIMENT_DIR / 'protocol/plate_layout/layout.tsv'

    # Load the layout
    plate_layout = pd.read_csv(plate_layout_path, sep='\t', index_col=0)
    plate_layout = plate_layout.fillna('Media control')

    # Initialize unit registry
    ureg = UnitRegistry()
    Q_ = ureg.Quantity

    # Load the protocol
    protocol = load_json(protocol_path)
    if protocol is None:
        raise ValueError("Protocol-JSON could not be loaded.")

    # Extract the names and doses of involved compounds
    compound_dict = extract_highest_doses(protocol, ureg, Q_)
    supplement_name = compound_dict['supplement'][0]
    negative_control_name = compound_dict['negative_control'][0]
    treatment_name = compound_dict['treatment'][0]
    treatment_name  = treatment_name.replace('Aluminium','Aluminum') # Stupid languages

    ## Not he nicest block ever, but works
    try:
        supplement_concentrations = find_rows_by_inchikey(EXPERIMENT_DIR / '../../data/compound_library/library.xlsx',
                                                          pcp.get_compounds(supplement_name, 'name')[0].inchikey)
        supplement_stock_concentration = get_closest_concentration(supplement_concentrations,
                                                                   compound_dict['supplement'][1], ureg)['Available concentrations']
        print('  -Rescue-agent found in compound library!')
    except:
        supplement_stock_concentration = input(f"  -Please provide {supplement_name} stock concentration (recommended at least {5*compound_dict['supplement'][1]}): ")

    try:
        negative_concentrations = find_rows_by_inchikey(EXPERIMENT_DIR / '../../data/compound_library/library.xlsx',
                                                        pcp.get_compounds(negative_control_name, 'name')[0].inchikey)
        negative_stock_concentration = get_closest_concentration(negative_concentrations,
                                                                 compound_dict['negative_control'][1], ureg)['Available concentrations']
        print('  -Negative control found in compound library!')
    except:
        negative_stock_concentration = input(f"  -Please provide {negative_control_name} stock concentration (recommended at least {5*compound_dict['negative_control'][1]}): ")

    try:
        treatment_concentrations = find_rows_by_inchikey(EXPERIMENT_DIR / '../../data/compound_library/library.xlsx',
                                                         pcp.get_compounds(treatment_name.replace(' derivative',''), 'name')[0].inchikey)
        treatment_stock_concentration = get_closest_concentration(treatment_concentrations,
                                                                  compound_dict['treatment'][1], ureg)['Available concentrations']
        print('  -Treatment compound found in compound library! \n')
    except:
        treatment_stock_concentration = input(f"  -Please provide {treatment_name} stock concentration (recommended at least {5*compound_dict['treatment'][1]}): ")


    # Initialize dosing table with necessary columns
    dosing_table = pd.DataFrame(columns=['Summary', 'Media (uL)', 'Inoculated Media (uL)',
                                         'Supplement (uL)', 'Negative Control (uL)','Treatment (uL)', 'MilliQ (uL)'])

    # Define fixed parameters
    V_total = Q_(config.WELL_VOLUME, config.WELL_VOLUME_UNIT)  # Total final volume per well
    D_pre_culture = config.DILUTION_FACTOR  # Pre-culture dilution factor
    V_pre_culture = (V_total / D_pre_culture).to('microliter').magnitude  # Calculated inoculation volume

    # Parse stock concentrations
    s_stock_type, s_C_stock, s_stock_percentage_type = parse_concentration(supplement_stock_concentration, Q_)
    t_stock_type, t_C_stock, t_stock_percentage_type = parse_concentration(treatment_stock_concentration, Q_)
    n_stock_type, n_C_stock, n_stock_percentage_type = parse_concentration(negative_stock_concentration, Q_)

    # Adjust stock_percentage handling
    if s_stock_type == 'percentage':
        s_stock_percentage = s_C_stock.magnitude
    else:
        s_stock_percentage = None

    if t_stock_type == 'percentage':
        t_stock_percentage = t_C_stock.magnitude
    else:
        t_stock_percentage = None

    if n_stock_type == 'percentage':
        n_stock_percentage = t_C_stock.magnitude
    else:
        n_stock_percentage = None

    for experiment in protocol['experiments']:
        experiment_type = experiment.get('type', '').lower()

        # Initialize volumes
        s_V1 = 0.0  # Supplement volume
        t_V1 = 0.0  # Treatment volume
        n_V1 = 0.0  # Negative control volume

        # Extract and parse concentrations
        s_conc_str = experiment.get('media_supplementation_doses', '0')
        t_conc_str = experiment.get('treatment_parameters', '0')

        if 'negative' in experiment_type:
            n_conc_str = experiment.get('media_supplementation_doses', '0')
            n_C_final_type, n_C_final, n_final_percentage_type = parse_concentration(n_conc_str, Q_)

        s_C_final_type, s_C_final, s_final_percentage_type = parse_concentration(s_conc_str, Q_)
        t_C_final_type, t_C_final, t_final_percentage_type = parse_concentration(t_conc_str, Q_)

        # Handle negative control
        if 'negative' in experiment_type:
            n_V1 = compute_volume_for_compound(
                negative_control_name,
                n_C_final,
                n_C_final_type,
                n_final_percentage_type,
                n_C_stock,
                n_stock_type,   # If you track negative stock type
                n_stock_percentage,  # If applicable
                V_total,
                Q_,
                ureg,
                print_name_for_error='negative control'
            )
        else:
            # Handle supplement
            s_V1 = compute_volume_for_compound(
                supplement_name,
                s_C_final,
                s_C_final_type,
                s_final_percentage_type,
                s_C_stock,
                s_stock_type,
                s_stock_percentage,
                V_total,
                Q_,
                ureg,
                print_name_for_error='supplement'
            )

        # Handle treatment
        t_V1 = compute_volume_for_compound(
            treatment_name.replace(' derivative',''),
            t_C_final,
            t_C_final_type,
            t_final_percentage_type,
            t_C_stock,
            t_stock_type,
            t_stock_percentage,
            V_total,
            Q_,
            ureg,
            print_name_for_error='treatment'
        )

        # Calculate remaining volumes
        V_rest = V_total - (
            V_pre_culture * ureg('uL')
            + t_V1 * ureg('uL')
            + s_V1 * ureg('uL')
            + n_V1 * ureg('uL')
        )
        V_additions = (t_V1 + s_V1 + n_V1) * ureg('uL')
        V_media = V_additions + (V_rest - V_additions) / 2
        V_milliQ = (V_rest - V_additions) / 2

        # Assertions, checks...
        assert V_milliQ > 0, 'Too low stock concentrations...'
        assert round(
            V_pre_culture
            + V_additions.magnitude
            + V_media.magnitude
            + V_milliQ.magnitude
        ) == V_total.magnitude, 'Incorrect total volume...'

        # Build new row
        new_row = {
            'Summary': experiment.get('summary', 'Unnamed Experiment'),
            'Media (uL)': V_media.magnitude,
            'Inoculated Media (uL)': V_pre_culture,
            'Supplement (uL)': s_V1,
            'Negative Control (uL)': n_V1,
            'Treatment (uL)': t_V1,
            'MilliQ (uL)': V_milliQ.magnitude
        }

        # Append to table
        dosing_table = pd.concat([dosing_table, pd.DataFrame([new_row])], ignore_index=True)

    # Add Media Control (no pre-culture, supplement, or treatment)
    media_control_row = {
        'Summary': 'Media control',
        'Media (uL)': (V_total / 2).magnitude,  # Complete to total volume (µL)
        'Inoculated Media (uL)': 0.0,
        'Supplement (uL)': 0.0,
        'Negative Control (uL)': 0.0,
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

    create_folder_structure(EXPERIMENT_DIR / 'protocol/hamilton', "channels", {
        'milliq': [], 
        'media': [],
        'supplement': [],
        'negative': [],
        'treatment': [],
        'inoculant': [],
    })

    # Step 3: Save each number group to a separate sheet in Excel
    with pd.ExcelWriter(EXPERIMENT_DIR / 'protocol/hamilton/runlist.xlsx') as writer:
        for num in range(1, 13):
            # Filter rows for the current number
            group_df = merged_df[merged_df['number'] == num].drop(columns=['letter', 'number'])

            # Write to a sheet named after the number (e.g., "Group 1", "Group 2", ...)
            group_df.to_excel(writer, sheet_name=f'{num}', index=False)

            # Write channel-value
            for i, ch in enumerate(['media','inoculant','supplement','negative','treatment','milliq']):
                channels = group_df.iloc[:,8+i]
                save_series_to_file(channels, EXPERIMENT_DIR / f"protocol/hamilton/channels/{ch}/{num}.txt")

    merged_df.to_excel(EXPERIMENT_DIR / 'protocol/hamilton/pipetting_layout.xlsx', index = False)

    # Group by "Summary" in df1 and aggregate the "well" values into a list
    result = merged_df.groupby('Summary')['well'].apply(list).reset_index()

    # Merge the result back into df1 to keep the original structure
    dosing_with_wells = pd.merge(dosing_table, result, on='Summary', how='left')
    dosing_with_wells.iloc[:, [1, 2, 3]] = dosing_with_wells.iloc[:, [1, 2, 3]].apply(pd.to_numeric, errors='coerce')

if __name__ == "__main__":

    import argparse
    parser = argparse.ArgumentParser(description="Hamilton protocol generator")
    parser.add_argument("--output_folder", required=True, type=str, help="Experiment folder")

    args = parser.parse_args()
    hamilton_protocol(args.output_folder)
