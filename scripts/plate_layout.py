#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:59:51 2024

@author: danbru
"""


# Function to load a JSON file
def load_json(file_path):
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return None

def save_json(data, file_path):
    try:
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)  # indent=4 is used for pretty-printing
        #print(f"Data successfully saved to {file_path}")
    except Exception as e:
        print(f"Error saving JSON file: {e}")

def run_minizinc_command(plate_file):
    # Use absolute paths for minizinc and files
    minizinc_path = "/Applications/MiniZincIDE.app/Contents/Resources/minizinc"
    mzn_file = os.environ["GEN_EXP_ROOT_DIR"] + "/plaid/plate-design.mzn"
    json_file = plate_file
    
    command = f"{minizinc_path} --solver Gecode {mzn_file} {json_file}"
    
    try:
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        if result.returncode == 0:
            return result.stdout
        else:
            print(f"Error: {result.stderr}")
            return None
    except Exception as e:
        print(f"Exception: {e}")
        return None

def output_to_dataframe(output_str):
    # Clean up the output to remove extra lines if needed
    output_str_cleaned = output_str.strip()  # Remove leading/trailing whitespaces

    # Use io.StringIO to simulate reading a CSV from a string
    data = io.StringIO(output_str_cleaned)

    # Read the string into a Pandas DataFrame
    df = pd.read_csv(data).iloc[:-1,:]
    
    return df

def visualize_grid(plate_layout):

    import numpy as np
    import pandas as pd
    from tabulate import tabulate
    
    # Create an empty grid (8 rows for A-H and 12 columns for 01-12)
    grid = pd.DataFrame(np.nan, index=list("ABCDEFGH"), columns=[f"{i:02d}" for i in range(1, 13)])

    # Map each unique value of CONCuM to a unique combination of S and T, handle "Media" as "Blank"
    unique_values = {val: ("Blank" if val == "Media" else f"S{idx // 2}T{idx % 2}") for idx, val in enumerate(plate_layout['CONCuM'].unique())}

    print(unique_values)  # Optional: Print unique mappings for reference

    # Populate the grid with the CONCuM values
    for _, row in plate_layout.iterrows():
        well = row['well']
        concum_value = row['CONCuM']
        
        # Get the value for the CONCuM (either "Blank" or SxTy)
        value = unique_values[concum_value]
        
        # Parse the well string (e.g., "A01" -> row "A", column "01")
        row_label = well[0]  # The first character (A-H)
        col_label = well[1:]  # The remaining part (01-12) 
        
        # Assign the value (SxTy or Blank) to the correct position in the grid
        grid.at[row_label, col_label] = value
    
    # Display the final grid
    print(tabulate(grid, headers='keys', tablefmt='psql'))
    #print(grid)


def plate_filler(plate_layout):

    # Define the list of all well positions (A01 to H12)
    all_wells = [f'{row}{col}' for row, col in itertools.product(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], [f'{i:02d}' for i in range(1, 13)])]
    # Find missing wells
    existing_wells = set(plate_layout['well'])
    missing_wells = set(all_wells) - existing_wells
    
    # Create placeholder rows for missing wells
    placeholder_rows = pd.DataFrame({'plateID': 'plate_1',
                                     'well': list(missing_wells),
                                     'cmpdname': 'Media',
                                     'CONCuM': 'Media control',
                                     'cmpdnum': 'Media',
                                     'VOLuL': 'nan'})
    
    # Append placeholder rows to the original dataframe
    df_filled = pd.concat([plate_layout, placeholder_rows], ignore_index=True)
    
    # Sort the DataFrame based on the well order (A01 to H12)
    df_filled['well'] = pd.Categorical(df_filled['well'], categories=all_wells, ordered=True)
    df_filled = df_filled.sort_values('well').reset_index(drop=True)

    return df_filled

#def generate_layout(json_data, reference_layout, output_path):
def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Hypothesis Generator with adjustable parameters.")
    
    # Define arguments with default values
    parser.add_argument('--folder', type=str, help='Working folder for experiment.')

    # Parse arguments
    args = parser.parse_args()
    
    # Return as a dictionary for easy access and logging
    return vars(args)

def main():
    
    # Parse arguments
    args = parse_arguments()

    # Log each variable setting
    for arg, value in args.items():
        print(f'Variable "{arg}" set to: {value}')

    # Example usage of parsed values
    folder = args['folder']
    os.chdir(folder)
    
    json_data = load_json('protocol/protocol.json')
    reference_layout = load_json('../../plaid/reference_plate.json')
    output_path = 'protocol/plate_layout/minizinc_reference.json'
    
    experiment_entry = [entry for entry in  json_data['experiments'] if 'type' in entry and 'experiment' in entry['type'].lower()][0]
    treatment_name = experiment_entry['treatment']
    supplementation_name = experiment_entry['media_supplementation']
    
    combinations = []
    
    df = pd.DataFrame(columns = ['Summary','Experiment Parameters'])
    
    for exp in json_data['experiments']:
        
     
        
        supplementation_dose = exp['media_supplementation_doses']
        treatment_dose = exp['treatment_parameters']
        
        if supplementation_dose == None:
            supplementation_dose = '0 mM'
        
        if treatment_dose == None:
            treatment_dose = '0 mM'
        
        supp = f'{supplementation_name}: '+supplementation_dose
        treatment = f'{treatment_name}: '+treatment_dose
                
        combinations.append(supp+' & '+treatment)
        
        new_row = {'Summary':exp['summary'], 
                   'Experiment parameters':supp+' & '+treatment}
        
        df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

        
    reference_layout['compound_concentrations'] = [len(combinations)]
    reference_layout['compound_concentration_names'] = [combinations]
    
    for repl in range(14,4,-1):
        
        reference_layout['compound_replicates'] = [repl]
        save_json(reference_layout, output_path)
        output = run_minizinc_command(output_path)
        
        if output == '=====UNSATISFIABLE=====\n':
            continue
        else:
            plate_layout = output_to_dataframe(output)
            print(plate_layout.shape[0])
            
            if 96-plate_layout.shape[0] < repl:
                continue
            
            print('Layout completed.')
            break
    
    plate_layout = plate_filler(plate_layout)
    plate_layout = plate_layout.merge(df, 
                                      left_on = 'CONCuM', 
                                      right_on = 'Experiment parameters', 
                                      how = 'left')[['plateID','well','CONCuM','Summary']]
    #visualize_grid(plate_layout)
    
    plate_save_path = "protocol/plate_layout/layout.tsv"
    plate_layout.to_csv(plate_save_path, sep = '\t')
    

import os
import itertools
import json
import subprocess
import pandas as pd
import io

if __name__ == "__main__":
    main()
