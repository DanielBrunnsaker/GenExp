#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 19 13:52:18 2024

@author: danbru
"""

import itertools
import json
import subprocess
import pandas as pd
import io

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
    mzn_file = "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/plaid/plate-design.mzn"
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
'''
def visualize_grid(plate_layout):
    import numpy as np
    import pandas as pd
    from tabulate import tabulate
    
    # Create an empty grid (8 rows for A-H and 12 columns for 01-12)
    grid = pd.DataFrame(np.nan, index=list("ABCDEFGH"), columns=[f"{i:02d}" for i in range(1, 13)])

    # Extract unique values for the first and second compound from the 'CONCuM' strings
    def parse_concums(concum):
        """Helper function to extract concentrations from the CONCuM string."""
        parts = concum.split(' & ')
        compound1 = parts[0].split(': ')[1]  # Get the concentration of the first compound (glutamine)
        compound2 = parts[1].split(': ')[1]  # Get the concentration of the second compound (nickel sulfate)
        return compound1, compound2

    # Get unique concentrations for the first and second compounds
    concum_parsed = plate_layout['CONCuM'].apply(parse_concums)
    unique_s_values = sorted(concum_parsed.apply(lambda x: x[0]).unique())
    unique_t_values = sorted(concum_parsed.apply(lambda x: x[1]).unique())

    # Create mappings for S (first compound) and T (second compound)
    s_mapping = {val: idx for idx, val in enumerate(unique_s_values)}
    t_mapping = {val: idx for idx, val in enumerate(unique_t_values)}

    print("S Mapping (First Compound):", s_mapping)  # Optional: Print S mapping
    print("T Mapping (Second Compound):", t_mapping)  # Optional: Print T mapping

    # Populate the grid with the corresponding SxTy values
    for _, row in plate_layout.iterrows():
        well = row['well']
        concum_value = row['CONCuM']
        
        # Parse the concentration values
        s_value, t_value = parse_concums(concum_value)
        
        # Get the S and T codes
        s_code = s_mapping[s_value]
        t_code = t_mapping[t_value]
        
        # Create the SxTy value (e.g., "S2T0")
        value = f"S{s_code}T{t_code}"
        
        # Parse the well string (e.g., "A01" -> row "A", column "01")
        row_label = well[0]  # The first character (A-H)
        col_label = well[1:]  # The remaining part (01-12) 
        
        # Assign the SxTy value to the correct position in the grid
        grid.at[row_label, col_label] = value
    
    # Display the final grid
    print(tabulate(grid, headers='keys', tablefmt='psql'))

'''
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




def generate_layout(json_data, reference_layout, output_path):

    #json_data = load_json('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/generated_outputs/20241021_1833_protocol_glutamine_25_0.5_1.0_0.2.json')
    #reference_layout = load_json('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/plaid/reference_plate.json')
    json_data = load_json(json_data)
    reference_layout = load_json(reference_layout)
    
    experiment_entry = [entry for entry in  json_data['experiments'] if 'type' in entry and 'experiment' in entry['type'].lower()][0]
    treatment_name = experiment_entry['treatment']
    supplementation_name = experiment_entry['media_supplementation']
    
    combinations = []
    
    df = pd.DataFrame(columns = ['Summary','Experiment Parameters'])
    
    for exp in json_data['experiments']:
        
        #supplementation_name = exp['media_supplementation']
        #treatment_name = exp['treatment']
        
        supplementation_dose = exp['media_supplementation_doses']
        treatment_dose = exp['treatment_parameters']
        
        if supplementation_dose == None:
            supplementation_dose = '0 mM'
        #    supplementation_name = 'control'
        
        if treatment_dose == None:
            treatment_dose = '0 mM'
        #    treatment_name = 'control'
        
        supp = f'{supplementation_name}: '+supplementation_dose
        treatment = f'{treatment_name}: '+treatment_dose
        
        #supp+'___'+treatment
        
        combinations.append(supp+' & '+treatment)
        
        
        new_row = {'Summary':exp['summary'], 
                   'Experiment parameters':supp+' & '+treatment}
        
        df = pd.concat([df, pd.DataFrame([new_row])], ignore_index=True)

        
    
    reference_layout['compound_concentrations'] = [len(combinations)]
    reference_layout['compound_concentration_names'] = [combinations]
    
    #for repl in [16, 14, 12, 10, 8, 6, 4]:
    for repl in range(14,4,-1):
        
        reference_layout['compound_replicates'] = [repl]
        #save_json(reference_layout, '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/plaid/reference_plate_etoh.json')
        #output = run_minizinc_command("/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/plaid/reference_plate_etoh.json")
        save_json(reference_layout, output_path)
        output = run_minizinc_command(output_path)
        
        if output == '=====UNSATISFIABLE=====\n':
            #print('Layout infeasible.')
            continue
        else:
            plate_layout = output_to_dataframe(output)
            print(plate_layout.shape[0])
            
            if 96-plate_layout.shape[0] < repl:
                continue
            
            print('Layout completed.')
    
            # Fill in plate layout
            #plate_layout = plate_filler(plate_layout)
            
            break
    
    plate_layout = plate_filler(plate_layout)
    plate_layout = plate_layout.merge(df, 
                                      left_on = 'CONCuM', 
                                      right_on = 'Experiment parameters', 
                                      how = 'left')[['plateID','well','CONCuM','Summary']]
    visualize_grid(plate_layout)
    
    return plate_layout
    




