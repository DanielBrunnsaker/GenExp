#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:59:51 2024

@author: danbru
"""

# Function to load a JSON file


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
        
        print(output_path)
        
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


    folder = args['folder'] # Ugly solution, fix at some point
    os.chdir(folder)
    plate_save_path = "protocol/plate_layout/layout.tsv"
    plate_layout.to_csv(plate_save_path, sep = '\t')
    

import os
import pandas as pd
import io
from utils.plate_utils import *

if __name__ == "__main__":
    main()
