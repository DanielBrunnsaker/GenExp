#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:59:51 2024

@author: danbru
"""

from pathlib import Path
import pandas as pd

from utils.plate_utils import *

def generate_layout(output_folder):
    
    EXPERIMENT_DIR = Path(output_folder) 
   
    json_data = load_json(EXPERIMENT_DIR / 'protocol/protocol.json')
    reference_layout = load_json(EXPERIMENT_DIR / '../../plaid/reference_plate.json')
    output_path = EXPERIMENT_DIR / 'protocol/plate_layout/minizinc_reference.json'
    
    experiment_entry = [entry for entry in  json_data['experiments'] if 'type' in entry and 'negative' in entry['type'].lower()][0]
    negcontrol_name = experiment_entry['media_supplementation'].replace(",","_")
    
    experiment_entry = [entry for entry in  json_data['experiments'] if 'type' in entry and 'experiment' in entry['type'].lower()][0]
    treatment_name = experiment_entry['treatment'].replace(",","_")
    supplementation_name = experiment_entry['media_supplementation'].replace(",","_")
    
    # If the treatment/supplement/control names have a comma or something in them, fix
    combinations = []
    
    df = pd.DataFrame(columns = ['Summary','Experiment Parameters'])
    
    for exp in json_data['experiments']:
        
        
        supplementation_dose = exp['media_supplementation_doses']
        treatment_dose = exp['treatment_parameters']
        
        if supplementation_dose == None:
            supplementation_dose = '0 mM'
        
        if treatment_dose == None:
            treatment_dose = '0 mM'
        
        supp = f"{exp['media_supplementation']}: "+supplementation_dose
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
            
            if 96-plate_layout.shape[0] < repl:
                continue
            break
    
    plate_layout = plate_filler(plate_layout)
    plate_layout = plate_layout.merge(df, 
                                      left_on = 'CONCuM', 
                                      right_on = 'Experiment parameters', 
                                      how = 'left')[['plateID','well','CONCuM','Summary']]

    plate_save_path = EXPERIMENT_DIR / "protocol/plate_layout/layout.tsv"
    plate_layout.to_csv(plate_save_path, sep = '\t')
    
if __name__ == "__main__":
    
    import argparse
    
    parser = argparse.ArgumentParser(description="Hypothesis generation")
    parser.add_argument("--output_folder", required=True, type=str, help="Experiment folder")
    
    args = parser.parse_args()
    generate_layout(args.output_folder)
