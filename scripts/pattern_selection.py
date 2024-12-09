
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug 30 15:24:19 2024

@author: danbru
"""


''' 
    Hypothesis selection -- this part needs to be automated!
    
    How do we select the hypothesis?
    
    Input: Coefficients, along with logic programs
    
    TODO:
        Filter the hypotheses that are not actionable
        Filter the hypotheses that are infeasible
        Filter away toxicants etc.
        How? Simple rule-based system?
    
    Options? 
        A. Automatically select the one with the highest magnitude?
        B. Present several and evaluate after somehow?
'''

def create_folder_structure_recursively(base_path, folder_structure):
    import os

    def create_folders_recursively(current_path, structure):
        for folder, sub_structure in structure.items():
            # Create the current folder
            folder_path = os.path.join(current_path, folder)
            if not os.path.exists(folder_path):
                os.makedirs(folder_path)
                #print(f"Created folder: {folder_path}")

            # If the current folder has subfolders, create them recursively
            if isinstance(sub_structure, dict):
                create_folders_recursively(folder_path, sub_structure)

    # Start creating folders from the base path
    create_folders_recursively(base_path, folder_structure)

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Hypothesis Generator with adjustable parameters.")
    
    # Define arguments with default values
    parser.add_argument('--target', type=str, default="alanine", help='Predictive target (e.g., amino acid).')
    parser.add_argument('--N', type=int, default=10, help='Number of patterns to pass to the initial prompt.')
    parser.add_argument('--alpha', type=float, default=0.1, help='Parameter for weighing the linear coefficients for pattern ranking.')

    # Parse arguments
    args = parser.parse_args()
    
    # Return as a dictionary for easy access and logging
    return vars(args)

def main():
    
    
    '''
    default values
    python h_generator.py --target "alanine" --N 10 --T 1.0 --alpha 1.0
        
    '''
    # Parse arguments
    args = parse_arguments()

    # Log each variable setting
    #for arg, value in args.items():
    #    print(f'Variable "{arg}" set to: {value}')

    # Example usage of parsed values
    target = args['target']
    N = args['N']
    alpha = args['alpha']
    
    
    # Load OpenAI-key
    key = open("../key.txt", "r")
    key = key.read()
    
    
    
    
    # Create a new folder?
    base_path = '../experiments' 
    main_folder = f'{target}_{alpha}_{N}_{datetime.now().strftime("%Y%m%d_%H%M")}'
    #subfolders_structure = {
    #    'hypothesis': ['generated_hypotheses', 'selected_hypothesis'],
    #    'protocol': ['hamilton', 'EVE', 'plate_layout'],
    #    'results': ['growth', 'metabolomics']
    #}
    
    
    base_path = '../experiments' 
    folder_structure = {
        main_folder: {
            "hypothesis": {
                "generated_hypotheses": {
                    "initial_stage": {},
                    "second_stage": {}
                },
                "selected_hypothesis": {}
            },
            "protocol": {
                "EVE": {},
                "hamilton": {},
                "plate_layout": {},
            },
            "results": {
                "Subsubfolder1": {}
            }
            
        }
    }
    
    create_folder_structure_recursively(base_path, folder_structure)
    #create_folder_structure(base_path, main_folder, subfolders_structure)
    
    
    relevance_scores = {
        # predicates
        'compound_name': 10, # Reward patterns with a specific condition in mind (if there is one)
        'condition': 10, # Reward patterns with a specific chemical in mind (if there is one)
        'participates_in_metabolism': 5, 
        'compound_modulates_target': 3, # add production? should make more sense as we are measuring accumulations
        'production of': 2,
        'consumption of': 1
    }
    
    # Load supports. Directly save as dict instead?
    id_dict = pd.read_csv('../data/uniprot_id_map.tsv', sep = '\t')
    id_dict['To'] = id_dict['To'].str.replace('sce:','')
    id_dict = id_dict.set_index('From')['To'].to_dict()
    orf_dict = pd.read_csv('../data/orf_dict.tsv', sep = '\t', header = None)[[1,3]]
    orf_dict = orf_dict.set_index(1)[3].to_dict()

    
    # Load reference data
    aaSet = pd.read_excel('../data/AA.xls', sheet_name = 'intracellular_concentration_mM').set_index('ORF').iloc[:,1:]
    #path = '../results/coefficients/20241025/'
    path = '../results/coefficients/20241114/'
    
    # Calculate scoring for hypotheses
    hypotheses, all_hypotheses = load_all_hypotheses(path, target, aaSet)    
    sorted_specifity = calculate_scores(all_hypotheses, target, alpha=alpha) # Calculates ordering of hypothesis based on selected alpha (i.e. penalize patterns not specific to one target)
    filtered_sorted_specifity_all = sorted_specifity[sorted_specifity[target] != 0]
    filtered_sorted_specifity = sorted_specifity[target][sorted_specifity[target] != 0]
    
    hypothesis_counter = 0
    full_prompt = ''
    for counter, hypothesis in enumerate(filtered_sorted_specifity.index):
        
        #ind_statement = ''
        directionality = 'higher' if np.sign(hypotheses.loc[hypothesis]) == 1.0 else 'lower'
        #break
        # Find potential alternative hypotheses?
        with open('../results/patterns/feature_dict.pickle', 'rb') as handle:
            duplicate_dict = pickle.load(handle)
        
        
        # Go through alternatives, and select the best one based on the scorelist given.
        if len(duplicate_dict[hypothesis]) != 0:
            score = -1
            for item in duplicate_dict[hypothesis]:
                temp_expl = find_row_with_string(f"../results/patterns/explanations/{item.split('_')[0]}.txt", f"({item.split('_')[1].split('p')[1]},")
                if score < calculate_string_score(temp_expl, relevance_scores):
                    score = calculate_string_score(temp_expl, relevance_scores)
                    current_clause = item
                    
        # In case of no alternative formulations
        else:
            current_clause = hypothesis
            
        first_example = current_clause.split('_')[0]
        h_explanation = find_row_with_string(f'../results/patterns/explanations/{first_example}.txt', f"({current_clause.split('_')[1].split('p')[1]},")
        clause = h_explanation.split(',(')[1][:-3].replace('gene(A):-','Cell(A):-')
        
        # Check for disallowed statements
        if any(clause.replace("'","") in s for s in disallowed_clauses):
            continue
        
        else:
            # Produce the custom output:
            #print(clause)
            hypothesis_counter += 1
            #prompt_part = f"{counter}. Cells with {directionality} than normal levels of intracellular {target} in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: "            
            prompt_part = f"{counter}. Cells with {directionality} levels of intracellular {target} in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: "            
            
            ind_statement = prompt_part + clause
            #full_prompt = full_prompt + prompt_part + clause + '. '
            
            hypothesis_text = prompt_gpt_for_hypothesis('../context/hypgen_per_pattern.txt', ind_statement, 0.5, key)
            
            # Save hypothesis, along with the pattern and the clause?
            #print(hypothesis_text)
            
            # Save details
            rev_clause = clause.replace("'","")
            details = {"logic_program": rev_clause, 
                       "observable": target, 
                       "qualifier": directionality, 
                       "coefficient": filtered_sorted_specifity[hypothesis], 
                       "initial_prompt": ind_statement}
            
            #details = f'{{"logic_program": "{rev_clause}", "observable": "{target}", "qualifier": "{directionality}", "coefficient": {filtered_sorted_specifity[hypothesis]}}}'
            
            # Save prompt../experiments/{main_folder}/
            with open(f'../experiments/{main_folder}/hypothesis/generated_hypotheses/initial_stage/hypothesis_{hypothesis_counter}.txt', 'w') as file:
                file.write(hypothesis_text)
                
            with open(f'../experiments/{main_folder}/hypothesis/generated_hypotheses/initial_stage/hypothesis_details_{hypothesis_counter}.json', 'w') as file:
                #json.dump(json.loads(details), file, indent=4)
                json.dump(details, file, indent=4)
            
            # Restrict to N generated clauses
            if hypothesis_counter == int(N):
                break
    
    
    print(f"{base_path}/{main_folder}")
    
    
    '''
    
    # Create a new folder?
    base_path = '../experiments' 
    main_folder = f'{target}_{alpha}_{N}_{datetime.now().strftime("%Y%m%d_%H%M")}'
    subfolders_structure = {
        'hypothesis': ['feasibility', 'patterns'],
        'protocol': ['hamilton', 'EVE', 'plate_layout'],
        'results': ['growth', 'metabolomics']
    }

    create_folder_structure(base_path, main_folder, subfolders_structure)
    
    # Filter the space of hyotheses provided by the prolog-programs using GPT.
    allowed_programs_result = safety_feasibility_prompt(full_prompt, key)
    
    # Attempt to load the JSON
    allowed_programs = json.loads(allowed_programs_result)
    
    # Write the content to the JSON file
    with open(f'../experiments/{main_folder}/hypothesis/feasibility/selection.json', 'w') as file:
        json.dump(allowed_programs, file, indent=4)
    
    # Save prompt
    with open(f'../experiments/{main_folder}/hypothesis/patterns/prompt.txt', 'w') as file:
        file.write(full_prompt)
        
    # Save config
    config_variables = {
        "folder": main_folder,
        "alpha": alpha,
        "N": N
    }
    
    text_template = '{{"folder": "{folder}", "alpha": {alpha}, "N": "{N}"}}'
    formatted_text = text_template.format(**config_variables)
        
    # Convert the formatted string to a JSON-compatible Python object
    json_data = json.loads(formatted_text)
        
    # Save the JSON data to a file
    with open(f'../experiments/{main_folder}/config.json', 'w') as file:
        json.dump(json_data, file, indent=4)  
    
    print(f"{base_path}/{main_folder}")
    '''

from datetime import datetime
import pandas as pd
import json
import pickle
import os
import numpy as np
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
from utils.hgen_support import *
from utils.gpt_support import *
#from layout_generator import *
#from dose_calculator2 import *

with open('../experiments/forbidden_patterns.txt', "r") as file:
    disallowed_clauses = file.readlines()
    
disallowed_clauses = [line.strip() for line in disallowed_clauses]

with open('../experiments/succesful_patterns.txt', "r") as file:
    successful_clauses = file.readlines()
    
successful_clauses = [line.strip() for line in successful_clauses]
disallowed_clauses = disallowed_clauses + successful_clauses

if __name__ == "__main__":
    main()





#os.chdir('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/scripts')


