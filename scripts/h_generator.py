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

def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Hypothesis Generator with adjustable parameters.")
    
    # Define arguments with default values
    parser.add_argument('--target', type=str, default="alanine", help='Predictive target (e.g., amino acid).')
    parser.add_argument('--N', type=int, default=25, help='Number of patterns to pass to the initial prompt.')
    parser.add_argument('--T', type=float, default=0.5, help='LLM temperature for hypothesis generation step.')
    parser.add_argument('--alpha', type=float, default=1.0, help='Parameter for weighing the linear coefficients for pattern ranking.')
    parser.add_argument('--beta', type=float, default=0.20, help='Parameter that penalizes the pattern ranking through the summation of alternative targets.')

    # Parse arguments
    args = parser.parse_args()
    
    # Return as a dictionary for easy access and logging
    return vars(args)

def main():
    
    
    '''

    default values
    python h_generator.py --target "alanine" --N 10 --T 1.0 --alpha 1.0 --beta 0.10
        
    '''
    # Parse arguments
    args = parse_arguments()

    # Log each variable setting
    for arg, value in args.items():
        print(f'Variable "{arg}" set to: {value}')

    # Example usage of parsed values
    target = args['target']
    N = args['N']
    T = args['T']
    alpha = args['alpha']
    beta = args['beta']
    
    
    # Load this in instead, as a JSON?
    relevance_scores = {
        # predicates
        'compound_name': 10, # Reward patterns with a specific condition in mind (if there is one)
        'condition': 10, # Reward patterns with a specific chemical in mind (if there is one)
        'interacts_with_metabolite': 5, 
        'biological_target': 3,
        ## Add perturbation?
    }
    
    # Load supports
    id_dict = pd.read_csv('../data/uniprot_id_map.tsv', sep = '\t')
    id_dict['To'] = id_dict['To'].str.replace('sce:','')
    id_dict = id_dict.set_index('From')['To'].to_dict()
    orf_dict = pd.read_csv('../data/orf_dict.tsv', sep = '\t', header = None)[[1,3]]
    orf_dict = orf_dict.set_index(1)[3].to_dict()

    key = open("../key.txt", "r")
    key = key.read()

    # Load reference data
    aaSet = pd.read_excel('../data/AA.xls', sheet_name = 'intracellular_concentration_mM').set_index('ORF').iloc[:,1:]
    path = '../results/coefficients/20241025/'
    
    hypotheses, all_hypotheses = load_all_hypotheses(path, target, aaSet)
    
    # Alpha = 1.0, beta = 0.25 seem reasonable starting points
    sorted_specifity = calculate_scores(all_hypotheses, target, alpha=alpha, beta=beta)
    filtered_sorted_specifity_all = sorted_specifity[sorted_specifity[target] != 0]
    filtered_sorted_specifity = sorted_specifity[target][sorted_specifity[target] != 0]
    
    print('\n Generating logic programs... \n')
    
    hypothesis_counter = 0
    full_prompt = ''
    for counter, hypothesis in enumerate(filtered_sorted_specifity.index):

        directionality = 'higher' if np.sign(hypotheses.loc[hypothesis]) == 1.0 else 'lower'
        
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
        if any(clause in s for s in disallowed_clauses):
            continue
        
        else:
            # Produce the custom output:
            hypothesis_counter += 1
            prompt_part = f"{counter}. Cells with {directionality} than normal levels of intracellular {target} in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: "
            
            #prompt_part = f"{hypothesis_counter}. Cells with {directionality} than normal levels of intracellular {target} in standard conditions (grown on minimal media without any amino acids) associate with the following phenotype, described as a prolog program: "
            full_prompt = full_prompt + prompt_part + clause + '. '
            
            # Restrict to N generated clauses
            if hypothesis_counter == int(N):
                break
    
    # Filter the space of hyotheses provided by the prolog-programs using GPT.
    allowed_programs_result = safety_feasibility_prompt(full_prompt, key)
    
    # Attempt to load the JSON
    allowed_programs = json.loads(allowed_programs_result)
    
    allowed_numbers = []
    for program in allowed_programs['allowed']:
        allowed_numbers.append(program['number'])
    allowed_numbers = sorted(allowed_numbers)
    
    
    #filtered_sorted_specifity_all
    revised_prompt = extract_statements(full_prompt, allowed_numbers)
    
    revised_prompt = ' '.join(revised_prompt)
    
    #print(revised_prompt)
    
    hypothesis_text = prompt_gpt_for_hypothesis(revised_prompt, T, key)
    print(hypothesis_text)
    
    print('\n ----------------------------------------- \n ')
    
    json_output, exp_plan = retry_experimental_plan(hypothesis_text, key,retries=3)
    
    if json_output == None:
        print('Faulty plan.')
        exit()
    
    pretty_json = json.dumps(json_output, indent=2)
    print(pretty_json)
    
    return_path = save_string_to_json('../experiments/generated_outputs', json_output, 'protocol', target, N, T, alpha, beta)
    save_string_to_json('../experiments/generated_outputs', allowed_programs, 'selection', target, N, T, alpha, beta)
    save_string_to_file('../experiments/generated_outputs/', hypothesis_text, 'hypothesis', target, N, T, alpha, beta)
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    # Present user with the options
    action = user_prompt()
    
    
    
    if action == "go_ahead":
        #print("Continuing with the rest of the script.")
        # Continue with the rest of the script
        
        # Construct the filename
        plate_layout = generate_layout(return_path, '../plaid/reference_plate.json', f"../experiments/generated_outputs/{timestamp}_layout_{target}_{N}_{T}_{alpha}_{beta}.json")
        plate_save_path = f"../experiments/generated_outputs/{timestamp}_layoutTable_{target}_{N}_{T}_{alpha}_{beta}.tsv"
        plate_layout.to_csv(plate_save_path, sep = '\t')
        
        #print(plate_layout)
        return_hamilton_concentrations(return_path, plate_save_path)
        
    elif action == "retry":
        print("Retrying the previous operation.")
        # You can loop back to retry the operation
        disallowed_clauses.append(json_output['logic_program'])
        main()  # Call the function again to simulate retry
        
    elif action == "cancel":
        print("Exiting the script.")
        # Exit or stop the execution
        exit()
    
    
from datetime import datetime
import pandas as pd
import json
import pickle
import os
import numpy as np

#os.chdir('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/scripts')

#os.chdir(os.getcwd())

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)

from hgen_support import *
from gpt_support import *
from layout_generator import *
from dose_calculator2 import *

# Set starting point of the disallowed formulations. Note that the starting point here are all logic programs 
# which have no meaning without a second, descriptive, predicate/atom.
disallowed_clauses = ["Cell(A):-phenotype(A,'increased chemical compound accumulation',B,C)",
              "Cell(A):-phenotype(A,'decreased chemical compound accumulation',B,C)",
              "Cell(A):-phenotype(A,'increased chemical compound excretion',B,C)",
              "Cell(A):-phenotype(A,'decreased chemical compound excretion',B,C)",
              "Cell(A):-phenotype(A,'increased resistance to chemicals',B,C)",
              "Cell(A):-phenotype(A,'decreased resistance to chemicals',B,C)",
              "Cell(A):-phenotype(A,'decreased chemical compound accumulation',B,C),compound_name(B,proton)",
              "Cell(A):-phenotype(A,'increased chemical compound accumulation',B,C),compound_name(B,proton)."]

if __name__ == "__main__":
    main()








