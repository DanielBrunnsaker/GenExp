#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:26:00 2024

@author: danbru
"""


def parse_arguments():
    import argparse
    parser = argparse.ArgumentParser(description="Hypothesis Generator with adjustable parameters.")
    
    # Define arguments with default values
    parser.add_argument('--folder', type=str, help='Working folder for experiment.')
    parser.add_argument('--T', type=float, default=0.7, help='Temperature for the hypothesis generation step.')
    parser.add_argument('--N', type=int, help='Number of hypotheses to choose from.')

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
    T = args['T']
    N = args['N']
    # Set working folder
    
    os.chdir(f"{folder}")
    # Load OpenAI-key
    key = open("../../key.txt", "r")
    key = key.read()
    
    complete_prompt = ''
    for nr_hyp in range(1,N+1):
        complete_prompt += f'####### GENERATED HYPOTHESIS {nr_hyp} #######\n'
        hyp_text = open(f"hypothesis/generated_hypotheses/initial_stage/hypothesis_{nr_hyp}.txt", "r").read()
        complete_prompt += hyp_text +'\n\n'
        
    final_selection_text = prompt_gpt_for_hypothesis('../../context/hypothesis_selection.txt', complete_prompt, T, key)
    print(final_selection_text)
    
    # Extract the hypothesis number
    N_extracted = int(re.search(r"HYPOTHESIS NUMBER #(\d+)", final_selection_text).group(1))
   
    selected_hypothesis = open(f"hypothesis/generated_hypotheses/initial_stage/hypothesis_{N_extracted}.txt", "r").read()
    hypothesis_details = load_json(f"hypothesis/generated_hypotheses/initial_stage/hypothesis_details_{N_extracted}.json")
    
    # Add the motivation behind selecting it
    hypothesis_details['motivation'] = final_selection_text
    with open('hypothesis/selected_hypothesis/hypothesis_details.json', 'w') as file:
        json.dump(hypothesis_details, file, indent=4)

    variants = '####### OPTION 1 #######\n\n' + selected_hypothesis
    with open('hypothesis/generated_hypotheses/second_stage/hypothesis_1.txt', 'w') as file:
        file.write(selected_hypothesis)
    
    for i, variants_iter in enumerate(range(2,4)):
        variants += f'\n\n####### OPTION {variants_iter} #######\n'
        new_variant = prompt_gpt_for_hypothesis('../../context/hypgen_per_pattern.txt', hypothesis_details['initial_prompt'], T, key)
        with open(f'hypothesis/generated_hypotheses/second_stage/hypothesis_{variants_iter}.txt', 'w') as file:
            file.write(new_variant)
        
        variants += new_variant
        
    with open('hypothesis/selected_hypothesis/variants.txt', 'w') as file:
        file.write(variants)
    
    final_selection_text = prompt_gpt_for_hypothesis('../../context/refine_plan.txt', variants, T, key)
    with open('hypothesis/selected_hypothesis/reasoning.txt', 'w') as file:
        file.write(final_selection_text)
    
    N_extracted = int(re.search(r"HYPOTHESIS NUMBER #(\d+)", final_selection_text).group(1))
    finalized_hypothesis = open(f"hypothesis/generated_hypotheses/second_stage/hypothesis_{N_extracted}.txt", "r").read()
    
    with open('hypothesis/selected_hypothesis/hypothesis.txt', 'w') as file:
        file.write(finalized_hypothesis)

        
    
    
    

import re
import json
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
from utils.hgen_support import *
from utils.gpt_support import *

if __name__ == "__main__":
    main()
