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

    # Parse arguments
    args = parser.parse_args()
    
    # Return as a dictionary for easy access and logging
    return vars(args)

def main():
    
    #script_dir = os.path.dirname(os.path.abspath(__file__))
    #os.chdir(script_dir)
    
    # Parse arguments
    args = parse_arguments()

    # Log each variable setting
    for arg, value in args.items():
        print(f'Variable "{arg}" set to: {value}')

    # Example usage of parsed values
    folder = args['folder']
    T = args['T']
    
    # Set working folder
    
    os.chdir(f"{folder}")
    #print(os.getcwd())
    # Load OpenAI-key
    key = open("../../key.txt", "r")
    key = key.read()
    
    with open('hypothesis/feasibility/selection.json', 'r') as file:
        allowed_programs = json.load(file)
        
    full_prompt = open("hypothesis/patterns/prompt.txt", "r")
    full_prompt = full_prompt.read()

    allowed_numbers = []
    for program in allowed_programs['allowed']:
        allowed_numbers.append(program['number'])
    allowed_numbers = sorted(allowed_numbers)
    
    revised_prompt = extract_statements(full_prompt, allowed_numbers)
    revised_prompt = ' '.join(revised_prompt)
    
    #hypothesis_text = prompt_gpt_for_hypothesis(revised_prompt, T, key)
    #hypothesis_text = prompt_gpt_for_hypothesis('../../context/hypgen_context.txt', revised_prompt, T, key)
    hypothesis_text = prompt_gpt_for_hypothesis('../../context/hypgen_context.txt', revised_prompt, T, key)
    
    #print(hypothesis_text)

    # Save prompt
    with open('hypothesis/generated_hypothesis.txt', 'w') as file:
        file.write(hypothesis_text)


import json
import os
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
from hgen_support import *
from gpt_support import *

if __name__ == "__main__":
    main()
