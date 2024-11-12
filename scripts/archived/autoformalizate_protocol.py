#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:48:05 2024

@author: danbru
"""
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

    # Load OpenAI-key
    key = open("../key.txt", "r")
    key = key.read()

    # Example usage of parsed values
    folder = args['folder']
    os.chdir(folder)
    
    # Load hypothesis text
    hypothesis_text = open("hypothesis/generated_hypothesis.txt", "r")
    hypothesis_text = hypothesis_text.read()
    
    # Autoformalize hypothesis text into valid experimental protocol
    json_output, exp_plan = retry_experimental_plan('../../context/expdesign_context.txt',hypothesis_text, key, retries=3)
    
    if json_output == None:
        print('Faulty plan.')
        exit()
    
    # Write the content to the JSON file
    with open('protocol/protocol.json', 'w') as file:
        json.dump(json_output, file, indent=4)
    

import os
from hgen_support import *
from gpt_support import *
import json

if __name__ == "__main__":
    main()