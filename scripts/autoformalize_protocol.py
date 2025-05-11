#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 10:48:05 2024

@author: danbru
"""

import os
from utils.hgen_support import *
from utils.gpt_support import *
import json
from pathlib import Path

BASE_DIR = Path(__file__).resolve().parent

def autoformalize(output_folder):
    
    # Load OpenAI-key
    key = open(BASE_DIR / "../key.txt", "r")
    key = key.read()
    
    EXPERIMENT_DIR = Path(output_folder)
    
    # Load hypothesis text
    hypothesis_text = open(EXPERIMENT_DIR / "hypothesis/selected_hypothesis/hypothesis.txt", "r")
    hypothesis_text = hypothesis_text.read()
    
    # Autoformalize hypothesis text into valid experimental protocol
    json_output, exp_plan, completion = retry_experimental_plan(EXPERIMENT_DIR / '../../context/expdesign_context.txt', hypothesis_text, key, retries=3)
    with open(Path(output_folder) / "versions/v_autoformalization.txt", "w") as file:
        file.write(completion.model)
      
    # Replace with the correct logic program
    with open(EXPERIMENT_DIR / 'hypothesis/selected_hypothesis/hypothesis_details.json', 'r') as file:
        details = json.load(file)
    
    json_output['logic_program'] = details['logic_program']
    
    # Write the content to the JSON file
    with open(EXPERIMENT_DIR / 'protocol/protocol.json', 'w') as file:
        json.dump(json_output, file, indent=4)
    

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description="Autoformalization")
    parser.add_argument("--output_folder", required=True, type=str, help="Experiment folder")
    
    args = parser.parse_args()
    
    autoformalize(args.output_folder)