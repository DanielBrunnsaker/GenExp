#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 11 13:57:16 2024

@author: danbru
"""



def run_script(command, input_text=None):
    """
    Helper function to run a command using subprocess. Allows passing input for stdin.
    """
    try:
        result = subprocess.run(command, check=True, shell=True, capture_output=True, text=True, input=input_text)
        return result.stdout.strip()  # Capture and return stdout output from the script
    except subprocess.CalledProcessError as e:
        print(f"Error running command: {e}")
        exit(1)

def main():
    
    script_dir = os.path.dirname(os.path.abspath(__file__))
    os.chdir(script_dir)
    GEN_EXP_ROOT_DIR = os.path.abspath("./..")
    os.environ["GEN_EXP_ROOT_DIR"] = GEN_EXP_ROOT_DIR
    
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Run all scripts in sequence with appropriate flags.")
    
    # Define flags that each script may require
    parser.add_argument('--target', required=True, help="Flag required for script1.")
    parser.add_argument('--alpha', required=True, help="Flag required for script2.")
    parser.add_argument('--N', required=True, help="Flag required for script3.")
    parser.add_argument('--T', required=False, help="Flag required for script4.")
    parser.add_argument('--volume', required=True, help="Flag required for script5.")
    
    parser.add_argument('--S', required=False, help="Additional input required by script5.")
    parser.add_argument('--Treatment', required=False, help="Additional input required by script5.")
    
    
    args = parser.parse_args()
    
    # Run Script 1 and capture the output folder name
    command1 = f"python pattern_selection.py --target {args.target} --alpha {args.alpha} --N {args.N}"
    #command1 = f"python pattern_selection.py --target glutamate --alpha 0.1 --N 10"
    output_folder = run_script(command1)
    output_folder = os.path.abspath(output_folder)
    print(f"Output folder from pattern_selector: {output_folder}")
    # Check if output_folder was successfully retrieved
    if not output_folder:
        print("Failed to retrieve output folder from pattern_selector")
        exit(1)

    ###########################################################################
    # Save config (moved here from `pattern_selection.py` script)
    main_folder = output_folder.rsplit("/", 1)[1]
    config_variables = {
        "folder": main_folder,
        "target": args.target,
        "alpha": args.alpha,
        "N": args.N,
        "T": args.T,
        "volume": args.volume
    }
    
    text_template = '{{"folder": "{folder}", "target": "{target}", "alpha": {alpha}, "N": "{N}", "T": "{T}", "volume": "{volume}"}}'
    formatted_text = text_template.format(**config_variables)
        
    # Convert the formatted string to a JSON-compatible Python object
    json_data = json.loads(formatted_text)
        
    # Save the JSON data to a file
    with open(f'../experiments/{main_folder}/config.json', 'w') as file:
        json.dump(json_data, file, indent=4) 
    ###########################################################################

    # Use the captured folder for subsequent scripts
    #command2 = f'python hypothesis_generation.py --T 0.7 --folder "{output_folder}"'
    print('Generating hypothesis... \n')
    command2 = f'python hypothesis_generation.py --T {args.T} --folder "{output_folder}"'
    run_script(command2)
    
    # Print hypothesis
    os.chdir(output_folder)
    hyp = open("hypothesis/generated_hypothesis.txt", "r")
    print('\n')
    print(hyp.read())
    print('\n')
    os.chdir(script_dir)
    
    print('Autoformalizing experimental plan... \n')
    command3 = f"python autoformalize_protocol.py --folder {output_folder}"
    run_script(command3)

    action = user_prompt()
    
    if action == "go_ahead":
        print('Proceeding with experimental design...')
        
    elif action == "retry":
        print("Retrying the previous operation. \n")
        
        # Specify the path to your file
        file_path = "../experiments/forbidden_patterns.txt"
        
        with open(f'{output_folder}/protocol/protocol.json', 'r') as file:
            data = json.load(file)
        
        with open(file_path, "r") as file:
            lines = file.readlines()
            lines = [line.strip() for line in lines]  # Strip newline characters
        
        new_line = data['logic_program'] # Replace with the text you want to add
        lines.append(new_line)
        
        with open(file_path, "w") as file:
            for line in lines:
                file.write(line + "\n")  # Write each item followed by a newline
        
        shutil.rmtree(output_folder)
        main()  # Call the function again to simulate retry
    
    elif action == "reprompt":
        print("Repromptng.")
        shutil.rmtree(output_folder)
        # Exit or stop the execution
        main()
    
    elif action == "cancel":
        print("Exiting the script.")
        shutil.rmtree(output_folder)
        # Exit or stop the execution
        exit()


    print('Designing plate layout... \n')
    command4 = f"python plate_layout.py --folder {output_folder}"
    run_script(command4)

    # Run script5 with optional input
    input_text_S = args.S if args.S else input("Please provide stock concentration for supplement 1 if applicable: ")
    input_text_T = args.Treatment if args.Treatment else input("Please provide stock concentration for treatment 1 if applicable: ")
    
    print('Generating hamilton runlists... \n')
    command5 = f"python hamilton_protocol.py --volume {args.volume} --folder {output_folder} --S {input_text_S} --Treatment {input_text_T}"
    run_script(command5)

    print("All steps completed successfully.")

import json
import subprocess
import argparse
import os
import shutil

script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
from hgen_support import *

if __name__ == "__main__":
    main()