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

def run_script_with_retry(command, input_text=None, retries=1):
    
    import time
    
    """
    Executes a command with a specified number of retries. Allows passing input for stdin.
    Captures and returns stdout output from the script.
    Args:
        command (str): The command to run.
        input_text (str, optional): Text to send to stdin.
        retries (int): Number of retries before giving up.
    Returns:
        str: The captured stdout output from the script.
    """
    for attempt in range(retries + 1):
        try:
            print(f"Running command: {command} (Attempt {attempt + 1})")
            result = subprocess.run(command, check=True, shell=True, capture_output=True, text=True, input=input_text)
            return result.stdout.strip()  # Capture and return stdout output
        except subprocess.CalledProcessError as e:
            print(f"Error occurred: {e}")
            if attempt < retries:
                print("Retrying...")
                time.sleep(2)  # Optional delay between retries
            else:
                print("All retries failed. Exiting...")
                exit(1) 


def main():
    
    '''
    Example: python hgen.py --target alanine --N 10 --alpha 0.5 --volume 270 --T 0.5 --variations 5

    '''
    
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Run all scripts in sequence with appropriate flags.")
    
    parser.add_argument('--target', required=True, help="Flag required for script1.")
    parser.add_argument('--alpha', required=True, help="Flag required for script2.")
    parser.add_argument('--variations', type=int, default=3, help='Flag required for script2.')
    parser.add_argument('--N', required=True, help="Flag required for script3.")
    #parser.add_argument('--T', required=False, help="Flag required for script4.")
    #parser.add_argument('--volume', required=True, help="Flag required for script5.")
    
    parser.add_argument('--supplement_stock', required=False, help="Additional input required by script5.")
    parser.add_argument('--treatment_stock', required=False, help="Additional input required by script5.")
    parser.add_argument('--negative_control_stock', required=False, help="Additional input required by script5.")
    
    args = parser.parse_args()
    
    
    config_file = {}
    config_file['target'] = args.target
    config_file['alpha'] = args.alpha
    config_file['variations'] = args.variations
    config_file['N'] = args.N
    
    config_file['T'] = config.LLM_TEMPERATURE
    config_file['volume'] = config.WELL_VOLUME
    #config_file['T'] = args.T
    #config_file['volume'] = args.volume
    
    
    '''
        Run hypothesis selection step
    '''
    
    command1 = f"python pattern_selection.py --target {args.target} --alpha {args.alpha} --N {args.N}"
    output_folder = run_script(command1)
    output_folder = os.path.abspath(output_folder.strip())
    print(f"Output folder from pattern_selector: {output_folder}")
    # Check if output_folder was successfully retrieved
    if not output_folder:
        print("Failed to retrieve output folder from pattern_selector")
        exit(1)


    '''
        Run hypothesis generation step
    '''
    # Use the captured folder for subsequent scripts
    print('Generating hypothesis... \n')
    command2 = f'python hypothesis_generation.py --T {args.T} --N {args.N} --variations {args.variations} --folder "{output_folder}"'
    run_script_with_retry(command2, retries=1)
    
    # Print hypothesis
    os.chdir(output_folder)
    hyp = open("hypothesis/selected_hypothesis/hypothesis.txt", "r")
    print('\n')
    print(hyp.read())
    print('\n')
    os.chdir(script_dir)
    
    '''
        Autoformalize 
    '''
    
    print('Autoformalizing experimental plan... \n')
    command3 = f"python autoformalize_protocol.py --folder {output_folder}"
    run_script_with_retry(command3, retries=1)

    action = user_prompt()
    
    if action == "go_ahead":
        print('Proceeding with experimental design...')
        
    elif action == "forbidding pattern":
        
        # Specify the path to your file
        file_path = "../experiments/forbidden_patterns.txt"
                
        with open(f'{output_folder}/hypothesis/selected_hypothesis/hypothesis_details.json', 'r') as file:
            data = json.load(file)
        
        with open(file_path, "r") as file:
            lines = file.readlines()
            lines = [line.strip() for line in lines]  # Strip newline characters
        
        new_line = f'{args.target}:'+data['logic_program'] # Replace with the text you want to add
        lines.append(new_line)
        
        with open(file_path, "w") as file:
            for line in lines:
                file.write(line + "\n")  # Write each item followed by a newline
        
        shutil.rmtree(output_folder)
        exit()
        #main()  # Call the function again to simulate retry
    
    elif action == "reprompt":
        print("Reprompting.")
        shutil.rmtree(output_folder)
        # Exit or stop the execution
        main()
    
    elif action == "cancel":
        print("Exiting the script.")
        shutil.rmtree(output_folder)
        # Exit or stop the execution
        exit()


    '''
        Set up technical details regarding experimental run 
    '''



    











    # Read the involved treatments
    protocol = load_json(f'{output_folder}/protocol/protocol.json')
    if protocol is None:
        raise ValueError("Protocol-JSON could not be loaded.")
        
    first_experiment = next((exp for exp in protocol['experiments'] if exp.get('type', '').lower() == 'experiment'), None)
    neg_experiment = [entry for entry in  protocol['experiments'] if 'type' in entry and 'negative' in entry['type'].lower()][0]

    supplement = first_experiment['media_supplementation']
    treatment = first_experiment['treatment']
    negative_control = neg_experiment['media_supplementation']

    # Add this one to the completed ones?
    successful_program = output_folder+'/hypothesis/selected_hypothesis/hypothesis_details.json'
    with open(successful_program, 'r') as file:
        lprog = json.load(file)['logic_program']
    
    if 'temperature' in lprog or '°' in lprog.lower(): # Change this so that it is more general
        print('Designing plate layout... \n')
        command4 = f"python plate_layout_environmental.py --folder {output_folder}"
        run_script(command4)
        
        input_text_S = args.supplement_stock if args.supplement_stock else input(f"Please provide {supplement} stock concentration, if applicable (If not in library, aim for minimum of 4.5 times the conc.): ")
        input_text_n = args.negative_control_stock if args.negative_control_stock else input(f"Please provide {negative_control} stock concentration, if applicable (If not in library, aim for minimum of 4.5 times the conc.): ")
        
        print('Generating hamilton runlists... \n')
        command5 = f'python hamilton_protocol.py --volume {args.volume} --folder {output_folder} --supplement_stock {input_text_S} --treatment_stock "0 mM" --negative_control_stock {input_text_n} --environmental 1'
        run_script(command5)
        
        config_file['supplement_stock'] = input_text_S
        config_file['negative_control_stock'] = input_text_n
        
    else:
        print('Designing plate layout... \n')
        command4 = f"python plate_layout.py --folder {output_folder}"
        run_script(command4)
        
        input_text_S = args.supplement_stock if args.supplement_stock else input(f"Please provide {supplement} stock concentration, if applicable (If not in library, aim for minimum of 4.5 times the conc.): ")
        input_text_n = args.negative_control_stock if args.negative_control_stock else input(f"Please provide {negative_control} stock concentration, if applicable (If not in library, aim for minimum of 4.5 times the conc.): ")
        input_text_T = args.treatment_stock if args.treatment_stock else input(f"Please provide {treatment} stock concentration, if applicable (If not in library, aim for minimum of 4.5 times the conc.): ")
        
        print('Generating hamilton runlists... \n')
        command5 = f"python hamilton_protocol.py --volume {args.volume} --folder {output_folder} --supplement_stock {input_text_S} --treatment_stock {input_text_T} --negative_control_stock {input_text_n}"
        run_script(command5)
        
        config_file['supplement_stock'] = input_text_S
        config_file['treatment_stock'] = input_text_T
        config_file['negative_control_stock'] = input_text_n
        
    with open('../experiments/succesful_patterns.txt', "r") as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]  # Strip newline characters
    
    lines.append(f"{args.target}:"+lprog)
    
    with open('../experiments/succesful_patterns.txt', "w") as file:
        for line in lines:
            file.write(line + "\n")  # Write each item followed by a newline    
    
    with open(f'{output_folder}/config.json', 'w') as f:
        json.dump(config_file, f)
    
    print(f"Experimental plan for '{lprog}' completed successfully.")
    
    
    
    

import json
import subprocess
import argparse
import os
import shutil
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
from utils.hgen_support import *
import config as config

if __name__ == "__main__":
    main()
    