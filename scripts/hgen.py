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
    Example: python hgen.py --target alanine --N 5 --alpha 0.25 --volume 270 --T 0.5

    '''
    
    #script_dir = os.path.dirname(os.path.abspath(__file__))
    #os.chdir(script_dir)
    
    # Set up argument parsing
    parser = argparse.ArgumentParser(description="Run all scripts in sequence with appropriate flags.")
    
    parser.add_argument('--target', required=True, help="Flag required for script1.")
    parser.add_argument('--alpha', required=True, help="Flag required for script2.")
    parser.add_argument('--N', required=True, help="Flag required for script3.")
    parser.add_argument('--T', required=False, help="Flag required for script4.")
    parser.add_argument('--volume', required=True, help="Flag required for script5.")
    
    parser.add_argument('--S', required=False, help="Additional input required by script5.")
    parser.add_argument('--Treatment', required=False, help="Additional input required by script5.")
    
    args = parser.parse_args()
    
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
    command2 = f'python hypothesis_generation.py --T {args.T} --N {args.N} --folder "{output_folder}"'
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
        
        new_line = data['logic_program'] # Replace with the text you want to add
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

    # Add this one to the completed ones?
    successful_program = output_folder+'/hypothesis/selected_hypothesis/hypothesis_details.json'
    with open(successful_program, 'r') as file:
        lprog = json.load(file)['logic_program']
    
    if 'temperature' in lprog or '°' in lprog.lower(): # Change this so that it is more general
        print('Designing plate layout... \n')
        command4 = f"python plate_layout_environmental.py --folder {output_folder}"
        run_script(command4)
        
        input_text_S = args.S if args.S else input("Please provide stock concentration for supplement 1 if applicable (If not in library, aim for minimum of 4.5 times the conc.): ")
        
        print('Generating hamilton runlists... \n')
        command5 = f'python hamilton_protocol.py --volume {args.volume} --folder {output_folder} --S {input_text_S} --Treatment "0 mM" --environmental 1'
        run_script(command5)
    else:
        print('Designing plate layout... \n')
        command4 = f"python plate_layout.py --folder {output_folder}"
        run_script(command4)
        
        input_text_S = args.S if args.S else input("Please provide stock concentration for supplement 1 if applicable (If not in library, aim for minimum of 4.5 times the conc.): ")
        input_text_T = args.Treatment if args.Treatment else input("Please provide stock concentration for treatment 1 if applicable (If not in library, aim for minimum of 4.5 times the conc.): ")
        
        print('Generating hamilton runlists... \n')
        command5 = f"python hamilton_protocol.py --volume {args.volume} --folder {output_folder} --S {input_text_S} --Treatment {input_text_T}"
        run_script(command5)
        
    with open('../experiments/succesful_patterns.txt', "r") as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]  # Strip newline characters
    
    lines.append(lprog)
    
    with open('../experiments/succesful_patterns.txt', "w") as file:
        for line in lines:
            file.write(line + "\n")  # Write each item followed by a newline    
    
    print(f"Experimental plan for '{lprog}' completed successfully.")

import json
import subprocess
import argparse
import os
import shutil
script_dir = os.path.dirname(os.path.abspath(__file__))
os.chdir(script_dir)
from utils.hgen_support import *

if __name__ == "__main__":
    main()
    