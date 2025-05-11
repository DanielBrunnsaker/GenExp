#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  5 11:54:45 2025

@author: danbru
"""

#from prefect import task, flow
import shutil
from pathlib import Path
import argparse
import os
import config

# Import functions from utility scripts
from pattern_selection import pattern_selection
from hypothesis_generation import hypothesis_generation
from autoformalize_protocol import autoformalize
from plate_layout import generate_layout
from hamilton_protocol import hamilton_protocol
from ms_runlist_generation import generate_mass_spec_metadata
from utils.hgen_support import create_experiment_folder, load_json


# Set a base directory for all relative paths used in the flow
BASE_DIR = Path(__file__).resolve().parent
EXPERIMENTS_DIR = BASE_DIR / "../experiments"
EXPERIMENTS_DIR = EXPERIMENTS_DIR.resolve()

#@task
def create_folder_structure(target, alpha, N, EXPERIMENTS_DIR):
    """ Script creating the folder structure of the experiment-files """
    return create_experiment_folder(target, alpha, N, EXPERIMENTS_DIR)

#@task
def print_header():
    ascii_art = r"""
               o 
               |
           .-------.       
          /         \      
        <|  ^     ^  |>  Beep boop! Let's do some science!
         |   [___]   |     _
         |___________|    | |
     \----/         \----/|-|
          |   EVE   |    /   \
          | _______ |   (_____) 
          ||       ||       
          ||       ||       
         (__)     (__)      
    """

    print(ascii_art)

#@task
def run_pattern_selection(target, alpha, N):
    """Run pattern selection step as a function call."""
    pattern_selection(target, alpha, N)

#@task
def run_hypothesis_generation(output_folder, T, N, variations):
    """Run hypothesis generation step as a function call."""
    hypothesis_generation(output_folder, T, N, variations)

#@task
def run_autoformalize(output_folder):
    """Autoformalize experimental plan as a function call."""
    autoformalize(output_folder)

#@task
def user_decision():
    """Prompt the user for action."""
    return input("Enter action (1 (continue), 2 (forbid pattern), 3 (cancel): ").strip()

#@task
def clean_up(output_folder):
    """Clean up temporary files."""
    shutil.rmtree(output_folder)

#@task
def finalize_layout(output_folder):
    """Generate the final plate configuration."""
    generate_layout(output_folder)

#@task
def design_dispensing_layout(output_folder):
    """Generate the final plate configuration."""
    hamilton_protocol(output_folder)

#@task
def prepare_mass_spectrometry_metadata(output_folder):
    """Generate the final plate configuration."""
    generate_mass_spec_metadata(output_folder)

#@task
def save_settings(output_folder):
    # Copy the current config file so that it is
    # clear what settings were used to generate the experiment
    shutil.copy2(BASE_DIR / 'config.py',
                 os.path.join(output_folder, 'versions/frozen_config.py'))
    shutil.copy2(BASE_DIR / '../context/hypgen_per_pattern_temp.txt',
                 os.path.join(output_folder, 'versions/hypgen_per_pattern_usedcontext.txt'))
    shutil.copy2(BASE_DIR / '../context/hypothesis_selection.txt',
                 os.path.join(output_folder, 'versions/hypothesis_selection_usedcontext.txt'))
    shutil.copy2(BASE_DIR / '../context/refine_plan.txt',
                 os.path.join(output_folder, 'versions/refine_plan_usedcontext.txt'))
    shutil.copy2(BASE_DIR / '../context/expdesign_context.txt',
                 os.path.join(output_folder, 'versions/expdesign_context_usedcontext.txt'))

#@task
def add_logic_program_to_list(list_of_programs, output_folder, target):

    # Add used logic program to the performed experiments list

    with open(Path(output_folder) / f'../{list_of_programs}.txt', "r") as file:
        lines = file.readlines()
        lines = [line.strip() for line in lines]  # Strip newline characters

    details = load_json(Path(output_folder) / 'protocol/protocol.json')
    lines.append(f"{target}:"+details['logic_program'])

    with open(Path(output_folder) / f'../{list_of_programs}.txt', "w") as file:
        for line in lines:
            file.write(line + "\n")  # Write each item followed by a newline


#@flow
def experiment_pipeline(target, alpha, N):
    """Main Prefect workflow for running the experiment pipeline."""

    print_header()

    # Create experiment folder first and some metadata regarding the experiment
    print("\nCreating experiment folder...")
    output_folder = create_folder_structure(target, alpha, N, EXPERIMENTS_DIR)
    print(f'✅ Folder created: {output_folder} \n')

    # Run a pattern selection process on the logic programs produced by aleph for given target
    print(f"Running pattern selection for {target}, alpha={alpha}, N={N}...")
    pattern_selection(target, alpha, N, output_folder)
    print("✅ Pattern selection complete! \n")

    # Run hypothesis generation using GPT4o (or other) alongside the provided logic programs
    print("Generating hypotheses...")
    hypothesis_generation(output_folder, N)
    print("✅ Hypothesis generation complete! \n")

    # Autoformalize experimental plan into a valid JSON-file
    # (note that it might have to retry if the JSON-file does not pass the validity check)
    print("Autoformalizing experimental plan...")
    autoformalize(output_folder)
    print("✅ Autoformalization complete! \n ")

    # Print the selected hypothesis in case one wants to do a sanity-check of
    # the produced experimental plan before proceeding
    if config.PRINT_HYPOTHESIS_DURING_WORKFLOW:
        print("Selected hypothesis: \n")
        selected_hypothesis = open(Path(output_folder) / 'hypothesis/selected_hypothesis/hypothesis.txt')
        print(selected_hypothesis.read())

    # Prompt for decision regarding how to proceed, if human-in-the-loop
    if config.ALLOW_FOR_CHOICE:

        print("\n Prompting for user decision...")
        action = user_decision()

        if action == "1":
            print("Proceeding with experimental design...")
        elif action == "2":
            # this option will remove the logic program from the pool of available programs
            print("Forbidding pattern. Cleaning up and exiting.")
            add_logic_program_to_list('forbidden_patterns', output_folder, target)
            clean_up(output_folder)
            return
        elif action == "3":
            print("Exiting the script.")
            clean_up(output_folder)
            return

    # Generate experimental layout and dispensing scheme using PLAID and
    # metadata derived from the autoformalized experimental protocol
    print("Generating plate layout...")
    finalize_layout(output_folder)
    design_dispensing_layout(output_folder)
    print("✅ Experimental plan finalized! \n")

    # Generate a randomized order mass-spec runlist and rapidfire settings
    # (e.g. cartridge, polarity, ...)
    if config.PERFORM_MASS_SPEC:
        print("Generating metadata for mass-spectrometry analysis. ")
        prepare_mass_spectrometry_metadata(output_folder)
        print("✅ Metadata successfully prepared! \n")


    # Save the settings for that run (e.g. llm-version, prompts/context, temperatures, and so forth)
    save_settings(output_folder)

    # Decide whether to classify a logic program/target-combination as completed
    # (thereby removing it from the pool of available programs)
    if not config.REUSE_LOGIC_PROGRAM:
        add_logic_program_to_list('succesful_patterns', output_folder, target)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run experiment pipeline using Prefect.")

    parser.add_argument('--target', required=True,
                        help="Target metabolite")
    parser.add_argument('--alpha', required=True, type=float,
                        help="Parameter that weighs the coefficients of the hypotheses based on its occurence for other targets")
    parser.add_argument('--N', required=True, type=int,
                        help="Number of logic programs to run through the selection process")
    args = parser.parse_args()

    experiment_pipeline(args.target, args.alpha, args.N)
