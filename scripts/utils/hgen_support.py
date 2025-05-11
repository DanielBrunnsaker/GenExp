#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:55:46 2024

@author: danbru
"""

import json
import re 
from datetime import datetime
import pandas as pd
from libchebipy._chebi_entity import ChebiEntity
import numpy as np

def create_experiment_folder(target, alpha, N, EXPERIMENTS_DIR): #genexp.py
    """
    Creates the required folder structure for the experiment.
    
    Args:
        target (str): The target molecule.
        alpha (float): The pattern ranking weight.
        N (int): The number of hypotheses.

    Returns:
        str: Path to the created experiment folder.
    """

    main_folder = EXPERIMENTS_DIR / f"{target}_{datetime.now().strftime('%Y%m%d%H%M')}"
    folder_structure = {
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
            "mass_spectrometry": {},
        },
        "results": {
            "metabolomics": {
                "processed": {},
                "raw": {},
                "tests": {
                    "models": {}}},
            "growth": {
                "processed": {},
                "raw": {},
                "tests": {
                    "models": {}}},
            "plots": {}
        },
        "versions": {}
    }

    def create_folders_recursively(base_path, structure):
        """ Recursively creates folder structure """
        for folder, sub_structure in structure.items():
            folder_path = base_path / folder
            folder_path.mkdir(parents=True, exist_ok=True)

            if isinstance(sub_structure, dict):
                create_folders_recursively(folder_path, sub_structure)

    create_folders_recursively(main_folder, folder_structure)

   # print(f"Experiment folder created at: {main_folder}")
    return str(main_folder)


def load_json(file_path): #genexp.py
    """
    Load a JSON file and return its content.
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return None
    
def get_entry_details(data, logic_program):
    # Loop through each entry in the 'allowed' list
    for entry in data['allowed']:
        # Check if the logic program matches
        if entry['logic_program'] == logic_program:
            # Return the number if a match is found
            return entry
    # Return None if no match is found
    return None

def extract_logic_program(text):
    
    # Updated regex pattern to capture logic program starting with "Cell(A):-exhibits_phenotype"
    pattern = r"Cell\(A\):-exhibits_phenotype.*?\."
    match = re.search(pattern, text, re.DOTALL)
    # If a match is found, return the logic program
    if match:
        return match.group(0)
    return None

def find_row_with_string(filename, target_string):
    with open(filename, 'r') as file:
        for line in file:
            if target_string in line:
                return line.strip()  # Return the line without leading/trailing whitespace
    return None  # If the string is not found in any line


def check_predicate_match(tuples_list, target_string):
    predicate_count = 0
    for item in tuples_list:
        if item[0] == target_string:
            predicate_count +=1
            
    return predicate_count

def flatten_tuple(tuple_input):
    flattened_list = []
    for item in tuple_input:
        if isinstance(item, tuple):
            flattened_list.extend(flatten_tuple(item))
        else:
            flattened_list.append(str(item))
    return flattened_list

def rank_strings_by_specificity(strings, relevance_order, relevance_scores):
    
    def calculate_specificity_score(string, relevance_scores):
        score = 0
        for idx, substring in enumerate(relevance_order):
            if substring in string:
                score += len(relevance_order) - idx  # Higher relevance substrings contribute more to the score
        return score

    rankings = [(string, calculate_specificity_score(string, relevance_scores)) for string in strings]
    rankings.sort(key=lambda x: x[1], reverse=True)  # Sort by specificity score in descending order
    return rankings

def calculate_specificity_score(string, relevance_scores):
        score = 0
        for substring, weight in relevance_scores.items():
            if substring in string:
                score += weight  # Add weight to score for each relevant substring found
        return score

def check_any_substring_in_target(df, column_name, target_string):
    
    # Iterate over the items in the specified column
    for item in df[column_name]:
        # Check if the current item is a substring of the target string
        if pd.notna(item) and str(item).lower() in target_string.lower():
            return True
    return False

def extract_inner_predicates(prolog_string):
        
    # Regular expression pattern to match the body of the feature
    body_pattern = r"feature\(\d+,\((.*)\)\)\."
    body_match = re.search(body_pattern, prolog_string)
    
    if body_match:
        body = body_match.group(1)
        # Regular expression pattern to match predicates within the body
        predicate_pattern = r"(\w+\([^\)]*\))"
        predicates = re.findall(predicate_pattern, body)
        return predicates
    return []

def find_identical_column(df, series):
    for column_name in df.columns:
        if df[column_name].equals(series):
            return column_name
    return None


def process_chebi_terms(text):
    
    # Find all terms starting with 'CHEBI:'
    chebi_terms = re.findall(r'CHEBI:\d+\w*', text)
    
    # Dictionary to hold original and processed terms
    term_mapping = {}
    
    # Apply the operation to each found term
    for term in chebi_terms:
        chebi_item = ChebiEntity(term)
        processed_term = chebi_item.get_name()
        term_mapping[term] = processed_term
    
    # Define a function to replace terms based on the dictionary
    def replace_term(match):
        return term_mapping[match.group(0)]
    
    # Replace each term in the original text with the processed term
    replaced_text = re.sub(r'CHEBI:\d+\w*', replace_term, text)
    return replaced_text


def normalize_column(column):
    min_val = min(column)
    max_val = max(column)
    return [(x - min_val) / (max_val - min_val) for x in column]


def calculate_scores(df, target_column, alpha):
    # Absolute values
    target_col_abs = df[target_column].abs()
    other_cols_abs_sum = df.drop(columns=[target_column]).abs().sum(axis=1) ** alpha
    
    # Logarithmic uniqueness score
    df['weighted_score'] = np.log(1 + target_col_abs) / np.log(1 + target_col_abs + other_cols_abs_sum)
    df_sorted = df.sort_values(by='weighted_score', ascending=False).drop('weighted_score', axis = 1)
    
    return df_sorted

def load_hypotheses(path, target):
         
     # Load results from the hypothesis generator
    hypotheses = pd.read_csv(f'{path}{target}_eCV_coefficients.csv', index_col = 0)
    hypotheses = hypotheses.mean(axis = 1)
    hypotheses = hypotheses/hypotheses.abs().max()
    hypotheses.columns = [target]
     
    return hypotheses,

def load_all_hypotheses(path, target, aaSet):
        
    st = pd.DataFrame(pd.read_csv(f'{path}/alanine_eCV_coefficients.csv', index_col = 0).mean(axis = 1), columns = ['alanine'])
    for trg in aaSet.columns[1:]:
        tmp = pd.DataFrame(pd.read_csv(f'{path}/{trg}_eCV_coefficients.csv', index_col = 0).mean(axis = 1), columns = [trg])
        st = st.merge(tmp, left_index = True, right_index = True)

    hypotheses = st[target]
    hypotheses = hypotheses/hypotheses.abs().max()
    hypotheses.columns = [target]
     
    return hypotheses, st/st.abs().max(axis = 0)


def calculate_string_score(string, word_scores):
    """
    Calculates the score of a string based on the provided word score dictionary.

    Parameters:
    string (str): The string to be scored.
    word_scores (dict): A dictionary where keys are words (substrings) and values are their respective scores.

    Returns:
    int: The total score of the string.
    """
    # Convert the string to lowercase for case-insensitive matching
    string = string.lower()
    
    # Initialize total score to 0
    total_score = 0
    
    # Iterate over each word in the dictionary
    for word, score in word_scores.items():
        # Count the occurrences of the word (substring) in the string
        count = string.count(word)
        
        # Add the score for each occurrence of the word
        total_score += count * score

    return total_score


