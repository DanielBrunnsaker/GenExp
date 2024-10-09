#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 19 15:55:46 2024

@author: danbru
"""


def save_string_to_json(path_dir, content, prefix, target, N, T, alpha, beta):
    
    import os
    import json
    from datetime import datetime
    
    """
    
    parser.add_argument('--target', type=str, help='Value for the "target" variable.')
    parser.add_argument('--N', type=str, help='Value for the "N" variable.')
    parser.add_argument('--T', type=str, help='Value for the "T" variable.')
    parser.add_argument('--alpha', type=str, help='Value for the "alpha" variable.')
    parser.add_argument('--beta', type=str, help='Value for the "beta" variable.')
    
    Saves the content as JSON to a file in directory B. The filename is generated using
    target, N, T, alpha, beta, and a timestamp.
    
    Args:
    path_dir (str): The directory where the file will be saved.
    content (dict): LLM output or data to be saved in JSON format.
    
    prefix (str): descriptor of the file.
    target (str): Which amino acid is the target of the hypothesis.
    N (int): N hypotheses given to the initial prompt.
    T (float): temperature of the hypothesis generation step.
    alpha (str): alpha-coefficient.
    beta (str): beta-coefficient.
    
    Returns:
    str: The path to the saved file.
    """
    # Get the current timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    
    # Construct the filename
    filename = f"{prefix}_{target}_{N}_{T}_{alpha}_{beta}_{timestamp}.json"
    
    # Combine the path and the filename
    file_path = os.path.join(path_dir, filename)
    
    # Write the content to the JSON file
    with open(file_path, 'w') as file:
        json.dump(content, file, indent=4)
    
    print(f"File saved to: {file_path}")
    return file_path

def save_string_to_file(path_dir, content, prefix, target, N, T, alpha, beta):
    
    import os
    from datetime import datetime
    
    """
    
    parser.add_argument('--target', type=str, help='Value for the "target" variable.')
    parser.add_argument('--N', type=str, help='Value for the "N" variable.')
    parser.add_argument('--T', type=str, help='Value for the "T" variable.')
    parser.add_argument('--alpha', type=str, help='Value for the "alpha" variable.')
    parser.add_argument('--beta', type=str, help='Value for the "beta" variable.')
    
    Saves the string A to a text file in directory B. The filename is generated using
    Variable1, Variable2, Variable3, and a timestamp.
    
    Args:
    path_dir (str): The directory where the file will be saved.
    
    content (str): LLM output for that specific step
    
    prefix (str): descriptor of the file
    target (str): Which amino acid is the target of the hypothesis.
    N (int): N hypotheses given to the initial prompt.
    T (float): temperature of the hypothesis generation step
    alpha (str): alpha-coefficient.
    beta (str): beta-coefficient.
    
    Returns:
    str: The path to the saved file.
    """
    # Get the current timestamp
    timestamp = datetime.now().strftime("%Y%m%d_%H%M")
    
    # Construct the filename
    filename = f"{prefix}_{target}_{N}_{T}_{alpha}_{beta}_{timestamp}.txt"
    
    # Combine the path and the filename
    file_path = os.path.join(path_dir, filename)
    
    # Write the content to the file
    with open(file_path, 'w') as file:
        file.write(content)
    
    print(f"File saved to: {file_path}")
    return file_path

        
def combinations_start_with_first(lst):
    
    from itertools import combinations

    first_item = lst[0]
    result = []
    
    for r in range(1, len(lst) + 1):  # lengths from 1 to N
        for comb in combinations(lst[1:], r - 1):  # generate combinations for the rest of the list
            result.append((first_item,) + comb)  # prepend the first item

    return result   

def find_row_with_string(filename, target_string):
    with open(filename, 'r') as file:
        for line in file:
            if target_string in line:
                return line.strip()  # Return the line without leading/trailing whitespace
    return None  # If the string is not found in any line

def split_program_into_atoms(program):
    
    import re
    
    # Regular expression pattern to match atoms
    atom_pattern = r"(phenotype_chemical|phenotype|chebi_name|chebi|condition)\(([^)]+)\)"

    # Find all matches of atoms in the program
    atoms = re.findall(atom_pattern, program)

    # Format and return the atoms
    formatted_atoms = [(predicate, arguments.split(',')) for predicate, arguments in atoms]
    return formatted_atoms

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
    
    import pandas as pd
    
    # Iterate over the items in the specified column
    for item in df[column_name]:
        # Check if the current item is a substring of the target string
        if pd.notna(item) and str(item).lower() in target_string.lower():
            return True
    return False

def extract_inner_predicates(prolog_string):
    
    import re
    
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
    
    import re
    from libchebipy._chebi_entity import ChebiEntity

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

# Function to calculate scores for each row based on target column and uniqueness
def calculate_scores(df, target_column, alpha=1, beta=1):
    # Target values: maximize absolute values in the target column
    target_values = df[target_column].abs()
    
    # Non-uniqueness penalty: sum of absolute values across other columns
    penalty_values = df.drop(columns=[target_column]).abs().sum(axis=1)
    
    # Calculate the score: prioritize high target values, penalize non-unique rows
    scores = alpha * target_values - beta * penalty_values
    
    # Add the scores to the DataFrame
    df['score'] = scores
    
    # Sort by the score column in descending order
    df_sorted = df.sort_values(by='score', ascending=False).drop(columns=['score'])
    
    return df_sorted


def load_hypotheses(path, target):
     
    import pandas as pd
    
     # Load results from the hypothesis generator
    hypotheses = pd.read_csv(f'{path}{target}_eCV_coefficients.csv', index_col = 0)
    hypotheses = hypotheses.mean(axis = 1)
    hypotheses = hypotheses/hypotheses.abs().max()
    hypotheses.columns = [target]
     
    return hypotheses,

def load_all_hypotheses(path, target, aaSet):
    
    import pandas as pd
    
    st = pd.DataFrame(pd.read_csv(f'{path}alanine_eCV_coefficients.csv', index_col = 0).mean(axis = 1), columns = ['alanine'])
    for trg in aaSet.columns[1:]:
        tmp = pd.DataFrame(pd.read_csv(f'{path}{trg}_eCV_coefficients.csv', index_col = 0).mean(axis = 1), columns = [trg])
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

def user_prompt():
    while True:
        # Present the user with options
        print("Choose an option:")
        print("[1] Go ahead with experiment")
        print("[2] Unsafe experiment, reprompt")
        print("[3] Cancel")
        
        # Get user input
        choice = input("Enter your choice (1, 2, or 3): ")

        # Handle the different choices
        if choice == "1":
            print("Proceeding with the script...")
            return "go_ahead"
        elif choice == "2":
            print("Retrying...")
            return "retry"
        elif choice == "3":
            print("Cancelling the script...")
            return "cancel"
        else:
            print("Invalid choice, please try again.")
            
            


