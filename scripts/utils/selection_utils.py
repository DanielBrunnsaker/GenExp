#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 13:13:37 2025

@author: danbru
"""
import pandas as pd
import json
from pathlib import Path

from utils.hgen_support import load_all_hypotheses, calculate_scores, find_row_with_string, calculate_string_score
import config

def add_negative_examples(input_path, output_path, neg_examples):

    # temporarily modify the context
    with open(input_path, "r") as file:
        content = file.read()

    # Replace <PLACEHOLDER> with the joined items
    replacement_text = ", ".join(neg_examples)
    updated_content = content.replace("<PLACEHOLDER>", replacement_text)

    # Write the updated content back to the file
    with open(output_path, "w") as file:
        file.write(updated_content)

def load_supporting_data(base_dir):
    """
    Load and return supporting data dictionaries.
    """
    # id_dict
    id_dict = pd.read_csv(base_dir / '../data/uniprot_id_map.tsv', sep='\t')
    id_dict['To'] = id_dict['To'].str.replace('sce:', '')
    id_dict = id_dict.set_index('From')['To'].to_dict()

    # orf_dict
    orf_dict = pd.read_csv(base_dir / '../data/orf_dict.tsv', sep='\t', header=None)[[1,3]]
    orf_dict = orf_dict.set_index(1)[3].to_dict()

    return id_dict, orf_dict

def load_and_filter_hypotheses(coefficients_path, target, aaSet, base_dir):
    """
    Load all hypotheses, then filter using metrics.
    Return the filtered hypotheses and raw hypotheses data.
    """
    # Load all hypotheses
    hypotheses, all_hypotheses = load_all_hypotheses(coefficients_path, target, aaSet)

    # Filter by R2 > 0
    metrics = pd.read_csv(base_dir / '../results/metrics/eCV_results_20241114.csv')
    mean_metrics = metrics.groupby('Amino acid').mean()
    filtered_hypotheses = all_hypotheses[mean_metrics.index[mean_metrics['R2'] > 0]]

    return hypotheses, filtered_hypotheses

def rank_hypotheses(all_hypotheses, target, alpha):
    """
    Rank and return the sorted specificity.
    """
    sorted_specificity = calculate_scores(all_hypotheses, target, alpha=alpha)
    # Filter out non-zero values for the target
    filtered_sorted_specificity = sorted_specificity[target][sorted_specificity[target] != 0]
    return filtered_sorted_specificity

def pick_best_clause(hypothesis, duplicate_dict, base_dir):
    """
    Determine the best clause among duplicates based on explanation scores.
    Return the clause identifier (e.g., 'X_12p34') to be used.
    """
    # If no duplicates, return as-is
    if not duplicate_dict[hypothesis]:
        return hypothesis

    best_clause = hypothesis
    best_score = -1
    for item in duplicate_dict[hypothesis]:
        # Extract explanation from disk
        try:
            temp_expl = find_row_with_string(
                f"{base_dir}/../results/patterns/explanations/{item.split('_')[0]}.txt",
                f"({item.split('_')[1].split('p')[1]},"
            )
            score = calculate_string_score(temp_expl, config.RELEVANCE_SCORES)
            if score > best_score:
                best_score = score
                best_clause = item
        except FileNotFoundError:
            # If explanation file doesn't exist, skip
            continue

    return best_clause

def prepare_clause_text(current_clause, base_dir):
    """
    Return the readable portion of the clause by extracting
    from the explanation text file and doing replacements.
    """
    # Extract file name for explanation
    first_example = current_clause.split('_')[0]
    # Explanation row
    h_explanation = find_row_with_string(
        f'{base_dir}/../results/patterns/explanations/{first_example}.txt',
        f"({current_clause.split('_')[1].split('p')[1]},"
    )
    # Clause
    clause = h_explanation.split(',(')[1][:-3].replace('gene(A):-', 'Cell(A):-')
    return clause

def save_hypothesis_details(
    logic_program, 
    target, 
    directionality, 
    coefficient, 
    initial_prompt, 
    output_folder, 
    hypothesis_counter
):
    """
    Save hypothesis metadata as a JSON file.
    """
    details = {
        "logic_program": logic_program.replace("'", ""),
        "observable": target,
        "qualifier": directionality,
        "coefficient": coefficient,
        "initial_prompt": initial_prompt
    }

    output_path = Path(output_folder).resolve() / 'hypothesis/generated_hypotheses/initial_stage'
    output_path.mkdir(parents=True, exist_ok=True)

    with open(output_path / f'hypothesis_details_{hypothesis_counter}.json', 'w') as file:
        json.dump(details, file, indent=4)