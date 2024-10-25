#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  9 14:31:39 2024

@author: danbru
"""

def rename_columns_with_ilp(df, gene):
    new_columns = [f'{gene}_ilp{i+1}' for i in range(len(df.columns))]
    df.columns = new_columns
    return df

def correlation_remover(dataset):
    '''
    Remove perfectly correlated features drom a dataframe
    
    Inputs:
        dataset: Dataframe containing your dataset
    
    Outputs:
        dataset: Dataframe with correlated features (except first occurence) removed
    
    '''
    print('Removing perfectly correlated features...')
    dataset = dataset.transpose().drop_duplicates(keep = 'first').transpose()
    return(dataset)


def define_aa(df, first_example):

    aaSet_reduced = df
    
    positives = list(aaSet_reduced.index)
    
    
    if first_example in positives:
        positives.remove(first_example)
    positives.insert(0, first_example)
    print('First example: '+positives[0])
    
    # Convert the list of entities to Datalog format
    positive_statements = [f"gene('{entity}')." for entity in positives]
    
    
    #file_path = '../prolog/metabolites_chi.f'
    file_path = '../prolog/patterns.f'
    
    with open(file_path, 'w') as file:
        file.write('\n'.join(positive_statements))
        
    return positives


def generate_frequent_features(target_folder):
    import subprocess
    print('Initializing Aleph for feature entailment...') 

    prolog_commands = [
        "[aleph_orig]",  # Run [aleph_orig]
        "read_all(patterns)",
        "induce_features",
        "saveQueries('generated_features/frequent_explanation.txt')",
        "show(features)",
        "stopQueriesSaving()",
        "saveQueries('generated_features/frequent_features.txt')",
        "show(train_pos)",
        "stopQueriesSaving()",
        "nl",
        "halt"
    ]
    
    # Now, using those, and the background file, generate features?
    # Make sure to set up basic KBs
    
    commands = ", ".join(prolog_commands)
    prolog_command = f"swipl -g '{commands}'"
    
    # Run the command using subprocess
    process = subprocess.Popen(prolog_command, cwd = target_folder, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, stdin=subprocess.PIPE)
    
    # Wait for the process to finish and get the output and errors
    stdout, stderr = process.communicate()  # Pass a newline to simulate pressing Enter

    if process.returncode == 0:
        print("Entailed - Prolog commands executed successfully")
    else:
        print(f"Error executing Prolog commands. Return code: {process.returncode}")


def create_datasets(pos_index, colname): 

    import pandas as pd
    
    pos = pd.read_csv('../prolog/generated_features/frequent_features.txt', sep = ' ', header = None).iloc[:,:-1]
    pos.index = pos_index
    pos = rename_columns_with_ilp(pos, colname)
    
    #dset_proc = correlation_remover(pos)
    dset_proc = pos
    #print(str(pos.shape[1]-dset_proc.shape[1])+' columns removed due to redundancy.')
    
    return dset_proc


def create_datasets_nored(pos_index, colname): 
    
    import pandas as pd

    pos = pd.read_csv('../prolog/generated_features/frequent_features.txt', sep = ' ', header = None).iloc[:,:-1]
    pos.index = pos_index
    pos = rename_columns_with_ilp(pos, colname)
    
    return pos

'''
def transform_data(X_train_complete, X_test_complete):
    
   
    preprocessing = Pipeline(
        steps=[
            ("variance_threshold", VarianceThreshold(threshold=0)),
            ("imputer", SimpleImputer(strategy = 'median')),
            ("std", StandardScaler())
        ]
    )
    
    prot_std = pd.DataFrame(
        data=preprocessing.fit_transform(X_train_complete),
        index=X_train_complete.index,
        columns=preprocessing[0].get_feature_names_out()
    )
    
    test_std = pd.DataFrame(
        data=preprocessing.transform(X_test_complete),
        index=X_test_complete.index,
        columns=preprocessing[0].get_feature_names_out()
    )

    return prot_std, test_std


def rank_rows_by_deviation(df):
    # Normalize the data using StandardScaler
    scaler = StandardScaler()
    df_normalized = pd.DataFrame(scaler.fit_transform(df), columns=df.columns)

    # Calculate the deviation of each row from the mean
    row_deviations = df_normalized.sub(df_normalized.mean(axis=0))

    # Calculate the absolute deviation for each row
    row_deviations_abs = row_deviations.abs()

    # Sum the absolute deviations for each row to get the total deviation
    row_deviations_total = row_deviations_abs.sum(axis=1)

    # Sort the DataFrame based on the total deviation (ascending order)
    df_ranked = df.iloc[row_deviations_total.argsort()[::-1]]

    return df_ranked
'''

def remove_or_update_duplicate_columns(df, duplicates_dict):
    # List to track columns that should be removed (duplicates only)
    columns_to_remove = []
    
    # Iterate through the columns in the DataFrame
    for i, col in enumerate(df.columns):
        # Skip columns that are already marked for removal
        if col in columns_to_remove:
            continue
        
        # Check if the current column is a duplicate of any column that is already a key in duplicates_dict
        is_duplicate = False
        for key in duplicates_dict:
            if key in df.columns and df[col].equals(df[key]):
                # If it matches an existing key, add this column to that key's list of duplicates
                duplicates_dict[key].append(col)
                columns_to_remove.append(col)  # Mark the duplicate for removal
                is_duplicate = True
                break
        
        # If the column is not a duplicate, ensure it is added to the dict with an empty list if it has no duplicates
        if not is_duplicate:
            duplicate_columns = [dup_col for dup_col in df.columns[i+1:] if df[col].equals(df[dup_col])]
            
            # Add the current column as a key in the dict, even if no duplicates are found
            if col not in duplicates_dict:
                duplicates_dict[col] = []
            
            # If duplicates are found, add them to the list of columns to remove
            if duplicate_columns:
                duplicates_dict[col].extend(duplicate_columns)
                columns_to_remove.extend(duplicate_columns)

    # Only remove columns that are true duplicates (do not remove the first occurrence)
    df = df.drop(columns=columns_to_remove)
    
    return df, duplicates_dict

def filter_lines(input_file, output_file):
    
    import re
    
    # Open the input file in read mode
    with open(input_file, 'r') as infile:
        # Read all lines into a list
        lines = infile.readlines()

    # List to store the feature numbers of lines that are kept
    kept_feature_numbers = []

    # Open the output file in write mode
    with open(output_file, 'w') as outfile:
        # Initialize a variable to track whether to skip the next line
        skip_next = False

        # Iterate through the lines with their index
        for i in range(len(lines)):
            # If skip_next is True, skip this iteration
            if skip_next:
                skip_next = False
                continue
            
            # Count the occurrences of the word "phenotype" in the current line (case insensitive)
            count = lines[i].lower().count('phenotype')

            # If "phenotype" appears more than once, set skip_next to True to skip the next line
            if count > 1:
                skip_next = True
                continue

            # Write the current line to the output file if it's not being skipped
            if i > 0 and lines[i-1].lower().count('phenotype') <= 1:
                outfile.write(lines[i-1])

                # Extract the feature number using regex and store it
                match = re.search(r'feature\((\d+),', lines[i-1])
                if match:
                    kept_feature_numbers.append(match.group(1))

            if i == len(lines) - 1 and count <= 1:
                outfile.write(lines[i])

                # Extract the feature number using regex and store it
                match = re.search(r'feature\((\d+),', lines[i])
                if match:
                    kept_feature_numbers.append(match.group(1))
    
    return kept_feature_numbers

def filter_lines(input_file, output_file):
    
    import re
    
    # Open the input file in read mode
    with open(input_file, 'r') as infile:
        # Read all lines into a list
        lines = infile.readlines()

    # List to store the feature numbers of lines that are kept
    kept_feature_numbers = []

    # Open the output file in write mode
    with open(output_file, 'w') as outfile:
        # Initialize a variable to track whether to skip the next line
        skip_next = False

        # Iterate through the lines with their index
        for i in range(len(lines)):
            # If skip_next is True, skip this iteration
            if skip_next:
                skip_next = False
                continue

            # Count the occurrences of the word "phenotype" in the current line (case insensitive)
            count = lines[i].lower().count('phenotype')

            # If "phenotype" appears more than once, set skip_next to True to skip the next line
            if count > 1:
                skip_next = True
                continue

            # Check for the "biological_target" without "interacts_with_metabolite" condition
            if 'drug_target' in lines[i] and 'interacts_with_metabolite' not in lines[i]:
                continue  # Skip this line

            # Write the current line to the output file if it's not being skipped
            if i > 0 and lines[i-1].lower().count('phenotype') <= 1:
                # Ensure it doesn't violate the biological_target condition
                if 'drug_target' not in lines[i-1] or 'interacts_with_metabolite' in lines[i-1]:
                    outfile.write(lines[i-1])

                    # Extract the feature number using regex and store it
                    match = re.search(r'feature\((\d+),', lines[i-1])
                    if match:
                        kept_feature_numbers.append(match.group(1))

            if i == len(lines) - 1 and count <= 1:
                # Ensure it doesn't violate the biological_target condition
                if 'drug_target' not in lines[i] or 'interacts_with_metabolite' in lines[i]:
                    outfile.write(lines[i])

                    # Extract the feature number using regex and store it
                    match = re.search(r'feature\((\d+),', lines[i])
                    if match:
                        kept_feature_numbers.append(match.group(1))
    
    return kept_feature_numbers


def filter_dataframe_columns(df, kept_feature_numbers):
    
    import re
    
    # Convert kept_feature_numbers to a set for faster lookups
    kept_feature_set = set(kept_feature_numbers)

    # Create a list to hold the names of columns to keep
    columns_to_keep = []

    # Iterate through the column names in the DataFrame
    for col in df.columns:
        # Use regex to extract the number after 'ilp' in the column name
        match = re.search(r'ilp(\d+)', col)
        if match:
            col_number = match.group(1)
            # If the number is in kept_feature_numbers, add the column to the list
            if col_number in kept_feature_set:
                columns_to_keep.append(col)

    # Filter the DataFrame to keep only the columns in columns_to_keep
    filtered_df = df[columns_to_keep]

    return filtered_df