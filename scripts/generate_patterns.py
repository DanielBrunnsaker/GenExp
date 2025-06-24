#!/usr/bin/env python3
# -*- coding: utf-8 -*-



#warnings.simplefilter("ignore", category=ConvergenceWarning)

def greedily_select_first_examples(df, n_FE):
    import pandas as pd
    
    # Combine into easily manageable phenotypic profiles
    df['Phenotypic_Profile'] = df['Gene > Phenotypes > Observable'] +'_'+ df['Gene > Phenotypes > Qualifier'] +'_'+ df['Gene > Phenotypes > Chemical'] +'_'+ df['Gene > Phenotypes > Condition']
    
    # Map genes to their phenotypic profiles
    gene_profiles = df.groupby('Gene > Systematic Name')['Phenotypic_Profile'].apply(set).to_dict()
    
    # Step 1: Initialize variables
    selected_genes = []
    covered_profiles = set()
    remaining_genes = set(gene_profiles.keys())
    
    # Step 2: Iteratively and greedily select genes with maximal phenotypic coverage
    while len(selected_genes) < n_FE and remaining_genes:
        max_new_profiles = -1
        best_gene = None
        
        for gene in sorted(remaining_genes):  
            new_profiles = gene_profiles[gene] - covered_profiles
            num_new_profiles = len(new_profiles)
            if num_new_profiles > max_new_profiles:
                max_new_profiles = num_new_profiles
                best_gene = gene
                
        if best_gene is None:
            break 
        
        selected_genes.append(best_gene)
        covered_profiles.update(gene_profiles[best_gene])
        remaining_genes.remove(best_gene)
    
        #print(best_gene, max_new_profiles)
    
    print("Selected Genes:", selected_genes)
    return selected_genes

def main():
    
    warnings.filterwarnings("ignore")
    os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses
    
    parser = argparse.ArgumentParser(description='Script with a command-line argument.')
    parser.add_argument('--N', type=str, help='Value for the "n_FE" variable.')
    args = parser.parse_args()
    if args.N is not None:
        n_FE = int(args.N)
        print(f'Variable "N" set to: {n_FE}')
    else:
        print('Variable "N" not provided.')
    
    print('new version')
    
    # set the terminal path
    target_folder = '../prolog'
    
    phenotypes = pd.read_csv('../data/phenotypes_manual_curation_v2.tsv', sep = '\t')
    
    # Remove some of the phenotypes, do this in the actual tsv before publication or similar
    
    subset = phenotypes[phenotypes['Gene > Phenotypes > Observable'] == 'chemical compound accumulation']
    terms_to_avoid = ['zinc cation', 'cadmium cation','calcium ion', 'cobalt cation','copper cation','gadolinium(3+)',
                      'iron cation','lead(2+)','magnesium cation','manganese cation','moleybdenum cation','nickel cation',
                      'phosphorus(1+)','potassium(1+)','proton','rubidium(1+)','sodium(1+)','sulfur(1+)']
    indices = []
    for term in terms_to_avoid:
        indices = indices + list(np.where(subset['Gene > Phenotypes > Chemical'] == term)[0])
    
    phenotypes = phenotypes.drop(indices)
    phenotypes = phenotypes[phenotypes['Gene > Phenotypes > Observable'] != 'cell size']
    phenotypes = phenotypes[phenotypes['Gene > Phenotypes > Observable'] != 'cell shape']
    
    #list_of_genes_to_try = list(phenotypes['Gene > Systematic Name'].value_counts().nlargest(n_FE).index)
    list_of_genes_to_try = greedily_select_first_examples(phenotypes, n_FE)
    # Example loop usage
    duplicates_dict = {}
    
    aap_set = pd.read_excel('../data/AA.xls', 
                            sheet_name='intracellular_concentration_mM', index_col = 0).iloc[:,1:]
    aap_set = aap_set.reset_index().groupby('ORF').mean()
    
    
    # Initialize full_frequent_dataset with no columns but the correct index
    full_frequent_dataset = pd.DataFrame(index=aap_set.index)
    
    for first_example in list_of_genes_to_try:  
        pos_index = define_aa(aap_set, first_example)
        generate_frequent_features(target_folder)
    
        xtrain_frequent = create_datasets(pos_index, first_example)
        
        input_file= '../prolog/generated_features/frequent_explanation.txt'
        output_file= '../prolog/generated_features/frequent_explanation_new.txt'
        
        # Removes double phenotypes
        kept_feature_numbers = filter_lines(input_file, output_file)
        xtrain = filter_dataframe_columns(xtrain_frequent, kept_feature_numbers)
        
        # Merge new data with the full dataset. Use pd.concat to accumulate all columns.
        merged_set = pd.merge(full_frequent_dataset, xtrain, left_index = True, right_index = True, how = 'left')
    
        # Remove/update duplicates
        merged_set, duplicates_dict = remove_or_update_duplicate_columns(merged_set, duplicates_dict)
            
        # Remove entries from the values that match the key
        for key in duplicates_dict:
            duplicates_dict[key] = [item for item in duplicates_dict[key] if item != key]
    
        # Update the full_frequent_dataset by concatenating it with the cleaned merged_set
        full_frequent_dataset = pd.concat([full_frequent_dataset, merged_set], axis=1)
        
        # Ensure no duplicate columns exist after concatenation
        full_frequent_dataset = full_frequent_dataset.loc[:,~full_frequent_dataset.columns.duplicated()]
        
        print(f"Total amount of unique features: {full_frequent_dataset.shape[1]}")
        
        with open('../prolog/generated_features/frequent_explanation_new.txt', 'r') as file:
            input_features = file.readlines()
        
        with open('../results/patterns/explanations/'+first_example+'.txt', 'w') as file:
            # Write each element of the list to a new line in the file
            for item in input_features:
                file.write(str(item) + '\n') 
                
                
        with open('../prolog/generated_features/frequent_explanation.txt', 'r') as file:
            complete_features = file.readlines()
        
        with open('../results/patterns/explanations/'+first_example+'_complete.txt', 'w') as file:
            # Write each element of the list to a new line in the file
            for item in complete_features:
                file.write(str(item) + '\n') 
                
        # Save the dict? 
        with open('../results/patterns/feature_dict.pickle', 'wb') as handle:
            pickle.dump(duplicates_dict, handle, protocol=pickle.HIGHEST_PROTOCOL)
    
        full_frequent_dataset.reset_index().to_feather('../results/patterns/datasets/frequent_20241114.feather')
    
    
    ## Rerun this?
    
    
    # Final verification step
    #assert len(duplicates_dict) == full_frequent_dataset.shape[1], "Mismatch between dict keys and DataFrame columns!"
    
import os
os.chdir(os.environ["GEN_EXP_ROOT_DIR"] + '/scripts')
import pandas as pd
from sklearn.feature_selection import VarianceThreshold
import warnings
import pickle
import re
import numpy as np
from pattern_generation_utils import *
import pandas as pd
import argparse

main()