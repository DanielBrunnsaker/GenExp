#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:41:26 2025

@author: danbru
"""
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'


import statsmodels.formula.api as smf
import statsmodels.api as sm

def create_design_table(df):
    design_df = df.copy()
    
    # Process Supplement Column
    unique_supplement_values = sorted(df['Supplement (uL)'].unique())
    supplement_mapping = {0: 'None', unique_supplement_values[1]: 'PosLow', unique_supplement_values[2]: 'PosHigh'}
    design_df['Supplement'] = df['Supplement (uL)'].map(supplement_mapping)
    
    # Process Negative Control Column
    max_negative = df['Negative Control (uL)'].max()
    design_df.loc[df['Negative Control (uL)'] == max_negative, 'Supplement'] = 'NegHigh'
    
    # Process Treatment Column
    treatment_mapping = {0: 'None', df['Treatment (uL)'].max(): 'Yes'}
    design_df['Treatment'] = df['Treatment (uL)'].map(treatment_mapping)
    
    # Keep only relevant columns
    design_df = design_df[['Well', 'Summary', 'Treatment', 'Supplement']]
    
    return design_df


def set_datatypes(data, testing_column):
    
    rls_data = data[[testing_column, 'Treatment', 'Supplement']]
    
    # Ensure variables are categorical?
    rls_data['Treatment'] = pd.Categorical(
        rls_data['Treatment'],
        categories=['None', 'Yes'],     
        ordered=True
    )
    rls_data['Supplement'] = pd.Categorical(
        rls_data['Supplement'],
        categories=['None', 'PosLow', 'PosHigh', 'NegHigh'], 
        ordered=True
    )
    
    return rls_data

def growth_testing(EXPERIMENT_DIR):

    # Load observables
    growth_data = pd.read_csv(EXPERIMENT_DIR / 'results/growth/processed/growth_parameters.tsv', sep = '\t', index_col = 0)
    
    # Create design table, for use with testing
    # Design table is created from the dispensing layout, seemed to be the most robust way of getting it right
    layout = pd.read_excel(EXPERIMENT_DIR / 'protocol/hamilton/pipetting_layout.xlsx')
    layout.rename(columns={'well':'Well'}, inplace=True)
    design_table = create_design_table(layout)
    data = design_table.merge(growth_data.drop('Summary', axis = 1), left_on = 'Well', right_index = True)
    data = data[data['Summary'] != 'Media control']
    
    
    for testing_column in ['AUC', 'mu', 'MaxOD']: # Whatever metric we want to test
    
        categorized_data = set_datatypes(data, testing_column)
        
        # Might have to make this conditional, have not tested if mu and max OD follows a gamma dist. Probably not.
        model_gamma = smf.glm(
            formula=f'{testing_column} ~ C(Treatment) * C(Supplement)',
            data=categorized_data,
            family=sm.families.Gamma(link=sm.families.links.Log())
        ).fit()
        
        print(f'\n ########### testing for {testing_column} ########### \n')
        print(model_gamma.summary()) # THis output should be enough i think
        
        # Extract model results
        results_df = pd.DataFrame({
            "Coefficient": model_gamma.params,
            "Std Error": model_gamma.bse,
            "z-value": model_gamma.tvalues,
            "p-value": model_gamma.pvalues
        })
        
        results_df.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}.tsv', sep = '\t')
        model_gamma.save(EXPERIMENT_DIR / f'results/growth/tests/models/{testing_column}.pickle')

