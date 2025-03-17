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
import numpy as np
from scipy.stats import shapiro
from statsmodels.stats.diagnostic import het_breuschpagan

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
    growth_data = pd.read_csv(EXPERIMENT_DIR / 'results/growth/processed/growth_parameters.tsv',
                              sep='\t', index_col=0)
    
    # Create design table from dispensing layout.
    layout = pd.read_excel(EXPERIMENT_DIR / 'protocol/hamilton/pipetting_layout.xlsx')
    layout.rename(columns={'well': 'Well'}, inplace=True)
    design_table = create_design_table(layout)
    data = design_table.merge(growth_data.drop('Summary', axis=1), left_on='Well', right_index=True)
    data = data[data['Summary'] != 'Media control']
    
    # Loop over each metric
    for testing_column in ['AUC', 'mu', 'MaxOD']:
        categorized_data = set_datatypes(data, testing_column)
        
        if testing_column == 'AUC':
            # For AUC, use a Gamma GLM with log link. Did some empirical testing here, and this one seemed the best.
            model = smf.glm(
                formula=f'{testing_column} ~ C(Treatment) * C(Supplement)',
                data=categorized_data,
                family=sm.families.Gamma(link=sm.families.links.Log())
            ).fit()
            model_type = "Gamma GLM (log link)"
        
        elif testing_column == 'mu':
            # For mu, use the Inverse Gaussian GLM with a log link. Did some empirical testing here, and this one seemed the best.
            model = smf.glm(
                formula='mu ~ C(Treatment) * C(Supplement)',
                data=categorized_data,
                family=sm.families.InverseGaussian(link=sm.families.links.Log())
            ).fit()
            model_type = "Inverse Gaussian GLM (log link)"
        
        elif testing_column == 'MaxOD':
            # For MaxOD, apply a log transformation and use OLS. Did some empirical testing here, and this one seemed the best.
            log_column = f'log_{testing_column}'
            categorized_data[log_column] = np.log(categorized_data[testing_column])
            model = smf.ols(
                formula=f'{log_column} ~ C(Treatment) * C(Supplement)',
                data=categorized_data
            ).fit()
            model_type = "OLS on log-transformed data"
        
        # Print final model summary for the current metric.
        print(f'\n########### Final Testing for {testing_column} using {model_type} ###########\n')
        print(model.summary())
        
        # Choose residuals for diagnostic tests.
        if testing_column in ['AUC', 'mu']:
            resid = model.resid_response  # For GLMs, using the response residuals.
            bp_resid = model.resid_pearson
        else:
            resid = model.resid  # For OLS.
            bp_resid = resid
        
        # Normality test (Shapiro-Wilk)
        shapiro_stat, shapiro_p = shapiro(resid)
        print("\nNormality Test (Shapiro-Wilk):")
        print(f"Statistic: {shapiro_stat:.4f}, p-value: {shapiro_p:.4f}")
        
        # Heteroscedasticity test (Breusch-Pagan)
        bp_stat, bp_pvalue, fvalue, f_pvalue = het_breuschpagan(bp_resid, model.model.exog)
        print("\nHeteroscedasticity Test (Breusch-Pagan):")
        print(f"LM Statistic: {bp_stat:.4f}, LM p-value: {bp_pvalue:.4f}")
        print(f"F-Statistic: {fvalue:.4f}, F p-value: {f_pvalue:.4f}")
        
        # Save model coefficients and the fitted model.
        results_df = pd.DataFrame({
            "Coefficient": model.params,
            "Std Error": model.bse,
            "t-value": model.tvalues,
            "p-value": model.pvalues
        })
        results_df.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}.tsv', sep='\t')
        model.save(EXPERIMENT_DIR / f'results/growth/tests/models/{testing_column}.pickle')


'''
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import statsmodels.formula.api as smf
import statsmodels.api as sm
from scipy.stats import shapiro, levene
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import numpy as np

def two_way_anova_testing(EXPERIMENT_DIR):
    """
    This function performs a two-way ANOVA using OLS,
    tests assumptions (normality of residuals and homogeneity of variances),
    and conducts a Tukey HSD post hoc test for pairwise comparisons.
    """
    # Load observables
    growth_data = pd.read_csv(EXPERIMENT_DIR / 'results/growth/processed/growth_parameters.tsv', sep='\t', index_col=0)
    
    # Create design table from the pipetting layout
    layout = pd.read_excel(EXPERIMENT_DIR / 'protocol/hamilton/pipetting_layout.xlsx')
    layout.rename(columns={'well':'Well'}, inplace=True)
    design_table = create_design_table(layout)
    
    # Merge design table with growth data
    data = design_table.merge(growth_data.drop('Summary', axis=1), left_on='Well', right_index=True)
    data = data[data['Summary'] != 'Media control']
    
    for testing_column in ['AUC', 'mu', 'MaxOD']:
        categorized_data = set_datatypes(data, testing_column)
        
        #categorized_data[testing_column] = np.log(categorized_data[testing_column])
        
        # Fit an OLS model for two-way ANOVA
        model_ols = smf.ols(formula=f'{testing_column} ~ C(Treatment) * C(Supplement)', data=categorized_data).fit()
        #model_ols = smf.ols(formula=f'{testing_column} ~ C(Control_Type) + C(Treatment) * C(Supplement) + C(Treatment):C(Control_Type)', data=categorized_data).fit()
        #formula=f'{testing_column} ~ C(Control_Type) + C(Treatment) * C(Supplement)', data=rls_data
        anova_results = sm.stats.anova_lm(model_ols, typ=2)  # Type II ANOVA
        
        print(f'\n ########### Two-Way ANOVA for {testing_column} ########### \n')
        print(anova_results)
        
        # Assumption testing: Normality of residuals using Shapiro-Wilk test
        shapiro_stat, shapiro_p = shapiro(model_ols.resid)
        print("Shapiro-Wilk test for normality of residuals:")
        print(f"Statistic = {shapiro_stat:.4f}, p-value = {shapiro_p:.4f}")
        
        # Assumption testing: Homogeneity of variances using Levene's test
        # Group by the combination of Treatment and Supplement
        groups = [model_ols.resid[group.index] 
                  for name, group in categorized_data.groupby(["Treatment", "Supplement"])]
        levene_stat, levene_p = levene(*groups)
        print("Levene's test for homogeneity of variances:")
        print(f"Statistic = {levene_stat:.4f}, p-value = {levene_p:.4f}")
        
        # Post hoc test: Tukey HSD for all pairwise comparisons
        # Create a combined group variable
        categorized_data["Group"] = categorized_data["Treatment"].astype(str) + "_" + categorized_data["Supplement"].astype(str)
        tukey_results = pairwise_tukeyhsd(
            endog=categorized_data[testing_column],
            groups=categorized_data["Group"],
            alpha=0.05
        )
        print("\nTukey HSD post hoc test results:")
        print(tukey_results.summary())
        
        # Save ANOVA results, assumption test outputs, and Tukey HSD summary
        anova_results.to_csv(EXPERIMENT_DIR / f'results/growth/tests/anova_{testing_column}.tsv', sep='\t')
        with open(EXPERIMENT_DIR / f'results/growth/tests/anova_assumptions_{testing_column}.txt', 'w') as f:
            f.write("Shapiro-Wilk test for normality of residuals:\n")
            f.write(f"Statistic = {shapiro_stat:.4f}, p-value = {shapiro_p:.4f}\n\n")
            f.write("Levene's test for homogeneity of variances:\n")
            f.write(f"Statistic = {levene_stat:.4f}, p-value = {levene_p:.4f}\n")
        with open(EXPERIMENT_DIR / f'results/growth/tests/tukey_{testing_column}.txt', 'w') as f:
            f.write(str(tukey_results.summary()))

'''

