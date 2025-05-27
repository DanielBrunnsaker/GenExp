#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 15:41:26 2025

@author: danbru
"""
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import statsmodels.formula.api as smf
import numpy as np
from tqdm import tqdm
import json
from sklearn.utils import resample


def create_design_table(df):
    """
    Create a design-table from the dispensing layouts.

    Parameters
    ----------
    df : dataframe
        dispensing layout.

    Returns
    -------
    design_df : dataframe
        a simple design matrix.

    """
    
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

def create_dose_table(df):
    """
    Takes in the dispensing layout and creates a table for use with testing.
    Includes the dose for the supplements, and flags for treatments and
    negative controls.

    Parameters
    ----------
    df : dataframe
        Dispensing layout generated for the hamilton protocols.

    Returns
    -------
    d: dataframe
        dataframe consisting of all relevant flags for testing.

    """
    
    d = df.copy()
    
    # 1) Rename raw columns to better names for this
    d = d.rename(columns={
        'media_supplementation':            'SupplementName',
        'media_supplementation_doses':      'SupplementDose_mM',
        'treatment':                        'TreatmentName',
        'treatment_parameters':             'TreatmentDose_uL',
        'Negative Control (uL)':            'NegativeControlVol_uL'
    })

    # 2) Parse numeric doses
    def parse_numeric(x):
        try:
            s = str(x)
            m = pd.Series([s]).str.extract(r'(\d+(?:\.\d+)?)', expand=False)[0]
            return float(m) if pd.notna(m) else 0.0
        except Exception:
            return 0.0

    d['SupplementDose_mM']   = d['SupplementDose_mM'].apply(parse_numeric)
    d['TreatmentDose_uL']    = d['TreatmentDose_uL'].apply(parse_numeric)
    d['NegativeControlVol_uL']= d['NegativeControlVol_uL'].fillna(0).astype(float)

    # 3) Flag negative control by extracting the dispensed volumes
    d['NegControlFlag'] = (d['NegativeControlVol_uL'] > 0).astype(int)
    d['Doses'] = d['SupplementDose_mM']
    
    # 4) Zero out supplement dose when neg control
    d.loc[d['NegControlFlag'] == 1, 'SupplementDose_mM'] = 0.0

    # 5) Clean TreatmentName and cast categorical
    d['TreatmentName'] = d['TreatmentName'].where(d['TreatmentDose_uL'] > 0, 'None').fillna('None')
    cats = ['None'] + [t for t in d['TreatmentName'].unique() if t != 'None']
    d['TreatmentName'] = pd.Categorical(d['TreatmentName'], categories=cats, ordered=True)

    # 6) Cast SupplementName categorical including 'None'
    d['SupplementName'] = d['SupplementName'].fillna('None').astype(str)
    supp_cats = ['None'] + [s for s in d['SupplementName'].unique() if s != 'None']
    d['SupplementName'] = pd.Categorical(d['SupplementName'], categories=supp_cats)

    # 7) Return cleaned table
    return d[['Well','Summary',
              'TreatmentName','TreatmentDose_uL',
              'SupplementName','SupplementDose_mM',
              'NegativeControlVol_uL','NegControlFlag','Doses']]


def bootstrap_summary(df, formula, B=1000):
    """
    Given a dataframe and a glm-formula, fit an OLS model, bootstrap it
    and compute relevant metrics.

    Parameters
    ----------
    df : dataframe
        dataframe with relevant experimental variables.
    formula : string
        the glm-formula.
    B : integer, optional
        the numbero of resamples. The default is 1000.

    Returns
    -------
    summary : df
        summary statistics.
    boot_df : df
        bootstrap estimates.

    """
    
    np.random.seed(0)
    
    # 1) Fit once to get param names, because we are lazy
    orig_mod    = smf.ols(formula, data=df).fit()
    params      = orig_mod.params
    
    # 2) Bootstrap
    boot_mat = np.zeros((B, len(params)))
    for i in tqdm(range(B), desc="Bootstrapping"):
        samp    = resample(df)
        mod_b   = smf.ols(formula, data=samp).fit()
        boot_mat[i,:] = mod_b.params.values
    
    boot_df = pd.DataFrame(boot_mat, columns=params.index)
    
    # 3) Summary stats
    boot_se   = boot_df.std(ddof=1)
    ci_low    = boot_df.quantile(0.025)
    ci_high   = boot_df.quantile(0.975)
    
    # 4) Empirical p-values
    p_emp = {
        name: ((boot_df[name].apply(np.sign) != np.sign(params[name])).sum() + 1) / (B + 1)
        for name in params.index
    }
    
    # 5) Build summary table
    summary = pd.DataFrame({
        'estimate':    params,
        'boot_se':     boot_se,
        'ci_2.5%':     ci_low,
        'ci_97.5%':    ci_high,
        'p_empirical': pd.Series(p_emp)
    })
    
    # 6) exp(beta) and %change
    expb      = np.exp(params)
    pct       = (expb - 1)*100
    summary['exp_beta'] = expb
    summary['%change']  = pct
    
    # 7) clean intercept
    summary.loc['Intercept', ['exp_beta','%change']] = [np.nan, np.nan]
    
    # 8) reorder
    summary = summary[[
        'estimate','boot_se','ci_2.5%','ci_97.5%',
        'exp_beta','%change','p_empirical'
    ]]
    
    return summary, boot_df

def set_datatypes(dose_df, testing_column, eps=1e-6):
    """
    Cleans up the dose-table and response variables. Also does some needed
    processing, such as log-transforming. Integrate this part into the 
    dose-table one to make it cleaner if there is time.

    Parameters
    ----------
    dose_df : dataframe
        dataframe created in previous steps, with relevant experimental
        variables and doses.
    testing_column : string
        the response variable we want to test. Usually AUC, Mu, OD...
    eps : float, optional
        Just a static offset to make sure we dont
        log zero. The default is 1e-6.

    Returns
    -------
    dataframe
        A corrected dataframe that is directly useable with our GLM setup.

    """
    
    df = dose_df.copy()

    # Log‐transform response
    df['log_resp'] = np.log(df[testing_column].astype(float) + eps)

    # Treatment as ordered Categorical
    df['Treatment'] = df['TreatmentName'].fillna('None').astype(str)
    treats = ['None'] + [t for t in df['Treatment'].unique() if t != 'None']
    df['Treatment'] = pd.Categorical(df['Treatment'], categories=treats, ordered=True)

    # Numeric dose and flag
    df['Supplementation (per mM)']     = df['SupplementDose_mM'].astype(float)
    df['Negative control'] = df['NegControlFlag'].astype(int)

    # 6) Return exactly the columns for modeling + this new factor
    return df[['log_resp', 'Treatment', 'Supplementation (per mM)', 'Negative control']]


def plot_ready_df(summary_df, boot_df, dose, treat_name):
    """
    Codeblock to just reformat the bootstrap summaries to allow for equimolar
    comparison to the negative controls. Thought it safer to estimate directly
    from the model instead of multiplying by the dose.

    Parameters
    ----------
    summary_df : dataframe
        bootstrap summary.
    boot_df : dataframe
        estimated bootstrap metrics.
    dose : float
        the dose used to correct for equimolar comparison.
    treat_name : string

    Returns
    -------
    dataframe
        a corrected dataframe that we can directly plot as a forest plot.

    """
    
    df = summary_df.copy()
    rows = []

    def make_scaled(name, term):
        beta = df.loc[term, 'estimate']
        se   = df.loc[term, 'boot_se']
        lo   = df.loc[term, 'ci_2.5%']
        hi   = df.loc[term, 'ci_97.5%']
        p    = df.loc[term, 'p_empirical']

        est_s = beta * dose
        se_s  = se * dose
        lo_s  = lo * dose
        hi_s  = hi * dose

        expb  = np.exp(est_s)
        pct   = (expb - 1) * 100

        return pd.Series({
            'estimate':    est_s,
            'boot_se':     se_s,
            'ci_2.5%':     lo_s,
            'ci_97.5%':    hi_s,
            'exp_beta':    expb,
            '%change':     pct,
            'p_empirical': p
        }, name=name)

    slope_term   = 'Q("Supplementation (per mM)")'
    inter_term   = f'C(Treatment)[T.{treat_name}]:Q("Supplementation (per mM)")'

    for idx in df.index:
        # Always add the original row
        rows.append(df.loc[idx].copy().rename(idx))

        if idx == slope_term:
            rows.append(make_scaled(f"Supplement at {dose} mM", slope_term))

        if idx == inter_term:
            rows.append(make_scaled(
                f"{treat_name}×(Supplement at {dose} mM)",
                inter_term
            ))

    result = pd.DataFrame(rows)
    return result

def growth_testing(EXPERIMENT_DIR):
    
    with open(EXPERIMENT_DIR / 'protocol/protocol.json', 'r') as f:
        data = json.load(f)
    
    # Convert experiments list into a DataFrame
    experiments_df = pd.json_normalize(data['experiments'])

    # Load observables
    growth_data = pd.read_csv(EXPERIMENT_DIR / 'results/growth/processed/growth_parameters.tsv',
                              sep='\t', index_col=0)
    
    # Create design table from dispensing layout.
    layout = pd.read_excel(EXPERIMENT_DIR / 'protocol/hamilton/pipetting_layout.xlsx')
    layout.rename(columns={'well': 'Well'}, inplace=True)
    layout = layout.merge(experiments_df, left_on = 'Summary', right_on='summary', how='left')
    
    # Need this for later
    design_table = create_design_table(layout)
    design_table.to_csv(EXPERIMENT_DIR / 'protocol/plate_layout/design_table.tsv', sep = '\t')
    
    dose_table = create_dose_table(layout)
    
    data = dose_table.merge(growth_data.drop('Summary', axis=1), left_on='Well', right_index=True)
    data = data[data['Summary'] != 'Media control']
    
    # Define formulas
    f_lin = 'log_resp ~ C(Treatment) * (Q("Supplementation (per mM)") + Q("Negative control"))'

    for testing_column in ['AUC']:
        
        dat = set_datatypes(data, testing_column, eps=1e-2)

        negative_dose = data['Doses'][data['NegControlFlag'] == 1].max()
        treat_name = data['TreatmentName'][data['TreatmentName'] != 'None'].unique()[0] # janky
        
        # bootstrap
        summary, boot_df = bootstrap_summary(dat, f_lin, B=5000)
        
        summary.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}_summary.csv')
        boot_df.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}_boot.csv')
    
        #print(f"\n===== Results for {testing_column} =====")
        #print(summary.to_string())
        
        df_forest = plot_ready_df(summary, boot_df, negative_dose, treat_name)
        df_forest.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}_forest.csv')
        
        print(f"\n===== Dose corrected results for {testing_column} =====")
        print(df_forest.to_string())
        
       