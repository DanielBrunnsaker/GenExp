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
from tqdm import tqdm
import json
from sklearn.utils import resample
from patsy import dmatrices
from pathlib import Path


def create_dose_table(df):
    """
    Create a dose table with clearer column names:
      - Well, Summary
      - TreatmentName (categorical)
      - TreatmentDose_uL (float)
      - SupplementName (categorical)
      - SupplementDose_mM (float)
      - NegativeControlVol_uL (float)
      - NegControlFlag (0/1)
    """
    d = df.copy()
    # 1) Rename raw columns to clear names
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
            # extract number
            s = str(x)
            m = pd.Series([s]).str.extract(r'(\d+(?:\.\d+)?)', expand=False)[0]
            return float(m) if pd.notna(m) else 0.0
        except Exception:
            return 0.0

    d['SupplementDose_mM']   = d['SupplementDose_mM'].apply(parse_numeric)
    d['TreatmentDose_uL']    = d['TreatmentDose_uL'].apply(parse_numeric)
    d['NegativeControlVol_uL']= d['NegativeControlVol_uL'].fillna(0).astype(float)

    # 3) Flag negative control by supplement name (e.g., 'NegHigh')
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
    Given df and a patsy formula, do:
      1. Fit OLS to extract params
      2. Bootstrap B resamples of that fit
      3. Compute bootstrap SE, 95% CI, empirical p-values
      4. Attach exp(beta) & % change, with NaN for the intercept
    """
    
    np.random.seed(0)
    
    # 1) Fit once to get param names
    orig_mod    = smf.ols(formula, data=df).fit()
    params      = orig_mod.params
    
    # 2) Bootstrap
    boot_mat = np.zeros((B, len(params)))
    for i in tqdm(range(B), desc="Bootstrapping"):
        samp    = resample(df)
        mod_b   = smf.ols(formula, data=samp).fit()
        boot_mat[i,:] = mod_b.params.values
    
    boot_df = pd.DataFrame(boot_mat, columns=params.index)
    
    # 3) Summaries
    boot_se   = boot_df.std(ddof=1)
    ci_low    = boot_df.quantile(0.025)
    ci_high   = boot_df.quantile(0.975)
    
    # 4) Empirical p-values
    p_emp = {
        name: ((boot_df[name].apply(np.sign) != np.sign(params[name])).sum() + 1) / (B + 1)
        for name in params.index
    }
    
    # 5) Build table
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
    Prepare dose_df for modeling, producing:
      - log_resp
      - Treatment (categorical None/Yes)
      - dose_mM (float)
      - neg_control (int flag)
      - SupplementLevel (categorical: None, Low, High, Negative)
    
    Automatically assigns 'Low' and 'High' based on the sorted unique positive doses.
    """
    df = dose_df.copy()

    # 1) Log‐transform response
    df['log_resp'] = np.log(df[testing_column].astype(float) + eps)

    # 2) Treatment as ordered Categorical
    df['Treatment'] = df['TreatmentName'].fillna('None').astype(str)
    treats = ['None'] + [t for t in df['Treatment'].unique() if t != 'None']
    df['Treatment'] = pd.Categorical(df['Treatment'], categories=treats, ordered=True)

    # 3) Numeric dose and flag
    df['Supplementation (per mM)']     = df['SupplementDose_mM'].astype(float)
    df['Negative control'] = df['NegControlFlag'].astype(int)

    # 4) Determine Low vs High thresholds
    # get unique positive doses (exclude zeros and negatives)
    pos_doses = sorted(df.loc[df['Supplementation (per mM)'] > 0, 'Supplementation (per mM)'].unique())
    if len(pos_doses) >= 2:
        low_val, high_val = pos_doses[0], pos_doses[-1]
    else:
        low_val, high_val = None, None

    # 5) Build 4‐level SupplementLevel
    def label_sup(row):
        if row['Negative control'] == 1:
            return 'Negative'
        dm = row['Supplementation (per mM)']
        if dm == 0:
            return 'None'
        if low_val is not None and np.isclose(dm, low_val):
            return 'Low'
        if high_val is not None and np.isclose(dm, high_val):
            return 'High'
        # fallback for unexpected
        return 'High' if dm > low_val else 'Low'

    df['SupplementLevel'] = df.apply(label_sup, axis=1)
    levels = ['None', 'Low', 'High', 'Negative']
    df['SupplementLevel'] = pd.Categorical(df['SupplementLevel'],
                                           categories=levels,
                                           ordered=False)

    # 6) Return exactly the columns for modeling + this new factor
    return df[['log_resp', 'Treatment', 'Supplementation (per mM)', 'Negative control', 'SupplementLevel']]

def plot_ready_df(summary_df, boot_df, dose, treat_name):
    """
    Same as before, but injects the two scaled rows immediately after
    the per-mM slope row and its interaction parent, preserving order.
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

    # Define the two parent terms
    slope_term   = 'Q("Supplementation (per mM)")'
    inter_term   = f'C(Treatment)[T.{treat_name}]:Q("Supplementation (per mM)")'

    for idx in df.index:
        # 1) Always add the original row
        rows.append(df.loc[idx].copy().rename(idx))

        # 2) Right after the slope term, inject "Supplement at dose"
        if idx == slope_term:
            rows.append(make_scaled(f"Supplement at {dose} mM", slope_term))

        # 3) Right after the interaction term, inject "Treatment×(Supplement at dose)"
        if idx == inter_term:
            rows.append(make_scaled(
                f"{treat_name}×(Supplement at {dose} mM)",
                inter_term
            ))

    # Reassemble into a DataFrame
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
    
    dose_table = create_dose_table(layout)
    
    data = dose_table.merge(growth_data.drop('Summary', axis=1), left_on='Well', right_index=True)
    data = data[data['Summary'] != 'Media control']
    
    # 1) Define formulas
    f_lin = 'log_resp ~ C(Treatment) * (Q("Supplementation (per mM)") + Q("Negative control"))'
    #f_cat = 'log_resp ~ C(Treatment)*C(SupplementLevel)'

    for testing_column in ['AUC']:
        
        dat = set_datatypes(data, testing_column, eps=1e-2)
        
        #mod_lin = smf.ols(f_lin, data=dat).fit()
        #mod_cat = smf.ols(f_cat, data=dat).fit()

        # select which formula to use based on the fit
        #chosen_formula = f_lin if 2*mod_lin.aic <= mod_cat.aic else f_cat
        
        #print("Selected:", chosen_formula)
        
        negative_dose = data['Doses'][data['NegControlFlag'] == 1].max()
        treat_name = data['TreatmentName'][data['TreatmentName'] != 'None'].unique()[0] # janky
        
        # bootstrap
        summary, boot_df = bootstrap_summary(dat, f_lin, B=5000)
        
        summary.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}_summary.csv')
        boot_df.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}_boot.csv')
    
        print(f"\n===== Results for {testing_column} =====")
        print(summary.to_string())
        
        df_forest = plot_ready_df(summary, boot_df, negative_dose, treat_name)
        df_forest.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}_forest.csv')
        
        print(f"\n===== Dose corrected results for {testing_column} =====")
        print(df_forest.to_string())
        
        #if chosen_formula != 'log_resp ~ C(Treatment)*C(SupplementLevel)':
        #dose_corrected_summary = corrected_summary(summary, dat)
        #else:
        #    dose_corrected_summary = unify_index_names(summary, dat)

        #print(f"\n===== Corrected for dose =====")
        #print(dose_corrected_summary.to_string())

        #dose_corrected_summary.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}_summary_namingcorr.csv')
