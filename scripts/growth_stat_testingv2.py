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
from tqdm import tqdm

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

def create_dose_table(df):
    """
    Create a dose table using standardized column names:
      - media_supplementation       : supplement name (e.g. None, 'X')
      - media_supplementation_doses : numeric dose for supplement (mM)
      - treatment                   : treatment name (e.g. None, 'T')
      - treatment_parameters        : numeric treatment dose
      - Well, Summary               : identifiers
    
    Returns a DataFrame with:
      - Well
      - Summary
      - Treatment (categorical)
      - treat_dose (float or None)
      - Supplement (categorical)
      - dose_mM (float or None)
    with NaNs replaced by None.
    """
    dose_df = df.copy()
    
    # Verify required columns
    required = ['Well', 'Summary', 'media_supplementation', 
                'media_supplementation_doses', 'treatment', 'treatment_parameters']
    missing = set(required) - set(dose_df.columns)
    if missing:
        raise KeyError(f"Missing required columns: {missing}")
    
    # Rename & select
    dose_df = dose_df[required].rename(columns={
        'media_supplementation': 'Supplement',
        'media_supplementation_doses': 'dose_mM',
        'treatment': 'Treatment',
        'treatment_parameters': 'treat_dose'
    })
    
    # Ensure categorical types
    dose_df['Supplement'] = pd.Categorical(dose_df['Supplement'], ordered=False)
    dose_df['Treatment'] = pd.Categorical(dose_df['Treatment'], ordered=False)
    
    # Replace NaN with None for better downstream compatibility
    dose_df = dose_df.where(pd.notnull(dose_df), None)
    
    return dose_df[['Well', 'Summary', 'Treatment', 'treat_dose', 'Supplement', 'dose_mM']]


def create_dose_table(df):
    """
    Build your dose table, using the explicit Negative Control column.
    Expects df (the merged design + growth_data) to contain:
      - 'Well', 'Summary'
      - 'media_supplementation', 'media_supplementation_doses'
      - 'treatment', 'treatment_parameters'
      - 'Negative Control (uL)'      <-- use this directly
    Returns columns:
      Well, Summary,
      Treatment (categorical),
      treat_dose (float),
      Supplement (categorical),
      dose_mM (float),
      neg_control (0/1)
    """
    d = df.copy()
    # rename
    d = d.rename(columns={
        'media_supplementation':       'Supplement',
        'media_supplementation_doses': 'dose_mM',
        'treatment':                   'Treatment',
        'treatment_parameters':        'treat_dose',
        'Negative Control (uL)':       'neg_vol'
    })
    # parse dose_mM and treat_dose as before
    d['dose_mM'] = (
        d['dose_mM'].fillna('0').astype(str)
         .str.extract(r'(\d+(?:\.\d+)?)', expand=False)
         .astype(float).fillna(0.0)
    )
    d['treat_dose'] = (
        d['treat_dose'].fillna('0').astype(str)
         .str.extract(r'(\d+(?:\.\d+)?)', expand=False)
         .astype(float).fillna(0.0)
    )
    # explicitly flag negative control wells
    d['neg_control'] = (d['neg_vol'].fillna(0).astype(float) > 0).astype(int)
    
    # null out Treatment when treat_dose==0
    d['Treatment'] = d['Treatment'].where(d['treat_dose'] > 0, 'None').fillna('None')
    # category
    treats = ['None'] + [t for t in d['Treatment'].unique() if t!='None']
    d['Treatment'] = pd.Categorical(d['Treatment'], categories=treats, ordered=True)
    # Supplement as category
    d['Supplement'] = pd.Categorical(d['Supplement'])
    
    # drop the floating `neg_vol`
    return d[['Well','Summary','Treatment','treat_dose',
              'Supplement','dose_mM','neg_control']]



def set_datatypes(data, testing_column):
    
    rls_data = data[[testing_column, 'Treatment', 'Supplement']]
    rls_data[testing_column] = np.log(rls_data[testing_column] + 0.01)
    
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

def bootstrap_summary(df, formula, B=5000, eps=1e-6):
    """
    Given df and a patsy formula, do:
      1. Fit OLS to extract params
      2. Bootstrap B resamples of that fit
      3. Compute bootstrap SE, 95% CI, empirical p-values
      4. Attach exp(beta) & % change, with NaN for the intercept
    """
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
    return summary


def set_datatypes(dose_df, testing_column, eps=1e-6):
    """
    Prepare the merged dose_table for modeling:
      - Log-transform response column testing_column
      - Ensure Treatment is categorical (None first)
      - Pass through dose_mM (float) and neg_control (int)
    
    Expects `dose_df` (from create_dose_table) to already contain:
      - testing_column (numeric raw response)
      - 'Treatment' (categorical with 'None' and real names)
      - 'dose_mM' (float)
      - 'neg_control' (0/1)
    """
    df = dose_df.copy()
    
    # 1) Log-transform the response
    df['log_resp'] = np.log(df[testing_column].astype(float) + eps)
    
    # 2) Ensure Treatment is categorical with 'None' first
    #    (If it's already categorical, we'll re-establish the correct order)
    unique_treats = ['None'] + [t for t in df['Treatment'].cat.categories if t != 'None']
    df['Treatment'] = pd.Categorical(df['Treatment'], 
                                     categories=unique_treats, 
                                     ordered=True)
    
    # 3) dose_mM and neg_control are already numeric—just pass them through
    df['dose_mM']     = df['dose_mM'].astype(float)
    df['neg_control'] = df['neg_control'].astype(int)
    
    # 4) Return only the four columns your model needs
    return df[['log_resp', 'Treatment', 'dose_mM', 'neg_control']]





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
    
    #design_table = create_design_table(layout)
    dose_table = create_dose_table(layout)
    
    # Save design table 
    #design_table.to_csv(EXPERIMENT_DIR / 'protocol/plate_layout/design_table.tsv', sep = '\t')
    
    #data = design_table.merge(growth_data.drop('Summary', axis=1), left_on='Well', right_index=True)
    #data = data[data['Summary'] != 'Media control']
    
    data = dose_table.merge(growth_data.drop('Summary', axis=1), left_on='Well', right_index=True)
    data = data[data['Summary'] != 'Media control']
    
    # Loop over each metric
    #for testing_column in ['AUC', 'mu', 'MaxOD']:
        
        
    # 0) configure
    B = 50
    
    
    for testing_column in ['AUC','MaxOD']:
        
        dat = set_datatypes(data, testing_column, eps=1e-2)
        
        #dat = set_datatypes(data, testing_column)
        # now dat has columns: log_resp, Treatment, dose_mM, neg_control
        
        # 1) Fit & bootstrap as before, but with this formula:
        formula = 'log_resp ~ C(Treatment) * (dose_mM + neg_control)'
        summary = bootstrap_summary(dat, formula, B=1000)
        print(summary.to_string())


        mod = smf.ols(formula, data=dat).fit()
        print(mod.params)    
    
    
    responses = ['AUC','MaxOD']
    
    # 1) build dose table
    #dose_table = create_dose_table(raw_df)
    
    for resp in responses:
        # 2) log-transform your response
        df = dose_table.copy()
        df['log_resp'] = np.log(df[resp] + 0.01)
        
        # 3) specify the model: log_resp ~ Treatment * dose_mM  (+neg_control if you have it)
        formula = 'log_resp ~ C(Treatment)*dose_mM'
        
        # 4) run bootstrap summary
        summary = bootstrap_summary(df, formula, B=B)
        
        print(f"\n===== Results for {resp} =====")
        print(summary)    
        
        
        
        
    
        
    from sklearn.utils import resample
    B = 1000    
    
    for testing_column in ['AUC', 'MaxOD', 'mu']:

        categorized_data = set_datatypes(data, testing_column)
        
        # 1) Dummy fit to extract length of parameters
        formula = f'{testing_column} ~ C(Treatment) * C(Supplement)'
        orig_mod = smf.ols(formula, data=categorized_data).fit()
        orig_params = orig_mod.params
        
        # 2) Bootstrap
        boot_params = np.zeros((B, len(orig_params)))
        for i in tqdm(range(B)):
            boot = resample(categorized_data)                      # sample rows with replacement
            mod_b = smf.ols(formula, data=boot).fit() # OLS fit on bootstrap sample
            boot_params[i, :] = mod_b.params.values
        
        boot_df = pd.DataFrame(boot_params, columns=orig_params.index)
        
        # 3) Summarize bootstrap distribution
        boot_se   = boot_df.std(ddof=1)
        ci_lower  = boot_df.quantile(0.025)
        ci_upper  = boot_df.quantile(0.975)
        
        # 4) Empirical two-sided p-values
        p_emp = {}
        for name in orig_params.index:
            vals = boot_df[name]
            p_emp[name] = ((np.sum(np.sign(vals) != np.sign(orig_params[name]))) + 1) / (B + 1)
        
        # 5) Build summary table
        summary = pd.DataFrame({
            'estimate':     orig_params,
            'exp_beta':     np.exp(orig_params),
            '%_change':     (np.exp(orig_params) - 1) * 100,
            'boot_se':      boot_se,
            'ci_2.5%':      ci_lower,
            'ci_97.5%':     ci_upper,
            'p_empirical':  pd.Series(p_emp)
        })
        
        print(f'----------------- {testing_column} -------------------')
        print(summary)
        
        
        # Save model coefficients and the fitted model.
        #results_df = pd.DataFrame({
        #    "Coefficient": model.params,
        #    "Std Error": model.bse,
        #    "t-value": model.tvalues,
        #    "p-value": model.pvalues
        #})
        summary.to_csv(EXPERIMENT_DIR / f'results/growth/tests/{testing_column}.tsv', sep='\t')
        #model.save(EXPERIMENT_DIR / f'results/growth/tests/models/{testing_column}.pickle')

    #print(model.summary())
        
        
        
        #if testing_column == 'AUC':
            # For AUC, use a Gamma GLM with log link. Did some empirical testing here, and this one seemed the best.
            # model = smf.glm(
             #    formula=f'{testing_column} ~ C(Treatment) * C(Supplement)',
              #  data=categorized_data,
                #family=sm.families.Gamma(link=sm.families.links.Log())
              #  family=sm.families.Tweedie(link=sm.families.links.Log(), var_power=1.5)
                #family=sm.families.Gaussian(link=sm.families.links.Identity())
            #).fit()
            
            #model = smf.glm(
            #    formula=f'{testing_column} ~ C(Treatment) * C(Supplement)',
            #    data=categorized_data,
            #    family=sm.families.Tweedie(
            #        link=sm.families.links.Log(),
            #        var_power=2    # let it estimate the mean–variance power p
            #    )
            #).fit(
            #    cov_type="HC3",  # robust SEs for any leftover heteroscedasticity
            #    maxiter=100
            #)
            
            
        
            
        '''
        def fit_best_tweedie(df, formula, p_grid=None):
            if p_grid is None:
                # search between Poisson-ish and IG-ish
                p_grid = np.linspace(1.0, 2.5, 31)
            
            best = {'p': None, 'disp': np.inf, 'model': None}
            for p in p_grid:
                fam = sm.families.Tweedie(link=sm.families.links.Log(), var_power=p)
                res = smf.glm(formula, df, family=fam).fit(disp=False)
                disp = res.pearson_chi2 / res.df_resid
                if abs(disp - 1) < abs(best['disp'] - 1):
                    best.update({'p': p, 'disp': disp, 'model': res})
            
            # finally refit that best model with HC3 standard errors
            fam = sm.families.Tweedie(link=sm.families.links.Log(), var_power=best['p'])
            final = smf.glm(formula, df, family=fam).fit(cov_type="HC3", maxiter=100)
            print(f"Selected Tweedie p ≈ {best['p']:.2f}, dispersion ≈ {best['disp']:.2f}")
            return final
        
        #model = fit_best_tweedie(categorized_data, f'{testing_column} ~ C(Treatment) * C(Supplement)')
        #model_type = "Tweedie GLM (log link)"
        
        model = smf.glm(
          f'{testing_column} ~ C(Treatment) * C(Supplement)',
          data=categorized_data,
          family=sm.families.Gaussian(link=sm.families.links.Identity())
        ).fit(cov_type="HC3")
        model_type = "Gaussian GLM (log transformed data)"
            # compute HC3‐robust covariances (you can also pick HC0, HC1, HC2…)
            #model = model._get_robustcov_results(cov_type="HC3")
            
            #print(model._get_robustcov_results(cov_type="HC3"))
            
            
        from statsmodels.regression.quantile_regression import QuantReg
        #import pandas as pd
        #import numpy as np
        
        # prepare design
        q_model = QuantReg.from_formula(f'{testing_column} ~ C(Treatment) * C(Supplement)', categorized_data)
        
        # fit median (q=0.5) – you can also try q=0.25, 0.75
        model = q_model.fit(q=0.5)
        #print(res_q50.summary())
            
        
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
        '''
        # Print final model summary for the current metric.
        #print(f'\n########### Final Testing for {testing_column} using {model_type} ###########\n')
        #print(model.summary())
        '''
        # Choose residuals for diagnostic tests.
        #if testing_column in ['AUC', 'mu']:
        resid = model.resid_response  # For GLMs, using the response residuals.
        bp_resid = model.resid_pearson
        #else:
        #    resid = model.resid  # For OLS.
        #    bp_resid = resid
        
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

