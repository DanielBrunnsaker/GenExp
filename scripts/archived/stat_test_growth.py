#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 20 10:33:50 2024

@author: danbru
"""

import pandas as pd
import statsmodels.api as sm
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import seaborn as sns
from scipy import stats
import numpy as np

# Sorry for lazy functions
def assign_proline_level(summary):
    if 'Low proline' in summary:
        return 'Low'
    elif 'High proline' in summary:
        return 'High'
    elif 'no supplementation' in summary or 'Control' in summary:
        return 'None'
    elif 'without proline' in summary or 'Control' in summary:
        return 'None'
    return 'Unknown'

def assign_lactic_acid(summary):
    if 'with lactic acid' in summary:
        return 'Yes'
    #elif 'without lactic acid' in summary or 'Control' in summary:
    #    return 'No'
    elif 'Lactic acid treatment without proline supplementation.' in summary:
        return 'Yes'
    else:
        return 'No'
    
    
# Sorry for lazy functions
def assign_supplementation_level(summary):
    if 'low' in summary:
        return 'Low'
    elif 'high' in summary:
        return 'High'
    elif 'no supplementation' in summary or 'Control' in summary:
        return 'None'
    elif 'without proline' in summary or 'Control' in summary:
        return 'None'
    return 'Unknown'

def assign_treatment(summary):
    if 'with lactic acid' in summary:
        return 'Yes'
    #elif 'without lactic acid' in summary or 'Control' in summary:
    #    return 'No'
    elif 'Lactic acid treatment without proline supplementation.' in summary:
        return 'Yes'
    else:
        return 'No'

    
def lr_test(full_model, restricted_model):
    
    # Log-likelihoods
    llf_full = full_model.llf         
    llf_restricted = restricted_model.llf  

    lr_stat = -2.0 * (llf_restricted - llf_full)
    df_diff = (full_model.df_model) - (restricted_model.df_model)

    # p-value from chi-square distribution
    p_value = 1 - stats.chi2.cdf(lr_stat, df_diff)

    return lr_stat, df_diff, p_value
    

def create_design_table(df):
    design_df = df.copy()
    
    # Process Supplement Column
    unique_supplement_values = sorted(df['Supplement (uL)'].unique())
    supplement_mapping = {0: 'None', unique_supplement_values[1]: 'Low', unique_supplement_values[2]: 'High'}
    design_df['Supplement'] = df['Supplement (uL)'].map(supplement_mapping)
    
    # Process Treatment Column
    treatment_mapping = {0: 'No', df['Treatment (uL)'].max(): 'Yes'}
    design_df['Treatment'] = df['Treatment (uL)'].map(treatment_mapping)
    
    # Process Negative Control Column
    negative_mapping = {0: 'No', df['Negative Control (uL)'].max(): 'Yes'}
    design_df['Negative'] = df['Negative Control (uL)'].map(negative_mapping)
    
    # Keep only relevant columns
    design_df = design_df[['Well','Summary','Treatment', 'Supplement', 'Negative']]
    
    return design_df



import pandas as pd
import numpy as np
#from scipy.integrate import simps

# Function to compute AUC for each well
def compute_auc(df):
    # Exclude the last column (experiment descriptor)
    time_series_data = df
    time_points = np.array(time_series_data.columns, dtype=float)  # Convert column names to time points

    auc_results = []
    for well_id, row in time_series_data.iterrows():
        growth_values = row.astype(float).values  # Growth data
        
        # Compute AUC using Trapezoidal Rule
        auc_trapz = np.trapz(growth_values, time_points)
        auc_results.append([well_id, auc_trapz])

    # Create DataFrame with AUC results
    auc_df = pd.DataFrame(auc_results, columns=["Well", "AUC_Trapezoidal"])
    auc_df.set_index("Well", inplace=True)
    
    return auc_df


#/results/growth/processed

#exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/spermine3/glutamate_0.0_1_20250128_1618'
#exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/caffeine/selected/arginine_202502131014'
#amiga_outputs = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/lactic acid/run/proline_0.5_5_20241219_1102/results/summary/lactic_acid_exp_corrected_summary.txt', sep = '\t')
#amiga_outputs = pd.read_csv(f'{exp_path}/results/summary/growth_data_corrected_summary.txt', sep = '\t')
# Remove outliers
#amiga_outputs = amiga_outputs[amiga_outputs['K_Error > 20%'] == False]

growth_data = pd.read_csv(f'{exp_path}/results/growth/processed/growth_data_corrected.txt', sep = '\t', index_col = 0)
growth_auc = compute_auc(growth_data)
#growth_auc = np.log(growth_auc)

layout = pd.read_excel(f'{exp_path}/protocol/hamilton/pipetting_layout.xlsx')
layout.rename(columns={'well':'Well'}, inplace=True)

design_table = create_design_table(layout)

data = design_table.merge(growth_auc, left_on = 'Well', right_index = True)
averages = data.groupby('Summary').mean()
testing_column = 'AUC_Trapezoidal'
#averages['tGr(s)'] = averages['t_gr']*3600

data = data[data['Summary'] != 'Media control']

#rls_data = data[['OD_Emp_AUC', 'Treatment', 'Supplement','Negative']]
rls_data = data[[testing_column, 'Treatment', 'Supplement']]

# 1) Ensure variables are categorical
rls_data['Treatment'] = pd.Categorical(
    rls_data['Treatment'],
    categories=['None', 'Yes'],     # 'No' is the baseline
    ordered=True
)
rls_data['Supplement'] = pd.Categorical(
    rls_data['Supplement'],
    categories=['None', 'PosLow', 'PosHigh', 'NegHigh'],  # 'None' is baseline
    ordered=True
)
# Instead, lets try for a model based on the gamma dist instead
model_gamma_cat = smf.glm(
    formula=f'{testing_column} ~ C(Treatment) * C(Supplement)',
    data=rls_data,
    family=sm.families.Gamma(link=sm.families.links.log())
).fit()

print(model_gamma_cat.summary()) # THis output should be enough i think





data = design_table.merge(mu_df, left_on = 'Well', right_index = True)
averages = data.groupby('Summary').mean()
testing_column = 'mu'
#averages['tGr(s)'] = averages['t_gr']*3600

data = data[data['Summary'] != 'Media control']

#rls_data = data[['OD_Emp_AUC', 'Treatment', 'Supplement','Negative']]
rls_data = data[[testing_column, 'Treatment', 'Supplement']]

# 1) Ensure variables are categorical
rls_data['Treatment'] = pd.Categorical(
    rls_data['Treatment'],
    categories=['None', 'Yes'],     # 'No' is the baseline
    ordered=True
)
rls_data['Supplement'] = pd.Categorical(
    rls_data['Supplement'],
    categories=['None', 'PosLow', 'PosHigh', 'NegHigh'],  # 'None' is baseline
    ordered=True
)

# Instead, lets try for a model based on the gamma dist instead
model_gamma_cat = smf.glm(
    formula=f'{testing_column} ~ C(Treatment) * C(Supplement)',
    data=rls_data,
    family=sm.families.Gamma(link=sm.families.links.log())
).fit()

print(model_gamma_cat.summary()) # THis output should be enough i think
















predict_df = pd.DataFrame({
    'Treatment': ['None', 'None', 'None', 'Yes', 'Yes', 'Yes'],
    'Supplement': ['None', 'PosHigh', 'NegHigh', 'None', 'PosHigh', 'NegHigh']
})

predict_df['predicted'] = model_gamma_cat.predict(predict_df)
print(predict_df)


print(rls_data[testing_column].describe())
print(rls_data['Treatment'].value_counts())
print(rls_data['Supplement'].value_counts())

# THis part is just to test ratio of likelihoods using the chi-squared distribution
# That is, we compare models with the factors removed to see if they are meaningful. 
# Could be useful depending on the hypithesis ine wants to prove, 
# but i think all of the unfo we need is likely cntained in the original model
# Might be useful for questions like;
# I.e. Does proline help at all?
# Is the interaction meaningful?


# Model with only Lactic Acid (no Proline factor included in the model)
model_no_proline = smf.glm(
    f'{testing_column} ~ C(Treatment)',
    data=rls_data,
    family=sm.families.Gamma(link=sm.families.links.log())
).fit()

lr_stat, df_diff, p_value = lr_test(full_model=model_gamma_cat,
                                    restricted_model=model_no_proline)

print("Likelihood Ratio Test for adding supplement (and interaction):")
print(f"LR statistic = {lr_stat:.3f}, df = {df_diff}, p-value = {p_value:.3e}")


# Model with only Lactic Acid (no Proline factor included in the model)
model_no_treatment = smf.glm(
    f'{testing_column} ~ C(Supplement)',
    data=rls_data,
    family=sm.families.Gamma(link=sm.families.links.log())
).fit()

lr_stat, df_diff, p_value = lr_test(full_model=model_gamma_cat,
                                    restricted_model=model_no_treatment)

print("Likelihood Ratio Test for adding treatment (and interaction):")
print(f"LR statistic = {lr_stat:.3f}, df = {df_diff}, p-value = {p_value:.3e}")


# Additive model (no interaction included in the model):
model_no_interaction = smf.glm(
    f'{testing_column} ~ C(Treatment) + C(Supplement)',
    data=rls_data,
    family=sm.families.Gamma(sm.families.links.log())
).fit()

lr_stat_int, df_diff_int, p_value_int = lr_test(
    full_model=model_gamma_cat,
    restricted_model=model_no_interaction
)

print("Likeligood Ratio Test for the interaction term:")
print(f"LR statistic = {lr_stat_int:.3f}, df = {df_diff_int}, p-value = {p_value_int:.3e}")





#### Interaction plot, i promise chatgpt did not generate this code

# Example data aggregation (replace 'data' with your actual DataFrame)
agg_data = rls_data.groupby(['Proline_Level_Ordinal', 'Lactic_Acid']).agg(
    mean_AUC=('OD_Emp_AUC', 'mean'),
    std_AUC=('OD_Emp_AUC', 'std')
).reset_index()

# Set aesthetic style
sns.set_theme(style="whitegrid", context="talk")

# Create the plot
fig, ax = plt.subplots(figsize=(9, 5))

# Define colors for the groups
colors = sns.color_palette('deep', n_colors=2)
labels = ['No Treatment', 'Lactic Acid Treatment']

# Plot for each group
for i, group in enumerate(agg_data['Lactic_Acid'].unique()):
    subset = agg_data[agg_data['Lactic_Acid'] == group]
    ax.plot(
        subset['Proline_Level_Ordinal'],
        subset['mean_AUC'],
        label=labels[i],
        color=colors[i],
        marker='o',
        linestyle='-',
        linewidth=2
    )
    ax.fill_between(
        subset['Proline_Level_Ordinal'],
        subset['mean_AUC'] - subset['std_AUC'],  # Lower bound
        subset['mean_AUC'] + subset['std_AUC'],  # Upper bound
        color=colors[i],
        alpha=0.2
    )

# Customize the axes
ax.set_xticks([0, 1, 2])
ax.set_xticklabels(['None', 'Low', 'High'], fontsize=12)
ax.set_xlabel('Proline Level', fontsize=14)
ax.set_ylabel('Observed AUC', fontsize=14)

# Add a legend
ax.legend(
    title='Treatment Group',
    title_fontsize=13,
    fontsize=12,
    loc='upper left',
    frameon=True
)

# Add gridlines and title
ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
ax.set_title('Effect of Proline on Lactic Acid treatment (Growth AUC)', fontsize=16, pad=20)

# Show the plot
plt.tight_layout()
plt.show()



# Test the dist

# Fit the gamma distribution
from scipy.stats import gamma, norm, expon, lognorm, beta, weibull_min, kstest, loggamma
# List of distributions to compare
distributions = {
    "Gamma": gamma,
    "LogGamma": loggamma,
    "Normal": norm,
    "Exponential": expon,
    "Lognormal": lognorm,
    "Beta": beta,
    "Weibull": weibull_min
}
results = []

# Fit each distribution and calculate metrics
for name, dist in distributions.items():
    try:
        # Fit the distribution to data
        params = dist.fit(data['OD_Emp_AUC'])

        # K-S test using the cumulative distribution function (CDF)
        ks_stat, p_value = kstest(data['OD_Emp_AUC'], dist.cdf, args=params)

        # Log-likelihood
        log_likelihood = np.sum(dist.logpdf(data['OD_Emp_AUC'], *params))

        # Store results
        results.append((name, ks_stat, p_value, log_likelihood, params))
    except Exception as e:
        print(f"Could not fit {name}: {e}")

# Sort results by K-S statistic
results.sort(key=lambda x: x[1])

# Print results
print(f"{'Distribution':<15}{'K-S Stat':<10}{'P-Value':<10}{'Log-Likelihood':<15}")
for name, ks_stat, p_value, log_likelihood, params in results:
    print(f"{name:<15}{ks_stat:<10.4f}{p_value:<10.4f}{log_likelihood:<15.2f}")

# Plot histogram and fitted PDFs
plt.hist(data, bins=30, density=True, alpha=0.6, label="Data")

x = np.linspace(data['OD_Emp_AUC'].min(), data['OD_Emp_AUC'].max(), 1000)
for name, _, _, _, params in results:
    dist = distributions[name]
    pdf = dist.pdf(x, *params)
    plt.plot(x, pdf, label=name)

plt.legend()
plt.title("Comparison of Fitted Distributions")
plt.xlabel("Value")
plt.ylabel("Density")
plt.show()


