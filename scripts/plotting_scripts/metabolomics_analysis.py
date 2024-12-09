#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 10:50:38 2024

@author: danbru
"""
import pandas as pd
import matplotlib.pyplot as plt
import math
import numpy as np
import seaborn as sns
from itertools import combinations
from scipy.stats import ttest_ind, mannwhitneyu
from statannotations.Annotator import Annotator

# Load in metabolomics data for FA first?
### Get OD

#/Users/danbru/Downloads/AUTOMATED___no barcode detected_.DAT
od_data = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_0.25_5_20241121_1256 copy/results/growth/FA_corrected.txt', sep = '\t', index_col = 0)
od_data = pd.DataFrame(od_data.iloc[:,-1])

#od_data.columns = [['Well','OD']]
#od_data['Well'] = od_data['Well'].replace(r'(\D)0+(\d+)', r'\1\2', regex=True)
od_data.index = od_data.index.str.replace(r'(\D)0+(\d+)', r'\1\2', regex=True)
od_data.columns = ['OD']
od_data[od_data['OD'] < 0] = 0


# Compare the active wells?
layout = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_0.25_5_20241121_1256 copy/protocol/plate_layout/layout.tsv', sep = '\t')
layout['well'] = layout['well'].str.replace(r'(\D)0+(\d+)', r'\1\2', regex=True)

# Only keep the relevant wells
layout = layout[layout['CONCuM'].isin(['L-glutamate: 0 mM & Formic acid: 0 mM', 'L-glutamate: 5 mM & Formic acid: 0 mM', 'L-glutamate: 15 mM & Formic acid: 0 mM'])]
group_order = ['L-glutamate: 0 mM & Formic acid: 0 mM', 'L-glutamate: 5 mM & Formic acid: 0 mM', 'L-glutamate: 15 mM & Formic acid: 0 mM']  # Replace with your desired order
custom_group_names = {
    "L-glutamate: 0 mM & Formic acid: 0 mM": "No Glutamate",
    "L-glutamate: 5 mM & Formic acid: 0 mM": "Low Glutamate",
    "L-glutamate: 15 mM & Formic acid: 0 mM": "High Glutamate"
}

TIC = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/Results/faPOSCHROM.tsv', sep = '\t')
TIC = TIC[TIC['FragmentIon'] == 'Summed'][['FileName','TotalArea']]


data = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/Results/faPOS.csv')
data['Well'] = data['Replicate'].str.extract(r'Expt1-([A-H]\d{1,2})')
data['Peptide'] = data['Peptide'].str.replace('L-Glutamic acid','L-Glutamic Acid')

# Map based on matching "Replicate" substring in "FileName"
data['TotalArea'] = data['Replicate'].apply(
    lambda replicate: TIC.loc[TIC['FileName'].str.contains(replicate), 'TotalArea'].values[0]
)

data['Normalized_Area'] = data.apply(
    lambda row: row['Area'] / row['TotalArea'] if pd.notnull(row['Area']) and row['TotalArea'] != 0 else np.nan,
    axis=1
)

# Map Biomass (OD) from od_data based on matching "Well"
data['Biomass'] = data['Well'].map(
    od_data['OD']
)

# Apply the combined normalization formula:
# Final Normalized Area = (Area - Blank) / (TIC * Biomass)
data['Normalized_Area'] = data.apply(
    lambda row: (row['Area'] / (row['TotalArea'] * row['Biomass']))
    if pd.notnull(row['Area']) and row['TotalArea'] > 0 and row['Biomass'] > 0 else np.nan,
    axis=1
)

# First, select an adduct for every amino acid
peptides_to_remove = data.groupby('Peptide')['Normalized_Area'].apply(
    lambda x: (x == 0).all() or x.isna().all()
)
data = data[~data['Peptide'].isin(peptides_to_remove[peptides_to_remove].index)]

# Filter out rows with zero or NA in the 'Area' column
df_filtered = data[data['Normalized_Area'].notna() & (data['Normalized_Area'] != 0)]

# Group by 'Peptide', 'Protein', 'Precursor Mz' and count non-zero, non-NA values in the 'Area' column
grouped = df_filtered.groupby(['Peptide', 'Protein', 'Precursor Mz'])['Normalized_Area'].count().reset_index()

# Find the combination with the maximum count for each unique 'Peptide'
#max_combinations = grouped.loc[grouped.groupby(['Peptide','Protein','Precursor Mz'])['Normalized_Area'].idxmax()]

# Find the combination with the maximum count for each unique 'Peptide'
max_combinations = grouped.loc[grouped.groupby('Peptide')['Normalized_Area'].idxmax()]

# Filter the original DataFrame to only include these combinations
reduced_data = data.merge(
    max_combinations[['Peptide', 'Protein', 'Precursor Mz']],
    on=['Peptide', 'Protein', 'Precursor Mz'],
    how='inner'
)


subset = reduced_data[reduced_data['Peptide'] == 'L-Glutamic Acid']

# Merge the filtered DataFrame with the layout DataFrame based on the 'Well' column
merged_df = reduced_data.merge(layout, left_on='Well', right_on='well')
# Remove all zero and nan rows?
merged_df = merged_df[merged_df['Normalized_Area'] > 0]



# Map the custom names to the 'CONCuM' column
merged_df['Experimental group'] = merged_df['CONCuM'].map(custom_group_names)


# Create the grouped barplot
plt.figure(figsize=(12, 6))
barplot = sns.barplot(
    x='Peptide',
    y='Normalized_Area',
    hue='Experimental group',
    data=merged_df,
    ci='sd',
    hue_order=[custom_group_names[g] for g in group_order]
)

# Customize x-axis and y-axis
plt.title('TIC-Normalized Area by Experimental Group for Each Amino Acid (Log Scale)')
plt.xlabel('Peptide')
plt.ylabel('Mean Area (Log Scale)')
plt.xticks(rotation=45, ha='right')  # Rotate and right-align x-ticks
plt.yscale('log')  # Apply log scale to y-axis

# Prepare pairs for statistical comparison
peptides = merged_df['Peptide'].unique()
pairs = []
for peptide in peptides:
    groups_in_peptide = merged_df[merged_df['Peptide'] == peptide]['Experimental group'].unique()
    group_pairs = list(combinations(groups_in_peptide, 2))
    pairs.extend([((peptide, pair[0]), (peptide, pair[1])) for pair in group_pairs])

# Define p-value thresholds and corresponding symbols
pvalue_thresholds = [
    (0.001, "****"),
    (0.01, "***"),
    (0.05, "**"),
    (0.1, "*")
]

# Initialize the Annotator
annotator = Annotator(
    barplot,
    pairs,
    data=merged_df,
    x='Peptide',
    y='Area',
    hue='Experimental group',
    hue_order=[custom_group_names[g] for g in group_order]
)

# Configure and Apply Annotations
annotator.configure(
    #test='t-test_ind',
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    alpha = 0.1,
    comparisons_correction=None,
    pvalue_thresholds=pvalue_thresholds,
    hide_non_significant=True,  # Hide non-significant comparisons
    line_offset=0.01,  # Adjust line offset to move lines closer to bars
    text_offset=0.01,  # Adjust text offset to move stars closer to lines
    line_height=0.01,  # Adjust line height
    #line_height_units='axes_fraction'  # Units for offsets ('axes_fraction', 'points', 'data')
)
annotator.apply_and_annotate()

plt.tight_layout()
plt.show()

# Create a DataFrame with the comparison results
from statsmodels.stats.multitest import multipletests

comparison_results = []

for peptide in peptides:
    
    peptide = 'L-Arginine'
    
    peptide_data = merged_df[merged_df['Peptide'] == peptide]
    groups = peptide_data['CONCuM'].unique()
    group_pairs = list(combinations(groups, 2))
    for pair in group_pairs:
        group1_data = peptide_data[peptide_data['CONCuM'] == pair[0]]['Normalized_Area']
        group2_data = peptide_data[peptide_data['CONCuM'] == pair[1]]['Normalized_Area']
        
        if group1_data.dropna().shape[0] == 0:
            continue
        
        #t_stat, p_value = ttest_ind(group1_data, group2_data, nan_policy='omit')
        t_stat, p_value = mannwhitneyu(group1_data, group2_data, nan_policy='omit')
        
        comparison_results.append({
            'Peptide': peptide,
            'Group1': custom_group_names.get(pair[0], pair[0]),
            'Group2': custom_group_names.get(pair[1], pair[1]),
            'P-Value': p_value
        })

comparison_df = pd.DataFrame(comparison_results)
# Apply BH correction within each peptide
comparison_df['Adjusted P-Value'] = comparison_df.groupby('Peptide')['P-Value'].transform(
    lambda x: multipletests(x, method='fdr_bh')[1]
)








# Biomass gradient!
TIC = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/Results/gradient_posCHROM.tsv', sep = '\t')
TIC = TIC[TIC['FragmentIon'] == 'Summed'][['FileName','TotalArea']]

data = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/Results/GradientPOS.csv')
data['Well'] = data['Replicate'].str.extract(r'test-([A-H]\d{1,2})')
data['Peptide'] = data['Peptide'].str.replace('L-Glutamic acid','L-Glutamic Acid')

# Map based on matching "Replicate" substring in "FileName"
data['TotalArea'] = data['Replicate'].apply(
    lambda replicate: TIC.loc[TIC['FileName'].str.contains(replicate), 'TotalArea'].values[0]
)

data['Normalized_Area'] = data.apply(
    lambda row: row['Area'] / row['TotalArea'] if pd.notnull(row['Area']) and row['TotalArea'] != 0 else np.nan,
    axis=1
)

# First, select an adduct for every amino acid
peptides_to_remove = data.groupby('Peptide')['Normalized_Area'].apply(
    lambda x: (x == 0).all() or x.isna().all()
)
data = data[~data['Peptide'].isin(peptides_to_remove[peptides_to_remove].index)]

# Filter out rows with zero or NA in the 'Area' column
df_filtered = data[data['Normalized_Area'].notna() & (data['Normalized_Area'] != 0)]

# Group by 'Peptide', 'Protein', 'Precursor Mz' and count non-zero, non-NA values in the 'Area' column
grouped = df_filtered.groupby(['Peptide', 'Protein', 'Precursor Mz'])['Normalized_Area'].count().reset_index()

# Find the combination with the maximum count for each unique 'Peptide'
#max_combinations = grouped.loc[grouped.groupby(['Peptide','Protein','Precursor Mz'])['Normalized_Area'].idxmax()]

# Find the combination with the maximum count for each unique 'Peptide'
max_combinations = grouped.loc[grouped.groupby('Peptide')['Normalized_Area'].idxmax()]

# Filter the original DataFrame to only include these combinations
reduced_data = data.merge(
    max_combinations[['Peptide', 'Protein', 'Precursor Mz']],
    on=['Peptide', 'Protein', 'Precursor Mz'],
    how='inner'
).dropna(subset = ['Normalized_Area','Well'])[['Peptide','Well', 'Area','Normalized_Area']]



df = reduced_data
# Extract the numeric column from the Well and sort the x-axis numerically
df['Column'] = df['Well'].str[1:].astype(int)

# Group by Peptide and Column, then calculate the average Normalized_Area
grouped = df.groupby(['Peptide', 'Column'])['Normalized_Area'].mean().reset_index()

# Pivot the data to structure it for plotting
pivot_df = grouped.pivot(index='Column', columns='Peptide', values='Normalized_Area')

# Plotting
plt.figure(figsize=(10, 6))
for peptide in pivot_df.columns:
    plt.plot(pivot_df.index, pivot_df[peptide], marker='o', label=peptide)

plt.title('Average Normalized Area Across Wells')
plt.xlabel('Well Column')
plt.ylabel('Average Normalized Area')
plt.xticks(sorted(pivot_df.index))  # Ensure the x-axis ticks are sorted numerically
plt.legend(title='Peptide')
plt.grid(True)
plt.show()






# Changing it to a featuretable?


# Compare the active wells?
layout = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_0.25_5_20241121_1256 copy/protocol/plate_layout/layout.tsv', sep = '\t')
layout['well'] = layout['well'].str.replace(r'(\D)0+(\d+)', r'\1\2', regex=True)

# Only keep the relevant wells
layout = layout[layout['CONCuM'].isin(['L-glutamate: 0 mM & Formic acid: 0 mM', 'L-glutamate: 5 mM & Formic acid: 0 mM', 'L-glutamate: 15 mM & Formic acid: 0 mM'])]
group_order = ['L-glutamate: 0 mM & Formic acid: 0 mM', 'L-glutamate: 5 mM & Formic acid: 0 mM', 'L-glutamate: 15 mM & Formic acid: 0 mM']  # Replace with your desired order
custom_group_names = {
    "L-glutamate: 0 mM & Formic acid: 0 mM": "No Glutamate",
    "L-glutamate: 5 mM & Formic acid: 0 mM": "Low Glutamate",
    "L-glutamate: 15 mM & Formic acid: 0 mM": "High Glutamate"
}

TIC = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/Results/faPOSCHROM.tsv', sep = '\t')
TIC = TIC[TIC['FragmentIon'] == 'Summed'][['FileName','TotalArea']]

data = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/Results/faPOS.csv')
data['Well'] = data['Replicate'].str.extract(r'Expt1-([A-H]\d{1,2})')
data['Peptide'] = data['Peptide'].str.replace('L-Glutamic acid','L-Glutamic Acid')
data = data[~data['Replicate'].str.startswith(('009', '010', '011', '012'))]


# First, select an adduct for every amino acid (pick our specific comparison)
peptides_to_remove = data[data['Well'].isin(list(layout['well']))].groupby('Peptide')['Area'].apply(
    lambda x: (x == 0).all() or x.isna().all()
)
data = data[~data['Peptide'].isin(peptides_to_remove[peptides_to_remove].index)]

# Filter out rows with zero or NA in the 'Area' column
df_filtered = data[data['Area'].notna() & (data['Area'] != 0)]

# Group by 'Peptide', 'Protein', 'Precursor Mz' and count non-zero, non-NA values in the 'Area' column
grouped = df_filtered.groupby(['Peptide', 'Protein', 'Precursor Mz'])['Area'].count().reset_index()

# Find the combination with the maximum count for each unique 'Peptide'
max_combinations = grouped.loc[grouped.groupby('Peptide')['Area'].idxmax()]

# Filter the original DataFrame to only include these combinations
reduced_data = data.merge(
    max_combinations[['Peptide', 'Protein', 'Precursor Mz']],
    on=['Peptide', 'Protein', 'Precursor Mz'],
    how='inner'
)[['Peptide','Replicate','Area']]

reduced_data['Area'] = reduced_data['Area'].fillna(0)

# Step 1: Define a function to extract the well identifier
def extract_well(replicate):
    if replicate.endswith('MAT1'):  # For blanks (MAT1 entries), use the first number
        #return 'blank-'+replicate.split('-')[0]
        return 'blank'
    else:  # For other entries, use the last part (e.g., A1, B7, H12)
        return replicate.split('-')[-1]

# Apply the function to create a 'Well' column
reduced_data['Well'] = reduced_data['Replicate'].apply(extract_well)

# Step 2: Handle duplicates by aggregating (e.g., summing the Area)
aggregated_data = reduced_data.groupby(['Well', 'Peptide'], as_index=False).agg({'Area': 'mean'})

# Step 3: Pivot the dataframe to wide format
wide_data = aggregated_data.pivot(index='Well', columns='Peptide', values='Area')

# Reset the column names to clean up the MultiIndex (optional)
wide_data.columns.name = None  # Remove the name 'Peptide'

# Display the result
print(wide_data)


# First, remove anything whose signal is lower than the blank?
blank_averages = wide_data.loc['blank']*0.5

# Step 3: Replace values lower than the blank average with NaN
def mask_below_blank_average(row):
    return row.where(row >= blank_averages, np.nan)

data_with_masked_values = wide_data.apply(mask_below_blank_average, axis=1)
data_with_masked_values[data_with_masked_values == 0] = 0
# Display the result
print(data_with_masked_values)



def total_sum_scaling(df):
    """
    Perform Total Sum Scaling (TSS) normalization on a DataFrame of metabolite intensities.
    Rows are metabolites, and columns are samples.
    """
    # Compute the total intensity (sum) for each sample
    sample_totals = df.sum(axis=0)
    
    # Scale each sample to make the total intensity equal across all samples
    normalized_df = df.div(sample_totals, axis=1) * sample_totals.mean()
    
    return normalized_df




pre_df = data_with_masked_values.transpose().drop('blank', axis = 1)
normalized_data = total_sum_scaling(pre_df).transpose()



data_long = normalized_data.reset_index().melt(
    id_vars='Well',  # Keep the 'Well' column as-is
    var_name='Peptide',  # Name for the peptide/metabolite column
    value_name='Area'  # Name for the intensity/area column
)

# Merge the filtered DataFrame with the layout DataFrame based on the 'Well' column
merged_df = data_long.merge(layout, left_on='Well', right_on='well')
# Remove all zero and nan rows?
merged_df = merged_df[merged_df['Area'] > 0]
merged_df['Normalized_Area'] = merged_df['Area']



# Map the custom names to the 'CONCuM' column
merged_df['Experimental group'] = merged_df['CONCuM'].map(custom_group_names)


# Create the grouped barplot
plt.figure(figsize=(12, 6))
barplot = sns.barplot(
    x='Peptide',
    y='Normalized_Area',
    hue='Experimental group',
    data=merged_df,
    ci='sd',
    hue_order=[custom_group_names[g] for g in group_order]
)

# Customize x-axis and y-axis
plt.title('Area by Experimental Group for Each Amino Acid (Normalized by total sum, Log Scale)')
plt.xlabel('Peptide')
plt.ylabel('Mean Area (Log Scale)')
plt.xticks(rotation=45, ha='right')  # Rotate and right-align x-ticks
plt.yscale('log')  # Apply log scale to y-axis

# Prepare pairs for statistical comparison
peptides = merged_df['Peptide'].unique()
pairs = []
for peptide in peptides:
    groups_in_peptide = merged_df[merged_df['Peptide'] == peptide]['Experimental group'].unique()
    group_pairs = list(combinations(groups_in_peptide, 2))
    pairs.extend([((peptide, pair[0]), (peptide, pair[1])) for pair in group_pairs])

# Define p-value thresholds and corresponding symbols
pvalue_thresholds = [
    (0.001, "****"),
    (0.01, "***"),
    (0.05, "**"),
    (0.1, "*")
]

# Initialize the Annotator
annotator = Annotator(
    barplot,
    pairs,
    data=merged_df,
    x='Peptide',
    y='Area',
    hue='Experimental group',
    hue_order=[custom_group_names[g] for g in group_order]
)

# Configure and Apply Annotations
annotator.configure(
    #test='t-test_ind',
    test='Mann-Whitney',
    text_format='star',
    loc='inside',
    alpha = 0.1,
    comparisons_correction=None,
    pvalue_thresholds=pvalue_thresholds,
    hide_non_significant=True,  # Hide non-significant comparisons
    line_offset=0.01,  # Adjust line offset to move lines closer to bars
    text_offset=0.01,  # Adjust text offset to move stars closer to lines
    line_height=0.01,  # Adjust line height
    #line_height_units='axes_fraction'  # Units for offsets ('axes_fraction', 'points', 'data')
)
annotator.apply_and_annotate()

plt.tight_layout()
plt.show()

# Create a DataFrame with the comparison results
from statsmodels.stats.multitest import multipletests

comparison_results = []

for peptide in peptides:
    
    #peptide = 'L-Methionine'
    peptide_data = merged_df[merged_df['Peptide'] == peptide]
    print(peptide)
    print(peptide_data['CONCuM'].value_counts())
    
    groups = peptide_data['CONCuM'].unique()
    group_pairs = list(combinations(groups, 2))
    for pair in group_pairs:
        group1_data = peptide_data[peptide_data['CONCuM'] == pair[0]]['Normalized_Area']
        group2_data = peptide_data[peptide_data['CONCuM'] == pair[1]]['Normalized_Area']
        
        if group1_data.dropna().shape[0] == 0:
            continue
        
        #t_stat, p_value = ttest_ind(group1_data, group2_data, nan_policy='omit')
        t_stat, p_value = mannwhitneyu(group1_data, group2_data, nan_policy='omit')
        
        comparison_results.append({
            'Peptide': peptide,
            'Group1': custom_group_names.get(pair[0], pair[0]),
            'Group2': custom_group_names.get(pair[1], pair[1]),
            'P-Value': p_value
        })

comparison_df = pd.DataFrame(comparison_results)
# Apply BH correction within each peptide
comparison_df['Adjusted P-Value'] = comparison_df.groupby('Peptide')['P-Value'].transform(
    lambda x: multipletests(x, method='fdr_bh')[1]
)








