#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Dec 12 16:57:13 2024

@author: danbru
"""

import pandas as pd
import matplotlib.pyplot as plt
from scipy.stats import ttest_ind
import re


def extract_last_three_digits(filename):
    # Find the pattern starting with "Inj" followed by digits
    match = re.search(r"Inj(\d+)", filename)
    if match:
        # Extract the digits after "Inj" and take the last three digits
        return match.group(1)[-3:]
    return None  # Return None if no match is found

import pandas as pd
import numpy as np


def select_adduct_with_most_nonzero_long_format(data):
    """
    Select the adduct with the most non-zero detections for each metabolite in long format data.

    Parameters:
        data (pd.DataFrame): A DataFrame in long format with columns "Peptide" (metabolite), 
                            "Precursor Mz" (adduct identifier), "Intensity" (intensity values), 
                            and "Sample" (sample identifiers).

    Returns:
        pd.DataFrame: A filtered DataFrame with one adduct per metabolite based on most non-zero detections.
    """
    # Count non-zero intensities for each adduct within each metabolite
    non_zero_counts = (
        data[data['Area'] > 0]
        .groupby(['Peptide', 'Precursor Mz'])
        .size()
        .reset_index(name='NonZeroCount')
    )

    # Find the adduct with the maximum non-zero detections for each metabolite
    best_adducts = (
        non_zero_counts.loc[non_zero_counts.groupby('Peptide')['NonZeroCount'].idxmax()]
    )

    # Filter the original data to include only the selected adducts
    filtered_data = data.merge(best_adducts[['Peptide', 'Precursor Mz']], on=['Peptide', 'Precursor Mz'])

    return best_adducts, filtered_data


import pandas as pd

def select_adduct_most_represented_across_groups(data):
    """
    Select the adduct with the highest mean occurrence across experimental groups for each metabolite.

    Parameters:
        data (pd.DataFrame): A DataFrame in long format with columns "Peptide" (metabolite), 
                             "Precursor Mz" (adduct identifier), "Intensity" (intensity values), 
                             "Sample" (sample identifiers), and "Summary" (experimental group).

    Returns:
        pd.DataFrame: A filtered DataFrame with one adduct per metabolite based on the highest mean occurrence per group.
    """
    # Count non-zero intensities per experimental group
    non_zero_counts = (
        data[data['Area'] > 0]
        .groupby(['Peptide', 'Precursor Mz', 'Summary'])
        .size()
        .reset_index(name='NonZeroCount')
    )
    
    # Compute the mean non-zero detections across groups for each adduct
    mean_non_zero_per_group = (
        non_zero_counts
        .groupby(['Peptide', 'Precursor Mz'])['NonZeroCount']
        .mean()
        .reset_index(name='MeanNonZeroCount')
    )
    
    # Find the adduct with the highest mean detection per group for each metabolite
    best_adducts = mean_non_zero_per_group.loc[
        mean_non_zero_per_group.groupby('Peptide')['MeanNonZeroCount'].idxmax()
    ]
    
    # Filter the original data to include only the selected adducts
    filtered_data = data.merge(best_adducts[['Peptide', 'Precursor Mz']], on=['Peptide', 'Precursor Mz'])
    
    return best_adducts, filtered_data


def select_adduct_with_most_nonzero_long_format_new(data):
    """
    Select the adduct with the most non-zero detections for each metabolite, ensuring presence across all experimental groups in the "Summary" column.

    Parameters:
        data (pd.DataFrame): A DataFrame in long format with columns "Peptide" (metabolite), 
                            "Precursor Mz" (adduct identifier), "Intensity" (intensity values), 
                            "Sample" (sample identifiers), and "Summary" (experimental groups).

    Returns:
        pd.DataFrame: A filtered DataFrame with one adduct per metabolite based on group-complete non-zero detections.
    """
    # Count non-zero intensities for each adduct within each group and peptide
    group_counts = (
        data[data['Area'] > 0]
        .groupby(['Peptide', 'Precursor Mz', 'Summary'])
        .size()
        .reset_index(name='NonZeroCount')
    )

    # Check if an adduct is present in all groups for each peptide
    group_presence = (
        group_counts
        .groupby(['Peptide', 'Precursor Mz'])['Summary']
        .nunique()
        .reset_index(name='GroupCount')
    )
    total_groups = data['Summary'].nunique()
    group_presence = group_presence[group_presence['GroupCount'] == total_groups]

    # Merge with original group counts and calculate total non-zero counts for valid adducts
    valid_adducts = (
        group_counts.merge(group_presence[['Peptide', 'Precursor Mz']], on=['Peptide', 'Precursor Mz'])
        .groupby(['Peptide', 'Precursor Mz'])['NonZeroCount']
        .sum()
        .reset_index()
    )

    # Select the adduct with the maximum non-zero detections for each peptide
    best_adducts = (
        valid_adducts.loc[valid_adducts.groupby('Peptide')['NonZeroCount'].idxmax()]
    )

    # Filter the original data to include only the selected adducts
    filtered_data = data.merge(best_adducts[['Peptide', 'Precursor Mz']], on=['Peptide', 'Precursor Mz'])

    return best_adducts, filtered_data

def probabilistic_quotient_normalization_long_format(data):
    
    """
    Perform Probabilistic Quotient Normalization (PQN) on long format data.

    Parameters:
        data (pd.DataFrame): A DataFrame in long format with columns "Sample", "Peptide", 
                            "Precursor Mz", and "Intensity".

    Returns:
        pd.DataFrame: A DataFrame with normalized intensities.
    """
    
    # Pivot data to wide format (samples as rows, peptides as columns)
    wide_data = data.pivot_table(index='Well', columns='Peptide', values='Area', aggfunc='first')

    # Step 1: Compute the reference spectrum (median across all samples for each peptide)
    reference_spectrum = wide_data.median(axis=0)

    # Step 2: Compute the ratio of each intensity to the reference spectrum
    ratios = wide_data.div(reference_spectrum, axis=1)

    # Step 3: Compute the median ratio (quotient) for each sample
    quotients = ratios.median(axis=1)

    # Step 4: Normalize each sample by its quotient
    normalized_wide_data = wide_data.div(quotients, axis=0)

    # Convert back to long format
    normalized_data = normalized_wide_data.reset_index().melt(id_vars='Well', var_name='Peptide', value_name='NormalizedArea')

    return normalized_data

def probabilistic_quotient_normalization_long_format(data):
    """
    Perform Probabilistic Quotient Normalization (PQN) on long format data.

    Parameters:
        data (pd.DataFrame): A DataFrame in long format with columns "Well" (Sample), 
                            "Peptide", and "Area" (Intensity).

    Returns:
        pd.DataFrame: A DataFrame with normalized intensities.
    """
    
    # Pivot data to wide format (samples as rows, peptides as columns)
    wide_data = data.pivot_table(index='Well', columns='Peptide', values='Area', aggfunc='first')

    # Step 1: Compute the reference spectrum (median across all samples for each peptide)
    reference_spectrum = wide_data.median(axis=0, skipna=True)  # Ensure missing values are ignored

    # Step 2: Compute the ratio of each intensity to the reference spectrum
    ratios = wide_data.div(reference_spectrum, axis=1)

    # Step 3: Compute the median ratio (quotient) for each sample
    quotients = ratios.median(axis=1, skipna=True)

    # Step 4: Prevent division by zero or NaN
    quotients.replace([0, None], 1, inplace=True)  # Ensure no division errors

    # Step 5: Normalize each sample by its quotient
    normalized_wide_data = wide_data.div(quotients, axis=0)

    # Convert back to long format
    normalized_data = normalized_wide_data.reset_index().melt(id_vars='Well', var_name='Peptide', value_name='NormalizedArea')

    return normalized_data


def calculate_p_values(normalized_data):
    """
    Calculate p-values for all pairwise comparisons within each Peptide group.

    Parameters:
        normalized_data (pd.DataFrame): A DataFrame with normalized intensities, "Peptide", and "Summary" columns.

    Returns:
        pd.DataFrame: A DataFrame with p-values for all pairwise comparisons within each Peptide.
    """
    p_values = []

    for peptide, group in normalized_data.groupby('Peptide'):
        summaries = group['Summary'].unique()
        for i, summary1 in enumerate(summaries):
            for summary2 in summaries[i + 1:]:
                group1 = group[group['Summary'] == summary1]['NormalizedArea']
                group2 = group[group['Summary'] == summary2]['NormalizedArea']
                stat, p_val = ttest_ind(group1, group2, equal_var=False, nan_policy='omit')
                p_values.append({
                    'Peptide': peptide,
                    'Group1': summary1,
                    'Group2': summary2,
                    'P-Value': p_val
                })

    return pd.DataFrame(p_values)


def select_adduct_most_represented_across_groups(data):
    """
    Select the adduct with the highest occurrence across experimental groups for each metabolite and provide a summary for all adducts.

    Parameters:
        data (pd.DataFrame): A DataFrame in long format with columns "Peptide" (metabolite), 
                             "Precursor Mz" (adduct identifier), "Intensity" (intensity values), 
                             "Sample" (sample identifiers), "Replicate" (well identifier),
                             and "Summary" (experimental group).

    Returns:
        tuple: A tuple containing:
            - pd.DataFrame: A filtered DataFrame with one adduct per metabolite based on the highest occurrence per group.
            - pd.DataFrame: A pivoted DataFrame with adducts as rows and experimental groups as columns, showing occurrence counts.
            - pd.DataFrame: A summary DataFrame containing all adducts with replicate counts per experimental group.
    """
    # Count unique replicates with non-zero intensities per experimental group
    unique_replicate_counts = (
        data[data['Area'] > 0]
        .groupby(['Peptide', 'Precursor Mz', 'Summary'])['Replicate']
        .nunique()
        .reset_index(name='UniqueReplicateCount')
    )
    
    # Compute the total unique replicate counts across groups for each adduct
    total_unique_replicate_per_group = (
        unique_replicate_counts
        .groupby(['Peptide', 'Precursor Mz'])['UniqueReplicateCount']
        .sum()
        .reset_index(name='TotalUniqueReplicateCount')
    )
    
    # Find the adduct with the highest unique replicate detection per group for each metabolite
    best_adducts = total_unique_replicate_per_group.loc[
        total_unique_replicate_per_group.groupby('Peptide')['TotalUniqueReplicateCount'].idxmax()
    ]
    
    # Filter the original data to include only the selected adducts
    filtered_data = data.merge(best_adducts[['Peptide', 'Precursor Mz']], on=['Peptide', 'Precursor Mz'])
    
    # Create a pivot table to show unique replicate counts per experimental group
    adduct_group_counts = unique_replicate_counts.pivot_table(
        index=['Peptide', 'Precursor Mz'],
        columns='Summary',
        values='UniqueReplicateCount',
        fill_value=0
    ).reset_index()
    
    # Create a summary DataFrame for all adducts with unique replicate counts
    all_adducts_summary = adduct_group_counts.sort_values(by=['Peptide'], ascending=[True])
    
    return best_adducts, filtered_data, adduct_group_counts, all_adducts_summary


def select_adduct_most_represented_across_groups(data, important_groups=None):
    """
    Select the adduct with the highest occurrence across specified experimental groups for each metabolite and provide a summary for all adducts.

    Parameters:
        data (pd.DataFrame): A DataFrame in long format with columns "Peptide" (metabolite), 
                             "Precursor Mz" (adduct identifier), "Intensity" (intensity values), 
                             "Sample" (sample identifiers), "Replicate" (well identifier),
                             and "Summary" (experimental group).
        important_groups (list, optional): A list of experimental groups to prioritize when selecting adducts.

    Returns:
        tuple: A tuple containing:
            - pd.DataFrame: A filtered DataFrame with one adduct per metabolite based on the highest occurrence in the specified groups.
            - pd.DataFrame: A pivoted DataFrame with adducts as rows and experimental groups as columns, showing occurrence counts.
            - pd.DataFrame: A summary DataFrame containing all adducts with replicate counts per experimental group.
    """
    # Count unique replicates with non-zero intensities per experimental group
    unique_replicate_counts = (
        data[data['Area'] > 0]
        .groupby(['Peptide', 'Precursor Mz', 'Summary'])['Replicate']
        .nunique()
        .reset_index(name='UniqueReplicateCount')
    )
    
    # If important groups are specified, filter only for those groups for adduct selection
    if important_groups is not None:
        filtered_selection = unique_replicate_counts[unique_replicate_counts['Summary'].isin(important_groups)]
    else:
        filtered_selection = unique_replicate_counts
    
    # Compute the total unique replicate counts across specified groups for each adduct
    total_unique_replicate_per_group = (
        filtered_selection
        .groupby(['Peptide', 'Precursor Mz'])['UniqueReplicateCount']
        .sum()
        .reset_index(name='TotalUniqueReplicateCount')
    )
    
    # Find the adduct with the highest unique replicate detection per group for each metabolite
    best_adducts = total_unique_replicate_per_group.loc[
        total_unique_replicate_per_group.groupby('Peptide')['TotalUniqueReplicateCount'].idxmax()
    ]
    
    # Filter the original data to include only the selected adducts
    filtered_data = data.merge(best_adducts[['Peptide', 'Precursor Mz']], on=['Peptide', 'Precursor Mz'])
    
    # Create a pivot table to show unique replicate counts per experimental group for all adducts
    adduct_group_counts = unique_replicate_counts.pivot_table(
        index=['Peptide', 'Precursor Mz'],
        columns='Summary',
        values='UniqueReplicateCount',
        fill_value=0
    ).reset_index()
    
    # Create a summary DataFrame for all adducts with unique replicate counts
    all_adducts_summary = adduct_group_counts.sort_values(by=['Peptide'], ascending=[True])
    
    return best_adducts, filtered_data, adduct_group_counts, all_adducts_summary




exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/lactic acid/run/proline_0.5_5_20241219_1102'

# Compare the active wells?
layout = pd.read_csv(f'{exp_path}/protocol/plate_layout/layout.tsv', sep = '\t')
layout['well'] = layout['well'].str.replace(r'(\D)0+(\d+)', r'\1\2', regex=True)
layout = layout[layout['CONCuM'] != 'Media control']

layout.loc[layout['Unnamed: 0'].max()+1] = [layout['Unnamed: 0'].max(), 'plate_1','Blank','Blank','Solvent Blank']  # adding a row

# Add a row to the layout, the blank


# Only keep the relevant wells
#layout = layout[layout['CONCuM'].isin(['L-glutamate: 0 mM & Formic acid: 0 mM', 'L-glutamate: 5 mM & Formic acid: 0 mM', 'L-glutamate: 15 mM & Formic acid: 0 mM'])]
group_order = layout['CONCuM'].unique()  # Replace with your desired order
#custom_group_names = {
#    "L-glutamate: 0 mM & Formic acid: 0 mM": "No Glutamate",
#    "L-glutamate: 5 mM & Formic acid: 0 mM": "Low Glutamate",
#    "L-glutamate: 15 mM & Formic acid: 0 mM": "High Glutamate"
#}


#116.071 for proline, 138.053 for control

data = pd.read_csv(f'{exp_path}/results/metabolomics/TransitionResults.csv')
data['Well'] = data['Replicate'].str.extract(r'-([A-H]\d{1,2})')
#data['Peptide'] = data['Peptide'].str.replace('L-Glutamic acid','L-Glutamic Acid')

# Add the blanks as well
# Change values in 'Well' column to 'Blank' if 'Replicate' contains 'MAT1'
data.loc[data['Replicate'].str.contains('MAT1', na=False), 'Well'] = 'Blank'


#data = data[data['Protein'].str.contains('AllCCS')]

data = data.merge(layout, left_on = 'Well', right_on = 'well')
data = data.drop_duplicates(subset = ['Peptide','Precursor Mz', 'Well'], keep = 'first')

best_adducts, filtered_data, adduct_group_counts, all_adducts_summary = select_adduct_most_represented_across_groups(data.dropna(), important_groups = ['Control group with no supplementation or treatment.',
                                                                                                                                                        'Low proline supplementation without lactic acid treatment.',
                                                                                                                                                        'High proline supplementation without lactic acid treatment.'])

best_adducts, filtered_data, adduct_group_counts, all_adducts_summary = select_adduct_most_represented_across_groups(data.dropna(), important_groups = ['Control group with no supplementation or treatment.'])
filtered_test_data = filtered_data[filtered_data['Peptide'].isin(['Proline'])]


filtered_data = filtered_data[filtered_data['Area'] > 0]
normalized_data = probabilistic_quotient_normalization_long_format(filtered_data).dropna().drop_duplicates()
normalized_data = normalized_data.merge(data[['Well','Summary']].drop_duplicates(), left_on = 'Well', right_on = 'Well') 

plot_data = normalized_data[['Well','Peptide','NormalizedArea','Summary']]


#plot_data = plot_data[plot_data['Peptide'].isin(['Glutamic acid','Proline','L-Arginine','L-Cysteine','L-Tryptophan',
#                                                 'L-Valine','Phenylalanine','Glutamine','Glutathione'])]
plot_data = plot_data[plot_data['Peptide'].isin(['Glutamic acid', 'L-Arginine', 'Proline'])]
summary_order = []

# Group by Peptide and Summary and calculate the mean normalized intensity
# Group by Peptide and Summary and calculate the mean normalized intensity
grouped_data = plot_data.groupby(['Peptide', 'Summary'])['NormalizedArea'].mean().reset_index()

# Create a complete set of Peptide and Summary combinations to ensure alignment
all_peptides = grouped_data['Peptide'].unique()
all_summaries = summary_order if summary_order else grouped_data['Summary'].unique()
complete_index = pd.MultiIndex.from_product([all_peptides, all_summaries], names=['Peptide', 'Summary'])
grouped_data = grouped_data.set_index(['Peptide', 'Summary']).reindex(complete_index).reset_index()

# Apply specific order for the Summary column if provided
if summary_order:
    grouped_data['Summary'] = pd.Categorical(grouped_data['Summary'], categories=summary_order, ordered=True)
    grouped_data = grouped_data.sort_values(['Peptide', 'Summary'])

# Create a bar plot
plt.figure(figsize=(12, 6))
peptides = grouped_data['Peptide'].unique()
x = np.arange(len(peptides))  # the label locations
width = 0.8 / len(all_summaries)  # the width of the bars

# Separate data by Summary groups
for i, summary in enumerate(all_summaries):
    summary_data = grouped_data[grouped_data['Summary'] == summary]
    plt.bar(x + (i - len(all_summaries) / 2) * width, summary_data['NormalizedArea'], width=width, label=summary)

# Add labels, title, and legend
plt.xlabel('Peptide')
plt.ylabel('Mean Normalized Area')
plt.title('Normalized Areas Grouped by Peptide and Summary')
plt.xticks(x, peptides, rotation=45, ha='right')
plt.legend(title='Summary')
plt.yscale('log')
plt.tight_layout()
plt.show()


# Calculate p values
pval_df = calculate_p_values(plot_data)
