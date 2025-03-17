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
    wide_data = data.pivot_table(index='Well', columns='Molecule Name', values='Area', aggfunc='first')

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
    normalized_data = normalized_wide_data.reset_index().melt(id_vars='Well', var_name='Molecule Name', value_name='NormalizedArea')

    return normalized_data

def qc_based_pqn_normalization(data):
    """
    Perform QC-driven Probabilistic Quotient Normalization (PQN) on long format data.
    
    Parameters:
        data (pd.DataFrame): DataFrame in long format with columns "Well", 
                             "Molecule Name", and "Area".
    
    Returns:
        pd.DataFrame: DataFrame with normalized intensities.
    """
    # Separate QC and study samples
    qc_data = data[data['Well'] == 'QC']
    study_data = data[data['Well'] != 'QC']
    
    # Pivot QC data to wide format (QCs as rows, peptides as columns)
    qc_wide = qc_data.pivot_table(index='Well', columns='Molecule Name', 
                                  values='Area', aggfunc='first')
    
    # Compute the QC reference spectrum (median across QC samples for each peptide)
    #qc_reference = qc_wide.median(axis=0, skipna=True)
    
    # Compute the QC reference spectrum (median across QC samples for each peptide)
    qc_reference = qc_wide.median(axis=0, skipna=True)
    # Optionally fill NaNs in reference with 1 to avoid dividing by NaN
    #qc_reference.fillna(1, inplace=True)
    
    # Pivot study data to wide format (samples as rows, peptides as columns)
    study_wide = study_data.pivot_table(index='Well', columns='Molecule Name', 
                                        values='Area', aggfunc='first')
    
    # Compute ratios of each study sample's intensity to the QC reference spectrum
    ratios = study_wide.div(qc_reference, axis=1)
    
    # Compute the median ratio (quotient) for each study sample
    quotients = ratios.median(axis=1, skipna=True)
    
    # Prevent division errors: replace zeros or NaNs with 1
    quotients.replace([0, None], 1, inplace=True)
    
    # Normalize each study sample by its quotient
    normalized_wide = study_wide.div(quotients, axis=0)
    
    # Convert back to long format
    normalized_data = normalized_wide.reset_index().melt(id_vars='Well', 
                                                         var_name='Molecule Name', 
                                                         value_name='NormalizedArea')
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

    for peptide, group in normalized_data.groupby('Molecule Name'):
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

def select_adduct_study_samples(data, restrict_with_blank=True, min_count=1):
    """
    Select the adduct with the highest occurrence across study samples (excluding QC and Blank) 
    for each metabolite, while optionally applying a restriction that only intensities above 
    the average blank value count. The returned filtered data will include all rows (study, QC, 
    and Blank) for the selected adducts.
    
    Parameters:
        data (pd.DataFrame): DataFrame in long format with columns:
            - "Molecule Name" (metabolite)
            - "Precursor Mz" (adduct identifier)
            - "Area" (intensity)
            - "Well" (sample type, e.g., study sample, QC, or Blank)
            - "Replicate" (replicate or well identifier)
        restrict_with_blank (bool, optional): If True, for each adduct the average blank intensity
            is computed and only study sample measurements that exceed this blank average are counted.
        min_count (int, optional): Minimum number of valid replicates required for an adduct to be considered.
    
    Returns:
        tuple: A tuple containing:
            - best_adducts (pd.DataFrame): One row per metabolite with the selected adduct.
            - filtered_data (pd.DataFrame): Data filtered to include all rows (study, QC, and Blank) 
              corresponding to the selected best adduct per metabolite.
            - adduct_counts (pd.DataFrame): A pivot table with counts of valid replicates per 
              metabolite and adduct (from study samples).
    """
    import pandas as pd

    # 1. Filter for study samples only (for selection purposes)
    study_data = data[~data['Well'].isin(['QC', 'Blank'])].copy()

    # 2. If restricting with blanks, compute blank averages per adduct.
    if restrict_with_blank:
        blank_data = data[data['Well'] == 'Blank'].fillna(0)
        # Compute the average blank intensity for each (Molecule Name, Precursor Mz)
        blank_avg = blank_data.groupby(['Molecule Name', 'Precursor Mz'])['Area'] \
                              .mean() \
                              .reset_index() \
                              .rename(columns={'Area': 'BlankAvg'})
        # Merge blank average with study data on Molecule Name and Precursor Mz
        study_data = study_data.merge(blank_avg, on=['Molecule Name', 'Precursor Mz'], how='left').fillna(0)
        # For study samples, only count measurements that exceed the blank average (if available)
        study_data = study_data[
            (study_data['BlankAvg'].isna()) | (study_data['Area'] > study_data['BlankAvg'])
        ]
    
    # 3. Count unique replicates with non-zero intensities for each adduct in study samples
    valid_counts = (
        study_data[study_data['Area'] > 0]
        .groupby(['Molecule Name', 'Precursor Mz'])['Replicate Name']
        .nunique()
        .reset_index(name='UniqueReplicateCount')
    )
    
    # Filter out adducts with counts below a minimum threshold.
    valid_counts = valid_counts[valid_counts['UniqueReplicateCount'] >= min_count]
    
    # 4. Create a pivot table summarizing counts for all adducts (based on study samples)
    adduct_counts = valid_counts.pivot_table(
        index=['Molecule Name', 'Precursor Mz'],
        values='UniqueReplicateCount',
        fill_value=0
    ).reset_index()
    
    # 5. For each metabolite, select the adduct with the highest count.
    best_adducts = valid_counts.loc[
        valid_counts.groupby('Molecule Name')['UniqueReplicateCount'].idxmax()
    ]
    
    # 6. Now, filter the ORIGINAL data to include all rows (study, QC, and Blank) 
    # for the selected best adducts.
    filtered_data = data.merge(
        best_adducts[['Molecule Name', 'Precursor Mz']],
        on=['Molecule Name', 'Precursor Mz'],
        how='inner'
    )
    
    return best_adducts, filtered_data, adduct_counts

def select_adduct_study_samples(
    data,
    restrict_with_blank=True,
    min_count=1,
    mh_weight_factor=1.25
):
    """
    Select the adduct with the highest occurrence across study samples (excluding QC and Blank) 
    for each metabolite, while optionally applying a restriction that only intensities above 
    the average blank value count. The returned filtered data will include all rows (study, QC, 
    and Blank) for the selected adducts.

    Additionally, M+H adducts are slightly "favored" by multiplying their replicate count
    by a specified mh_weight_factor (default 1.1).

    Parameters:
    -----------
        data (pd.DataFrame): DataFrame in long format with columns:
            - "Molecule Name" (metabolite)
            - "Precursor Mz" (adduct identifier, e.g. "M+H", "M+Na", etc.)
            - "Area" (intensity)
            - "Well" (sample type, e.g., 'QC', 'Blank', or a study sample)
            - "Replicate" (replicate or well identifier)
        restrict_with_blank (bool, optional): If True, for each adduct the average blank intensity
            is computed and only study sample measurements that exceed this blank average are counted.
        min_count (int, optional): Minimum number of valid replicates required for an adduct to be considered.
        mh_weight_factor (float, optional): Factor by which to multiply replicate counts for "M+H" adducts.
            e.g. 1.1 => 10% bonus for M+H coverage.

    Returns:
    --------
        tuple of (best_adducts, filtered_data, adduct_counts):
            best_adducts : pd.DataFrame
                One row per metabolite with the selected adduct (including WeightedCount).
            filtered_data : pd.DataFrame
                Rows (study, QC, Blank) corresponding to these selected best adducts.
            adduct_counts : pd.DataFrame
                Pivot table with counts of valid replicates (UniqueReplicateCount) per 
                (Molecule, Precursor Mz) (based on study samples).
    """
    import pandas as pd

    # 1. Filter for study samples only (for selection purposes)
    study_data = data[~data['Well'].isin(['QC', 'Blank'])].copy()

    # 2. If restricting with blanks, compute blank averages per adduct.
    if restrict_with_blank:
        blank_data = data[data['Well'] == 'Blank'].fillna(0)
        # Compute the average blank intensity for each (Molecule Name, Precursor Mz)
        blank_avg = (
            blank_data
            .groupby(['Molecule Name', 'Precursor Mz'])['Area']
            .mean()
            .reset_index()
            .rename(columns={'Area': 'BlankAvg'})
        )
        # Merge blank average with study data
        study_data = (
            study_data
            .merge(blank_avg, on=['Molecule Name', 'Precursor Mz'], how='left')
            .fillna(0)
        )
        # Only count measurements that exceed blank average (if available)
        study_data = study_data[
            (study_data['BlankAvg'].isna()) | (study_data['Area'] > study_data['BlankAvg'])
        ]
    
    # 3. Count unique replicates with non-zero intensities for each adduct in study samples
    valid_counts = (
        study_data[study_data['Area'] > 0]
        .groupby(['Molecule Name', 'Precursor Mz'])['Replicate Name']
        .nunique()
        .reset_index(name='UniqueReplicateCount')
    )

    # 4. Filter out adducts with counts below min_count
    valid_counts = valid_counts[valid_counts['UniqueReplicateCount'] >= min_count]

    # 5. Create a pivot table summarizing coverage for all adducts
    adduct_counts = valid_counts.pivot_table(
        index=['Molecule Name', 'Precursor Mz'],
        values='UniqueReplicateCount',
        fill_value=0
    ).reset_index()

    # -------- New Part: Weight M+H Adducts Slightly Higher -------------
    def apply_mh_weight(row):
        # If "Precursor Mz" exactly matches "M+H" (case-sensitive),
        # multiply coverage by mh_weight_factor. Otherwise, leave as is.
        if row['Precursor Adduct'].contains('M+H'):
            return row['UniqueReplicateCount'] * mh_weight_factor
        else:
            return row['UniqueReplicateCount']

    valid_counts['WeightedCount'] = valid_counts.apply(apply_mh_weight, axis=1)

    # 6. For each metabolite, select the adduct with the highest WeightedCount
    best_adducts = valid_counts.loc[
        valid_counts.groupby('Molecule Name')['WeightedCount'].idxmax()
    ]

    # 7. Now, filter the ORIGINAL data to include all rows for these best adducts
    filtered_data = data.merge(
        best_adducts[['Molecule Name', 'Precursor Mz']],
        on=['Molecule Name', 'Precursor Mz'],
        how='inner'
    )

    return best_adducts, filtered_data, adduct_counts

def impute_lod_half_per_peptide_if_under_50pct(df, threshold=0.5):
    """
    For each Peptide group:
      1. Check the fraction of values that are 0 or NaN.
      2. If < treshold, compute LOD (minimum positive value) and impute 0/NaN with LOD/2.
      3. If >= treshold, leave the values as-is.
      
    Parameters:
      df (pd.DataFrame): Must have columns ["Peptide", "Area"] (and possibly others).
      
    Returns:
      pd.DataFrame: A copy of df with imputed values where appropriate.
    """
    def fill_group(grp):
        # Calculate fraction of zero or NaN
        fraction_zero_na = ((grp['Area'] == 0) | (grp['Area'].isna())).mean()
        
        if fraction_zero_na < threshold:
            # Compute LOD as the smallest positive intensity
            lod = grp.loc[grp['Area'] > 0, 'Area'].min()
            # Define fallback to avoid errors if no positive values exist
            lod_val = lod / 2 if pd.notna(lod) else 1e-6
            
            # Replace 0 with NaN, then fill
            grp['Area'] = grp['Area'].replace(0, np.nan)
            grp['Area'] = grp['Area'].fillna(lod_val)
        
        # If fraction_zero_na >= 0.5, do nothing (no imputation)
        return grp

    # Group by Peptide and apply the function
    df_imputed = df.groupby('Molecule Name', group_keys=False).apply(fill_group)
    return df_imputed



def qc_based_pqn_normalization(data, min_peptide_valid_fraction=0.5, min_valid_ratios=10):
    """
    Perform QC-driven Probabilistic Quotient Normalization (PQN) on long format data,
    handling missing values by filtering or imputation.
    
    Parameters:
        data (pd.DataFrame): DataFrame in long format with columns "Well", 
                             "Molecule Name", and "Area".
        min_peptide_valid_fraction (float): Minimum fraction of samples in which a peptide 
                                              must be observed to be included.
        min_valid_ratios (int): Minimum number of valid peptide ratios required for a sample.
    
    Returns:
        pd.DataFrame: DataFrame with normalized intensities.
    """
    # Separate QC and study samples
    qc_data = data[data['Well'] == 'QC']
    study_data = data[data['Well'] != 'QC']
    
    # Determine common peptides between QC and study data
    common_peptides = set(qc_data['Molecule Name']).intersection(study_data['Molecule Name'])
    if not common_peptides:
        raise ValueError("No common peptides found between QC and study samples.")
    
    qc_data = qc_data[qc_data['Molecule Name'].isin(common_peptides)]
    study_data = study_data[study_data['Molecule Name'].isin(common_peptides)]
    
    # Pivot QC data to wide format
    qc_wide = qc_data.pivot_table(index='Well', columns='Molecule Name', 
                                  values='Area', aggfunc='first')
    
    # Compute the QC reference spectrum (median across QC samples for each peptide)
    qc_reference = qc_wide.median(axis=0, skipna=True)
    qc_reference = qc_reference.dropna()  # Remove peptides with no QC signal
    
    # Pivot study data to wide format and restrict to peptides with a valid QC reference
    study_wide = study_data.pivot_table(index='Well', columns='Molecule Name', 
                                        values='Area', aggfunc='first')
    study_wide = study_wide[qc_reference.index]
    
    # Filter out peptides with insufficient valid values in study samples
    min_valid = min_peptide_valid_fraction * study_wide.shape[0]
    valid_peptides = study_wide.columns[study_wide.notna().sum() >= min_valid]
    study_wide = study_wide[valid_peptides]
    qc_reference = qc_reference[qc_reference.index.isin(valid_peptides)]
    
    # Optionally, impute missing values with the median of each peptide
    study_wide = study_wide.apply(lambda x: x.fillna(x.median()), axis=0)
    
    # Compute ratios of each study sample's intensity to the QC reference spectrum
    ratios = study_wide.div(qc_reference, axis=1)
    
    # Check that each sample has enough valid ratios
    valid_counts = ratios.notna().sum(axis=1)
    insufficient_samples = valid_counts[valid_counts < min_valid_ratios].index
    if not insufficient_samples.empty:
        print("Warning: The following samples have too few valid peptide ratios and may be unreliable:", insufficient_samples)
        # Optionally, remove these samples:
        ratios = ratios.drop(index=insufficient_samples)
        study_wide = study_wide.drop(index=insufficient_samples)
    
    # Compute the median ratio (quotient) for each sample, with skipna=True
    quotients = ratios.median(axis=1, skipna=True)
    quotients = quotients.replace(0, 1).fillna(1)
    
    # Normalize each sample by its quotient
    normalized_wide = study_wide.div(quotients, axis=0)
    
    # Convert back to long format
    normalized_data = normalized_wide.reset_index().melt(id_vars='Well', 
                                                         var_name='Molecule Name', 
                                                         value_name='NormalizedArea')
    return normalized_data

def detect_outlier_wells_by_group(data, group_col='Summary', metric='TIC', threshold=3):
    """
    Detect outlier wells based on a summary metric within each experimental group.
    
    Parameters:
        data (pd.DataFrame): DataFrame with at least the columns 'Well', 'Area', 
                             and an experimental group column (default 'Summary').
        group_col (str): Column name for the experimental group.
        metric (str): Metric to use ('TIC' or 'median').
                      'TIC' uses total ion current, while 'median' uses the median intensity per well.
        threshold (float): Threshold in terms of MAD to flag outliers.
    
    Returns:
        tuple: A tuple containing:
            - outliers_by_group (dict): Dictionary with experimental group as keys and lists 
              of outlier wells as values.
            - all_outliers (list): A flat list of all unique wells flagged as outliers.
    """
    outliers_by_group = {}
    all_outliers = set()
    
    # Process each group separately
    for grp, grp_data in data.groupby(group_col):
        if metric == 'TIC':
            summary = grp_data.groupby('Well')['Area'].sum()
        elif metric == 'median':
            summary = grp_data.groupby('Well')['Area'].median()
        else:
            raise ValueError("Metric must be 'TIC' or 'median'")
        
        med = summary.median()
        mad = np.median(np.abs(summary - med))
        
        # Flag wells deviating more than threshold * MAD from the median
        outliers = summary[np.abs(summary - med) > threshold * mad].index.tolist()
        outliers_by_group[grp] = outliers
        
        # Add to the overall set of outliers
        all_outliers.update(outliers)
    
    return outliers_by_group, list(all_outliers)


def summarize_all_adducts_no_scoring(data):
    """
    Summarize ALL (Molecule, Adduct) pairs with these columns:
      1) Molecule Name
      2) Precursor Adduct
      3) QCCV (coefficient of variation in QC)
      4) Coverage (in total) (# of study replicates where Area>0)
      5) Mean intensity in study samples (including zeros)
      6) Mean intensity in blank (including zeros)

    This function does NOT filter out or score any adducts. 
    You can apply thresholds or scoring later if desired.

    Parameters
    ----------
    data : pd.DataFrame
        Must have columns:
          ["Molecule Name", "Adduct", "Area", "Well", "Replicate"] (or "Replicate Name").
        'Well' should allow us to identify "Blank" vs "QC" vs "Study" samples.

    Returns
    -------
    pd.DataFrame with columns:
      [
        "Molecule Name", 
        "Precursor Adduct", 
        "QCCV", 
        "Coverage (in total)",
        "Mean intensity in study samples",
        "Mean intensity in blank"
      ]
    """

    df = data.copy()

    # 1. Identify sample type
    #    Adjust if your well naming is different (case-insensitive for 'blank' or 'qc').
    df['SampleType'] = 'Study'
    df.loc[df['Well'].str.lower().str.contains('blank'), 'SampleType'] = 'Blank'
    df.loc[df['Well'].str.lower().str.contains('qc'), 'SampleType'] = 'QC'

    group_cols = ["Molecule Name", "Precursor Adduct"]

    # A) Calculate QCCV = stdev / mean for QC replicates
    qc_stats = (
        df[df['SampleType'] == 'QC']
        .groupby(group_cols)['Area']
        .agg(['mean', 'std'])
        .rename(columns={'mean': 'MeanQC', 'std': 'StdQC'})
        .reset_index()
    )
    qc_stats['QCCV'] = qc_stats['StdQC'] / qc_stats['MeanQC']
    qc_stats.drop(columns=['MeanQC', 'StdQC'], inplace=True, errors='ignore')

    # B) Coverage in study samples = # of replicates with Area > 0
    coverage_stats = (
        df[df['SampleType'] == 'Study']
        .query("Area > 0")
        .groupby(group_cols)['Replicate Name']
        .nunique()
        .rename("Coverage")
        .reset_index()
    )

    # C) Mean intensity in study samples (including zeros)
    mean_study = (
        df[df['SampleType'] == 'Study'].fillna(0)
        .groupby(group_cols)['Area']
        .mean()  # includes zeros if those replicates exist in the dataset
        .rename("MeanStudy")
        .reset_index()
    )

    # D) Mean intensity in blanks (including zeros)
    mean_blank = (
        df[df['SampleType'] == 'Blank'].fillna(0)
        .groupby(group_cols)['Area']
        .mean()
        .rename("MeanBlank")
        .reset_index()
    )

    # E) Merge all summaries together
    merged = (
        qc_stats
        .merge(coverage_stats, on=group_cols, how='outer')
        .merge(mean_study,     on=group_cols, how='outer')
        .merge(mean_blank,     on=group_cols, how='outer')
    )

    # Fill missing values: 
    # - If no QC data => QCCV = NaN
    # - If coverage is missing => 0
    # - If no study/blank => 0 for means
    merged['QCCV']      = merged['QCCV'].fillna(np.nan)
    merged['Coverage']  = merged['Coverage'].fillna(0)
    merged['MeanStudy'] = merged['MeanStudy'].fillna(0)
    merged['MeanBlank'] = merged['MeanBlank'].fillna(0)

    # Rename "Adduct" -> "Precursor Adduct" to match requested column name
    #merged.rename(columns={"Adduct": "Precursor Adduct"}, inplace=True)

    # 2. Reorder & rename columns for final output
    final_cols = [
        "Molecule Name",
        "Precursor Adduct",
        "QCCV",
        "Coverage",
        "MeanStudy",
        "MeanBlank"
    ]
    final_df = merged[final_cols].copy()
    final_df.rename(
        columns={
            "Coverage": "Coverage (in total)",
            "MeanStudy": "Mean intensity in study samples",
            "MeanBlank": "Mean intensity in blank",
        },
        inplace=True
    )

    return final_df


def select_adduct_study_samples(
    data,
    restrict_with_blank=True,
    min_count=1,
    mh_weight_factor=1
):
    """
    Select the adduct with the highest occurrence across study samples (excluding QC and Blank) 
    for each metabolite, while optionally applying a restriction that only intensities above 
    the average blank value count.

    M+H is given a slight "boost" (via mh_weight_factor) when choosing the best adduct.

    Parameters:
    -----------
        data (pd.DataFrame): DataFrame in long format with columns:
            - "Molecule Name" (metabolite)
            - "Precursor Adduct" (adduct identifier, e.g., "M+H", "M+Na", etc.)
            - "Area" (intensity)
            - "Well" (sample type, e.g., 'QC', 'Blank', or a study sample)
            - "Replicate" (replicate or well identifier)

        restrict_with_blank (bool, optional): If True, for each adduct the average blank intensity
            is computed and only study sample measurements that exceed this blank average are counted.

        min_count (int, optional): Minimum number of valid replicates required for an adduct 
            to be considered.

        mh_weight_factor (float, optional): Factor by which to multiply replicate counts for "M+H" 
            adducts, giving "M+H" a slight advantage in the final selection.

    Returns:
    --------
        tuple of (best_adducts, filtered_data, adduct_counts):
            best_adducts : pd.DataFrame
                One row per metabolite with the chosen adduct (including WeightedCount).
            filtered_data : pd.DataFrame
                All rows (study, QC, Blank) corresponding to these selected best adducts.
            adduct_counts : pd.DataFrame
                Pivot table with the raw replicate counts (UniqueReplicateCount) per 
                (Molecule, Precursor Adduct).
    """
    import pandas as pd

    # 1. Filter for study samples only (for selection)
    study_data = data[~data['Well'].isin(['QC', 'Blank'])].copy()

    # 2. If restricting with blanks, compute average blank intensities per adduct
    if restrict_with_blank:
        blank_data = data[data['Well'] == 'Blank'].fillna(0)
        blank_avg = (
            blank_data
            .groupby(['Molecule Name', 'Precursor Adduct'])['Area']
            .mean()
            .reset_index()
            .rename(columns={'Area': 'BlankAvg'})
        )
        # Merge blank average into study data
        study_data = (
            study_data
            .merge(blank_avg, on=['Molecule Name', 'Precursor Adduct'], how='left')
            .fillna(0)
        )
        # Only count measurements exceeding the blank average
        study_data = study_data[
            (study_data['BlankAvg'].isna()) | (study_data['Area'] > study_data['BlankAvg'])
        ]

    # 3. Count unique replicates with non-zero intensities for each adduct in study samples
    valid_counts = (
        study_data[study_data['Area'] > 0]
        .groupby(['Molecule Name', 'Precursor Adduct'])['Replicate Name']
        .nunique()
        .reset_index(name='UniqueReplicateCount')
    )

    # 4. Filter out adducts below the min_count threshold
    valid_counts = valid_counts[valid_counts['UniqueReplicateCount'] >= min_count]

    # 5. Create pivot of raw coverage (for user reference)
    adduct_counts = (
        valid_counts
        .pivot_table(
            index=['Molecule Name', 'Precursor Adduct'],
            values='UniqueReplicateCount',
            fill_value=0
        )
        .reset_index()
    )

    # ----- Weighting M+H Adducts -----
    def apply_mh_weight(row):
        if 'M+H' in row['Precursor Adduct']:
            return row['UniqueReplicateCount'] * mh_weight_factor
        else:
            return row['UniqueReplicateCount']

    valid_counts['WeightedCount'] = valid_counts.apply(apply_mh_weight, axis=1)

    # 6. Select the adduct with the highest WeightedCount per metabolite
    best_adducts = valid_counts.loc[
        valid_counts.groupby('Molecule Name')['WeightedCount'].idxmax()
    ]

    # 7. Filter original data to include the best adducts (study, QC, and Blank)
    filtered_data = data.merge(
        best_adducts[['Molecule Name', 'Precursor Adduct']],
        on=['Molecule Name', 'Precursor Adduct'],
        how='inner'
    )

    return best_adducts, filtered_data, adduct_counts

import pandas as pd
import numpy as np
import statsmodels.api as sm

import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.interpolate import interp1d
import pandas as pd
import numpy as np
import statsmodels.api as sm
from scipy.interpolate import interp1d

def qc_rlsc_loess_per_adduct(
    df,
    replicate_col="Replicate Name",
    molecule_col="Molecule Name",
    adduct_col="Precursor Adduct",
    intensity_col="Area",
    sampletype_col="Well",
    qc_label="QC",
    frac=0.3
):
    """
    1) Sorts the DataFrame alphabetically by replicate_col to infer injection order.
    2) For each (molecule, adduct) combination, fits a LOESS curve to QC intensities 
       vs. injection order (if >=3 QC data points).
    3) Corrects the entire group's intensities by dividing by the LOESS-fitted drift,
       scaled so QCs have a median intensity equal to their LOESS median.

    Parameters
    ----------
    df : pd.DataFrame
        Must have at least the following columns (names adjustable via parameters):
          - replicate_col    (e.g., "Replicate Name")
          - molecule_col     (e.g., "Molecule Name")
          - adduct_col       (e.g., "Precursor Adduct")
          - intensity_col    (e.g., "Area")
          - sampletype_col   (e.g., "Well" with values like "QC", "Blank", or sample IDs)
        Will also contain zeros and possibly NaNs in 'intensity_col'.

    replicate_col : str
        The column that identifies each injection/vial. Alphabetical sort = run order.
    
    molecule_col : str
        The column for the metabolite or compound name (e.g., "Molecule Name").
    
    adduct_col : str
        The column for the adduct identifier (e.g., "M+H", "M+Na", etc.).

    intensity_col : str
        The column containing the raw intensity/area measurements (many zeros/NaNs possible).
    
    sampletype_col : str
        The column indicating sample type. We look for rows where sampletype_col == qc_label for QCs.

    qc_label : str
        The value in 'sampletype_col' indicating a QC injection.

    frac : float
        Fraction of data in LOESS window (0 < frac <= 1). 
        Larger => smoother fit, smaller => more local.

    Returns
    -------
    pd.DataFrame
        A copy of the input with:
          - "InjectionOrder" (1-based integer)
          - "Loess_Fit" (the predicted drift for each row)
          - "Intensity_Corrected" (LOESS-corrected intensities)
    """

    # ---- 1. Sort by replicate name to infer run order ----
    df = df.copy()
    df.sort_values(by=replicate_col, inplace=True)
    df["InjectionOrder"] = range(1, len(df) + 1)

    # Prepare new columns for the output
    df["Loess_Fit"] = np.nan
    df["Intensity_Corrected"] = np.nan

    # We'll group by the combination: (molecule_col, adduct_col)
    group_cols = [molecule_col, adduct_col]
    grouped = df.groupby(group_cols)

    corrected_subdfs = []

    for (molecule, adduct), subdf in grouped:
        subdf = subdf.copy()

        # Identify which rows are QC
        qc_mask = (subdf[sampletype_col] == qc_label)
        qc_df = subdf[qc_mask].copy()

        # If we have fewer than 3 QC data points with meaningful intensities,
        # we can't reliably fit LOESS. We'll skip correction.
        # But let's also remove NaNs or 0 from the QC intensities before checking.
        qc_nonzero = qc_df.dropna(subset=[intensity_col])
        qc_nonzero = qc_nonzero[qc_nonzero[intensity_col] > 0]

        if len(qc_nonzero) < 3:
            # Not enough QC data to do LOESS. Just leave intensities as-is.
            subdf["Loess_Fit"] = np.nan
            subdf["Intensity_Corrected"] = subdf[intensity_col]
            corrected_subdfs.append(subdf)
            continue

        # Sort QCs by injection order, in case it's not sorted
        qc_nonzero.sort_values(by="InjectionOrder", inplace=True)

        # We'll do LOESS on injection order vs. intensity
        x_qc = qc_nonzero["InjectionOrder"].values
        y_qc = qc_nonzero[intensity_col].values

        # Fit LOESS
        # returns array of shape (n,2): [:,0] = x_sorted, [:,1] = y_smoothed
        loess_fit_array = sm.nonparametric.lowess(
            endog=y_qc,
            exog=x_qc,
            frac=frac,
            return_sorted=True
        )

        # Interpolate so we can predict for *all* injection orders
        f_interp = interp1d(
            loess_fit_array[:, 0],
            loess_fit_array[:, 1],
            kind='linear',
            fill_value="extrapolate"
        )

        # Predict for entire subdf (QC+Blanks+Samples)
        all_x = subdf["InjectionOrder"].values
        predicted = f_interp(all_x)
        subdf["Loess_Fit"] = predicted

        # We typically scale so the median QC intensity = median of the LOESS fit
        # i.e. if LOESS-Fit on QCs has median M, that becomes the "reference" level
        qc_loess_vals = f_interp(qc_df["InjectionOrder"].values)
        median_qc_loess = np.median(qc_loess_vals)

        # If median_qc_loess is extremely small (e.g. ~0), let's avoid dividing by zero
        if median_qc_loess <= 0:
            # It's possible that QC intensities are extremely small or basically zero
            # We'll skip scaling or just do raw ratio to predicted
            median_qc_loess = 1.0  # so we won't blow up in the next step

        # Correction formula:
        #  Intensity_Corrected = raw_intensity * (median_qc_loess / predicted_loess)
        # So that the corrected QC intensities revolve around 'median_qc_loess'
        subdf["Intensity_Corrected"] = subdf[intensity_col] * (median_qc_loess / subdf["Loess_Fit"])

        corrected_subdfs.append(subdf)

    # Combine results from all groups
    df_corrected = pd.concat(corrected_subdfs, axis=0, ignore_index=True)

    return df_corrected

def pqn_after_loess(df, 
                    feature_cols=["Molecule Name","Precursor Adduct"], 
                    intensity_col="Intensity_Corrected", 
                    sampletype_col="Sample Type", 
                    sample_label="Sample"):
    """
    Perform PQN (Probabilistic Quotient Normalization) on LOESS-corrected intensities
    for all 'Sample' rows, ignoring QCs/Blanks.

    Assumes 'df' already has a column with LOESS-corrected intensities.

    Parameters
    ----------
    df : pd.DataFrame
        Required columns:
          - feature_cols (list of column names) that identify each feature 
            (e.g. ["Molecule Name","Precursor Adduct"])
          - intensity_col (default "Intensity_Corrected")
          - sampletype_col (default "Sample Type"), which should indicate sample/QC/blank
        We'll only do PQN for rows where sampletype_col == sample_label.
    feature_cols : list
        Column names that together uniquely identify the feature (e.g. metabolite + adduct).
    intensity_col : str
        Column name of the corrected intensity to be PQN-normalized.
    sampletype_col : str
        Column that indicates whether row is QC, Blank, or Sample.
    sample_label : str
        The value in sampletype_col that means "study sample" (default "Sample").

    Returns
    -------
    pd.DataFrame
        Copy of input with an additional column "Intensity_PQN" for the final normalized intensities.
    """

    df = df.copy()

    # 1. Identify study sample rows
    sample_mask = (df[sampletype_col] == sample_label)
    sample_df = df[sample_mask].copy()
    
    # 2. Pivot the data so that each row is a sample, each column is a (feature) => intensity
    #    We'll combine feature_cols into a single key if multiple columns
    df["Feature_ID"] = df[feature_cols].astype(str).agg("|".join, axis=1)  # e.g. "MoleculeA|M+H"
    sample_df["Feature_ID"] = sample_df[feature_cols].astype(str).agg("|".join, axis=1)

    pivoted = sample_df.pivot_table(
        index="Replicate Name", 
        columns="Feature_ID", 
        values=intensity_col, 
        aggfunc='mean'  # or sum, if duplicates
    )

    # pivoted is now shape (n_samples x n_features)

    # 3. Compute reference spectrum (e.g. median across samples for each feature)
    reference_spectrum = pivoted.median(axis=0)  # yields a Series of length = # features

    # 4. Compute ratio matrix R_{s,f} = I_{s,f} / Reference_f
    ratio = pivoted.div(reference_spectrum, axis=1)

    # 5. For each sample, find the median ratio across features => scale factor
    #    We skip features with NaN if the sample had zero or missing intensities
    sample_scale_factors = ratio.median(axis=1)

    # 6. PQN scaling: I_{s,f}^{(PQN)} = I_{s,f}^{(corrected)} / scale_factor(s)
    pqn_matrix = pivoted.div(sample_scale_factors, axis=0)

    # 7. We'll melt back to long-form to merge with original df
    pqn_long = pqn_matrix.reset_index().melt(
        id_vars="Replicate Name",
        var_name="Feature_ID",
        value_name="Intensity_PQN"
    )

    # 8. Merge PQN results back into the original DataFrame
    #    For non-sample rows (QCs, Blanks), we won't have PQN values => fill with NaN
    df_merged = pd.merge(
        df,
        pqn_long,
        on=["Replicate Name", "Feature_ID"],
        how="left"
    )

    return df_merged



# Example usage:
# corrected_df = qc_rlsc_loess(my_df, frac=0.3)
# corrected_df.head()





#exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/proline_202412191102_deprecated'
exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_202501281618'
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
#data = pd.read_csv(f'/Users/danbru/Downloads/Transition Results newlist.csv')
#data = pd.read_csv('/Users/danbru/Downloads/Transition Results (1).csv')
#data = pd.read_csv('/Users/danbru/Downloads/glutamate_202501281618/autonoms_output/Positive/output_report.tsv', sep = '\t').drop(['Height',
#                                                                                                                                  'Precursor Charge',
#                                                                                                                                  'Collisional Cross Section'], axis = 1)
data = pd.read_csv('/Users/danbru/Downloads/output_report.tsv', sep = '\t').drop(['Height','Precursor Charge','Collisional Cross Section'], axis = 1)




#data = pd.read_csv(f'{exp_path}/results/metabolomics/TransitionResults.csv')
data['Well'] = data['Replicate Name'].str.extract(r'-([A-H]\d{1,2})')
#data['Peptide'] = data['Peptide'].str.replace('L-Glutamic acid','L-Glutamic Acid')

# Add the blanks as well
# Change values in 'Well' column to 'Blank' if 'Replicate' contains 'MAT1'
data.loc[data['Replicate Name'].str.contains('MAT1', na=False), 'Well'] = 'Blank'
data.loc[data['Replicate Name'].str.contains('MAT2', na=False), 'Well'] = 'QC'
data = data.dropna(subset = 'Well')






#outlier_wells = detect_outlier_wells(data, metric='TIC', threshold=2)
#outlier_wells.remove('QC')
#outlier_wells.remove('Blank')

#print("Outlier wells:", outlier_wells)

#data = data[~data['Well'].isin(outlier_wells)]

#data = data[data['Protein'].str.contains('AllCCS')]

data = data.merge(layout, left_on = 'Well', right_on = 'well', how = 'outer')

outliers_per_group, outliers = detect_outlier_wells_by_group(data, group_col='Summary', metric='TIC', threshold=3)
data = data[~data['Well'].isin(outliers)]
data = data[data['Well'] != 'A5']
data = data[data['Well'] != 'A9']


#data = data.drop_duplicates(subset = ['Molecule Name','Precursor Mz', 'Well'], keep = 'first')
#filtered_data = data[data['Area'] > 0.0]
#normalized_data = probabilistic_quotient_normalization_long_format(filtered_data).dropna().drop_duplicates()


data_corrected = qc_rlsc_loess_per_adduct(
    data,
    replicate_col="Replicate Name",
    molecule_col="Molecule Name",
    adduct_col="Precursor Adduct",
    intensity_col="Area",
    sampletype_col="Well",
    qc_label="QC",
    frac=0.3
)

data_corrected['Area'] = data_corrected['Intensity_Corrected']
normalized_data = quantile_normalization(
    data,
    well_col="Well",
    replicate_col="Replicate Name",
    molecule_col="Molecule Name",
    adduct_col="Precursor Adduct",
    intensity_col="Area"
)





pqn_after_loess(data_corrected, 
                    feature_cols=["Molecule Name","Precursor Adduct"], 
                    intensity_col="Area", 
                    sampletype_col="Well", 
                    sample_label="Sample")




summary = summarize_all_adducts_no_scoring(data_corrected)
best_adducts, filtered_data, adduct_counts = select_adduct_study_samples(data_corrected, restrict_with_blank=True, min_count=1)




summary = summarize_all_adducts_no_scoring(data)
best_adducts, filtered_data, adduct_counts = select_adduct_study_samples(data, restrict_with_blank=True, min_count=1)



#best_adducts, filtered_data, adduct_group_counts, all_adducts_summary = select_adduct_most_represented_across_groups(data, important_groups = None)


#best_adducts, filtered_data, adduct_group_counts, all_adducts_summary = select_adduct_most_represented_across_groups(data.dropna(), important_groups = ['Control group with no supplementation or treatment.',
#                                                                                                                                                        'Low proline supplementation without lactic acid treatment.',
#                                                                                                                                                        'High proline supplementation without lactic acid treatment.'])
#best_adducts, filtered_data, adduct_group_counts, all_adducts_summary = select_adduct_most_represented_across_groups(data, important_groups = ['Control group with no supplementation or treatment.',
#                                                                                                                                                        'Low proline supplementation without lactic acid treatment.',
#                                                                                                                                                        'High proline supplementation without lactic acid treatment.'])
#import pandas as pd
#import numpy as np


# Impute with LoD/2
#imputed_data = impute_lod_half_per_peptide_if_under_50pct(filtered_data, 0.3)
#imputed_data
#imputed_data = filtered_data
#test_data = filtered_data[filtered_data['Peptide'] == 'Proline']
# 1. Percentage of zero values per metabolite
#percent_zeros = filtered_data.groupby('Peptide')['Area'].apply(
#    lambda x: (x == 0).sum() / len(x) * 100
#).reset_index(name='PercentZeros')

# 2. Percentage of missing (NaN) values per metabolite
#percent_missing = filtered_data.groupby('Peptide')['Area'].apply(
#    lambda x: x.isna().sum() / len(x) * 100
#).reset_index(name='PercentMissing')#

# Merge these results into a single DataFrame
#summary_df = pd.merge(percent_zeros, percent_missing, on='Peptide')
#summary_df = 


#print(summary_df)


#best_adducts, filtered_data, adduct_group_counts, all_adducts_summary = select_adduct_most_represented_across_groups(data.dropna(), important_groups = ['Control group with no supplementation or treatment.'])
#filtered_test_data = filtered_data[filtered_data['Peptide'].isin(['Proline'])]




#normalized_data = qc_based_pqn_normalization(filtered_data, min_peptide_valid_fraction=0.5, min_valid_ratios=10)

imputed_data = impute_lod_half_per_peptide_if_under_50pct(filtered_data, 0.3)
#filtered_data = imputed_data[imputed_data['Area'] > 0.0]
#filtered_data['Area'] = filtered_data['Area']/filtered_data['Total Ion Current Area']
#normalized_data = probabilistic_quotient_normalization_long_format(filtered_data).dropna().drop_duplicates()
#filtered_data = filtered_data[filtered_data['Area'] > 0.0]
normalized_data = qc_based_pqn_normalization(imputed_data)


normalized_data = qc_based_pqn_normalization(filtered_data)
normalized_data = normalized_data[normalized_data['NormalizedArea'] > 0.0]


#normalized_data = qc_based_pqn_normalization(filtered_data)
normalized_data = normalized_data.merge(data[['Well','Summary']].drop_duplicates(), left_on = 'Well', right_on = 'Well') 
#normalized_data = normalized_data.merge(layout[['well','Summary']], left_on = 'Well', right_on = 'well')


#normalized_data = normalized_data[normalized_data['Peptide'] == 'Proline']

plot_data = normalized_data[['Well','Molecule Name','NormalizedArea','Summary']]


#plot_data = plot_data[plot_data['Peptide'].isin(['Glutamic acid','Proline','L-Arginine','L-Cysteine','L-Tryptophan',
#                                                 'L-Valine','Phenylalanine','Glutamine','Glutathione'])]
#plot_data = plot_data[plot_data['Peptide'].isin(['L-Glutamic Acid', 'L-Arginine', 'L-Proline'])]

plot_data = plot_data[plot_data['Molecule Name'].isin(['Glutamic Acid'])]
#plot_data = plot_data[plot_data['Molecule Name'].isin(['Glutamic Acid'])]




summary_order = []

# Group by Peptide and Summary and calculate the mean normalized intensity
# Group by Peptide and Summary and calculate the mean normalized intensity
grouped_data = plot_data.groupby(['Molecule Name', 'Summary'])['NormalizedArea'].mean().reset_index()

# Create a complete set of Peptide and Summary combinations to ensure alignment
all_peptides = grouped_data['Molecule Name'].unique()
all_summaries = summary_order if summary_order else grouped_data['Summary'].unique()
complete_index = pd.MultiIndex.from_product([all_peptides, all_summaries], names=['Molecule Name', 'Summary'])
grouped_data = grouped_data.set_index(['Molecule Name', 'Summary']).reindex(complete_index).reset_index()

# Apply specific order for the Summary column if provided
if summary_order:
    grouped_data['Summary'] = pd.Categorical(grouped_data['Summary'], categories=summary_order, ordered=True)
    grouped_data = grouped_data.sort_values(['Molecule Name', 'Summary'])

# Create a bar plot
plt.figure(figsize=(12, 6))
peptides = grouped_data['Molecule Name'].unique()
x = np.arange(len(peptides))  # the label locations
width = 0.8 / len(all_summaries)  # the width of the bars

# Separate data by Summary groups
for i, summary in enumerate(all_summaries):
    summary_data = grouped_data[grouped_data['Summary'] == summary]
    plt.bar(x + (i - len(all_summaries) / 2) * width, summary_data['NormalizedArea'], width=width, label=summary)

# Add labels, title, and legend
plt.xlabel('Molecule')
plt.ylabel('Mean Normalized Area')
plt.title('Normalized Areas Grouped by Peptide and Summary')
plt.xticks(x, peptides, rotation=45, ha='right')
plt.legend(title='Summary')
plt.yscale('log')
plt.tight_layout()
plt.show()


# Calculate p values
pval_df = calculate_p_values(plot_data)
pval_df = pval_df[pval_df.iloc[:, 1:3].apply(lambda x: x.astype(str).str.contains('Baseline', case=False, na=False)).any(axis=1)]


