#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 14:51:58 2025

@author: danbru
"""

import re
from pathlib import Path

import numpy as np
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler

from scipy.stats import f, chi2


def select_lowest_variation_qc(df, sample_labels):
    """
    For each molecule selects the feature that has the lowest RSD

    Parameters:
        df (pd.DataFrame): DataFrame containing the features.
        sample_labels (pd.Series): Series of sample labels corresponding to the rows of df.
        min_coverage (float): Minimum fraction of non-missing values in non-blank samples required.
        
    Returns:
        pd.DataFrame: A DataFrame containing only the selected features.
    """
    
    # Create masks for QC samples and non-blank samples
    qc_mask = sample_labels == "QC"
    
    selected_features = []
    molecule_dict = {}

    # Group features by molecule name (everything before " [")
    for feat in df.columns:
        molecule = feat.split(" [")[0]
        molecule_dict.setdefault(molecule, []).append(feat)
        
    selected_features = []
    for molecule, feats in molecule_dict.items():
        
        # compute mean and SD in QCs
        stats = {}
        for feat in feats:
            vals = df.loc[qc_mask, feat].dropna()
            
            if len(vals) < 0.5*len(df.loc[qc_mask, feat]): # Ensure some coverage in QCs
                continue
            
            mean = np.mean(vals)
            sd   = np.std(vals, ddof=1)
            rsd  = (sd / mean) * 100.0
            stats[feat] = rsd
    
        # pick lowest %RSD
        if stats:
            best_feat = min(stats, key=stats.get)
            selected_features.append(best_feat)

    return df[selected_features]


def perform_pqn_normalization(data, sample_labels):
    
    """
    Perform PQN normalization with DataFrame outputs, using QC samples as reference if available.

    Returns a normalized DataFrame of biological samples and a boolean mask.
    """

    # Ensure inputs are pandas objects
    df = data.copy()
    labels = pd.Series(sample_labels, index=df.index)

    # Masks
    qc_mask = labels.eq("QC")
    bio_mask = ~labels.isin(["QC", "Blank"])

    qc_df = df.loc[qc_mask]
    total_int = qc_df.sum(axis=1)
    med = total_int.median()
    mad = (total_int - med).abs().median() or 1e-9
    good_qc = qc_df.loc[((total_int - med).abs() / mad) < 3]

    # Identify bad features
    frac_bad = ((good_qc == 0) | good_qc.isna()).mean()
    bad_feats = frac_bad[frac_bad > 0.5].index.tolist()
    print(f"{len(bad_feats)} features have >50% zeros or NaNs in QC:")
    for feat in bad_feats:
        print("  –", feat)

    # Filter features
    keep = frac_bad.le(0.5)
    df = df.loc[:, keep]
    good_qc = good_qc.loc[:, keep]

    ref = good_qc.median()
     
    # Protect ref, should never happen though
    ref = ref.replace({0: np.nan}).fillna(1.0)

    # PQN normalization
    quotients = df.div(ref, axis=1)
    scaling = quotients.median(axis=1).replace({0: np.nan}).fillna(1.0)
    norm_df = df.div(scaling, axis=0)

    # Final biological subset and mask
    final_df = norm_df.loc[bio_mask]
    valid_mask = bio_mask

    return final_df, valid_mask

def perform_pqn_normalization(data, sample_labels):
    """
    Perform PQN normalization, skipping NaNs, filtering by QC, and normalizing only biological samples.

    Parameters
    ----------
    data : pd.DataFrame
        Samples × features intensity matrix.
    sample_labels : list-like or pd.Series
        Labels for each sample; should include "QC" and optionally "Blank".

    Returns
    -------
    final_df : pd.DataFrame
        PQN-normalized intensities for biological samples.
    valid_mask : pd.Series (bool)
        Mask indicating biological samples in the original data.
    """
    # 1) Setup
    df = data.copy()
    labels = pd.Series(sample_labels, index=df.index)

    qc_mask = labels.eq("QC")
    bio_mask = ~labels.isin(["QC", "Blank"])

    # 2) Select good QCs based on total intensity
    qc_df = df.loc[qc_mask]
    total_int = qc_df.sum(axis=1)
    med = total_int.median()
    mad = (total_int - med).abs().median() if (total_int - med).abs().median() > 0 else 1e-9
    good_qc = qc_df.loc[((total_int - med).abs() / mad) < 3]

    # 3) Drop features missing or zero in >50% of good QCs
    frac_bad = (good_qc.isna() | (good_qc == 0)).mean(axis=0)
    keep = frac_bad <= 0.5
    df = df.loc[:, keep]
    good_qc = good_qc.loc[:, keep]

    # 4) Compute reference spectrum (nan-aware)
    ref = good_qc.median(axis=0, skipna=True)
    # Drop features with zero median reference
    ref = ref.replace(0, np.nan).dropna()
    df = df.loc[:, ref.index]

    # 5) Compute quotients and scaling factors
    quotients = df.div(ref, axis=1)
    # Scaling only on biological samples, nan-aware median
    scaling = quotients.loc[bio_mask].apply(lambda row: np.nanmedian(row.values), axis=1)
    # Protect against zero or NaN scaling
    scaling = scaling.replace(0, np.nan).fillna(1.0)

    # 6) Normalize
    norm_df = df.div(scaling, axis=0)

    # 7) Return normalized biological samples and mask
    final_df = norm_df.loc[bio_mask]
    valid_mask = bio_mask

    return final_df, valid_mask



def apply_blank_threshold(df, sample_labels, N=3, sn_ratio = 1):

    blank_mask = sample_labels == "Blank"
    blank_mask[0] = ~blank_mask[0] # Just to remove the first blank of the series
    
    bio_mask = (sample_labels != "Blank") & (sample_labels != "QC")

    bio_idx = df.index[bio_mask]
    blank_idx = df.index[blank_mask]
        
    # Compute the median intensity of each feature across biological and blank samples
    mean_blank = df.loc[blank_idx].fillna(0).mean(skipna=True)
    mean_bio = df.loc[bio_idx].fillna(0).mean(skipna=True)
        
    # Identify features where the biological mean is lower than the blank mean
    remove_features = mean_bio < mean_blank*sn_ratio
    print(f'{len(df.columns[remove_features])} features removed due to low S/N ratio.')
    df_filtered = df.drop(columns=df.columns[remove_features])
        
    return df_filtered


def detect_pca_outliers(data, labels, var_explained=0.95, conf_level=0.95, res_level=0.95):
    
    # 1) Impute & scale…
    data = do_rf_imputation(data)
    X = StandardScaler().fit_transform(data)

    # 2) PCA
    pca = PCA(n_components=var_explained)
    scores = pca.fit_transform(X)
    lambdas = pca.explained_variance_
    n_pc = pca.n_components_
    N = X.shape[0]

    # 3) Hotelling T²
    T2 = np.sum((scores**2) / lambdas, axis=1)
    F_val = f.ppf(conf_level, dfn=n_pc, dfd=N-n_pc)
    T2_cutoff = (n_pc*(N-1)/(N-n_pc))*F_val

    # 4) get something similar to DModX via chi²
    SPE = np.sum((X - pca.inverse_transform(scores))**2, axis=1)
    df_spe = X.shape[1] - n_pc  
    theta = SPE.mean() / df_spe  # rough avg eigenvalue
    SPE_cutoff = theta * chi2.ppf(res_level, df_spe)
    DModX = np.sqrt(SPE)
    DModX_cutoff = np.sqrt(SPE_cutoff)

    # 5) Flag & return…
    is_outlier = (T2 > T2_cutoff) | (DModX > DModX_cutoff)
    mask_inlier = ~is_outlier
    outlier_info = pd.DataFrame({
        'T2': T2, 'T2_cutoff': T2_cutoff,
        'DModX': DModX, 'DModX_cutoff': DModX_cutoff,
        'Outlier': is_outlier
    }, index=data.index)
    
    return (data.loc[mask_inlier],
            labels.loc[mask_inlier],
            mask_inlier,
            outlier_info)

# Update only rows with MAT1 or MAT2 in the material tag
def assign_exp_group(row):
    if row['well'] == "MAT1":
        return "Blank"
    elif row['well'] == "MAT2":
        return "QC"
    else:
        # Leave the existing value unchanged.
        return row["Experimental group"]

def drop_high_nan_columns(df, threshold=1/3):
    nan_ratio = df.replace(0.0, np.nan).isna().mean() # This to ensure that ultra low intensities also count as missing in thsi case
    filtered_df = df.loc[:, nan_ratio <= threshold]
    return filtered_df

def get_wells_below_threshold(filepath, threshold=0.600):
    
    wells_below = []
    current_well = None  

    # Mapping from row number to well letter
    row_map = {i: chr(64 + i) for i in range(1, 9)}

    # Regex patterns:
    well_pattern = re.compile(r'<RapidFire>:\s+BatchThread:\s+process_well\(\)\s+\(\s*(\d+),\s*(\d+)\s*\)')
    timing_pattern = re.compile(r',\s*t\s*=\s*([\d\.]+)')

    with open(filepath, 'r') as f:
        for line in f:
            
            # Check if this line starts a new well processing block.
            well_match = well_pattern.search(line)
            if well_match:
                
                # Before starting a new block, reset any pending current_well.
                current_well = (int(well_match.group(1)), int(well_match.group(2)))
                continue

            # Check if this line marks the end of the block with a timing value.
            timing_match = timing_pattern.search(line)
            if timing_match and current_well is not None:
                try:
                    t_value = float(timing_match.group(1))
                except ValueError:
                    t_value = None

                # Only record the well if:
                # 1. The well's row is not 1000 (i.e. a QC or blank injection).
                # 2. The timing value exists and is below the threshold.
                
                row_num, col_num = current_well
                if row_num != 1000 and t_value is not None and t_value <= threshold:
                    well_letter = row_map.get(row_num, f'Row{row_num}')
                    well_name = f"{well_letter}{col_num}"
                    wells_below.append(well_name)

                # Reset current_well for the next block.
                current_well = None

    return wells_below

def do_rf_imputation(df):

    # Do imputation
    imputed = df.copy()
   
    # Select numeric columns for imputation
    numeric_cols = imputed.select_dtypes(include=[np.number]).columns
    
    imputer = IterativeImputer(
        estimator=RandomForestRegressor(n_estimators=100, random_state=42),
        max_iter=10,
        initial_strategy="mean",
        imputation_order="ascending",
        skip_complete=True,
        random_state=0
    )
    
    # Perform imputation on numeric columns
    imputed[numeric_cols] = imputer.fit_transform(imputed[numeric_cols])
    
    return imputed

def merge_columns_by_value_agreement(df: pd.DataFrame, match_thresh: float = 0.9) -> pd.DataFrame:

    cols = list(df.columns)
    n = len(df)
    visited = set()
    groups = []

    # Identify groups
    for i, c1 in enumerate(cols):
        if c1 in visited:
            continue
        group = [c1]
        visited.add(c1)
        for c2 in cols[i+1:]:
            if c2 in visited:
                continue
            # fraction of rows where both non-null and equal, relative to total rows
            mask = df[c1].notna() & df[c2].notna()
            frac_match = (df.loc[mask, c1] == df.loc[mask, c2]).sum() / n
            if frac_match >= match_thresh:
                group.append(c2)
                visited.add(c2)
        groups.append(group)

    # Print merged groups
    for group in groups:
        if len(group) > 1:
            print(f"Merged columns: {', '.join(group)} into { '/'.join(group) }")

    # Build merged DataFrame
    merged = pd.DataFrame(index=df.index)
    for group in groups:
        if len(group) == 1:
            merged[group[0]] = df[group[0]]
        else:
            name = "/".join(group)
            merged[name] = df[group].mean(axis=1)
    return merged


##############################################
# 3. GLOBAL PREPROCESSING PIPELINE
##############################################

def global_preprocess(merged_df, exp_path, params):
    """
    Identical to your original pipeline. The only addition:
    - We explicitly preserve 'ReplicateName' as a separate column
      (so it does NOT get treated as numeric).
    """

    # Drop columns that are entirely NaN
    merged_df = merged_df.dropna(axis=1, how='all')
    
    # 1) Subset to QC rows
    qc_df = merged_df[merged_df['Experimental group'] == 'QC']
    frac_bad = qc_df[qc_df.replace(0, np.nan).select_dtypes(include=[np.number]).columns].isna().mean()
    cols_to_drop = frac_bad[frac_bad > 0.4].index.tolist()
    
    print(f'Following columns dropped due to low presence in QCs: \n {cols_to_drop}')
    
    merged_df = merged_df.drop(columns=cols_to_drop)


    # Extract numeric features and sample labels
    drop_cols = ["Experimental group", "ReplicateName",'Total Ion Current Area']
    features = merged_df.drop(columns=drop_cols, errors="ignore")
    numeric_features = features.select_dtypes(include=[np.number])
    sample_labels = merged_df["Experimental group"]

   
    # 1) Blank-based filtering
    blank_filtered = apply_blank_threshold(numeric_features, 
                                           sample_labels,
                                           params['n_blanks'], 
                                           params['sn_ratio'])
    
    # 2) Remove blank samples from further analysis
    non_blank_mask = sample_labels != "Blank"
    blank_filtered = blank_filtered.loc[non_blank_mask]
    sample_labels = sample_labels.loc[non_blank_mask]
    replicate_series = merged_df.loc[non_blank_mask, "ReplicateName"]

    # 3) Drop columns with too many NaNs
    blank_filtered = drop_high_nan_columns(blank_filtered, 
                                           threshold=params['missing_threshold'])
    
    # 3.5) # impute with LoD/2 for peaks that are detected, but set at 0. Important to do this AFTER missingness removal!
    blank_filtered = blank_filtered.apply(lambda col: col.replace(0, col[col != 0].min() / 2) if (col != 0).any() else col)

    # 4) Select the best adducts
    reduced_features = select_lowest_variation_qc(blank_filtered, 
                                                  sample_labels)
    
    # 4.5) # Merge the metabolites which we cant reliably separate with our resolution
    reduced_features = merge_columns_by_value_agreement(reduced_features, match_thresh = 2/3)

    # 4.75) Remove samples with more than % missing, to avoid super low concentration samples
    initial_rows = reduced_features.shape[0]
    reduced_features = reduced_features[reduced_features.isnull().mean(axis=1) <= params['missing_per_sample']]
    print(f'{initial_rows - reduced_features.shape[0]} injections removed due to feature missingness.')

    # 5) PQN normalization with QC quality check
    pqn_arr, valid_mask = perform_pqn_normalization(reduced_features, sample_labels)
    
    # Create a new DataFrame for normalized data using the filtered index
    filtered_index = reduced_features.index[valid_mask]
    pqn_df = pd.DataFrame(pqn_arr, columns=reduced_features.columns, index=filtered_index)
    
    cleaned_features, cleaned_labels, valid_mask, outlier_info = detect_pca_outliers(pqn_df, sample_labels.loc[filtered_index],
                            var_explained=params['pca_num_pcs'],
                            conf_level=params['pca_threshold'],
                            res_level=params['res_threshold'])
    
    outlier_info.to_csv(Path(exp_path) / 'results/metabolomics/processed/outlier_report.tsv', sep = '\t')
    
    
    print(f'{pqn_df.shape[0]-cleaned_features.shape[0]} outliers removed.')

    # Update sample_labels and replicate_series accordingly
    filtered_sample_labels = sample_labels.loc[filtered_index].loc[valid_mask]
    replicate_series = merged_df.loc[filtered_index, "ReplicateName"].loc[valid_mask]

    # Return final numeric feature DataFrame, final labels, plus replicate names
    return cleaned_features, filtered_sample_labels, replicate_series




def main_method(exp_folder):
    
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/partially_completed/glutamine_202504251440' # glutamine_acetate
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/partially_completed/lysine_202504291748' # lysine sucrose
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/partially_completed/aminoadipate_202504291411' # aminoadipate
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756' # FA
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539' # Caffeine
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618' # Spermine
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407' # Lactic acid
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655' # LiCl
    exp_paths = [
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamine_202504251440', # glutamine_acetate
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/lysine_202504291748', # lysine sucrose
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/aminoadipate_202504291411' # aminoadipate
    ]
    for exp_folder in exp_paths:
    
        EXPERIMENT_DIR = Path(exp_folder)
        
        # Remember to move these to global config later
        params = {
            'min_qc': 5, # Just to make sure we have enough QCs to actually do PQN
            'n_blanks': 6, # N closest blanks in the runorder, used to do the blank-filtering
            'pca_num_pcs': 0.95, # How many PCA dimensions to use for outlier removal
            'pca_threshold': 0.95, # T2 threshold for outlier removal
            'res_threshold': 0.99, # threshold for unexplained residual removal
            'missing_threshold': 1/3, # inverse of the one above, dummy
            'sip_threshold': 0.600, # Threshold for outlier removal based on rapidfire sensor values
            'sn_ratio': 1,
            'missing_per_sample': 0.5 # missingness threshold per sample. i.e. if one sample has more than X% missing peaks, it is likely not a reliable one
        }
    
        # Paths
        data_path = EXPERIMENT_DIR / 'results/metabolomics/raw/output_report.tsv'
        metadata_path = EXPERIMENT_DIR / 'protocol/plate_layout/layout.tsv'
    
        # Load raw data
        raw_df = pd.read_csv(data_path, sep='\t')
        metadata = pd.read_csv(metadata_path, sep = '\t', index_col = 0).drop(columns=['plateID','CONCuM'])
        metadata.columns = ['well','Experimental group']
            
        # Create single "Feature" that describes the adduct
        raw_df["Feature"] = raw_df["Molecule Name"] + " [" + raw_df["Precursor Adduct"] + "]"
    
        # Pivot so that each Replicate Name becomes a row
        wide_df = raw_df.pivot(index="Replicate Name", columns="Feature", values="Area")
        wide_df["ReplicateName"] = wide_df.index
        wide_df = wide_df.reset_index(drop=True)  # So index is just 0..N
        
        # Formatting shit so we can merge with metadata
        metadata['well'] = metadata['well'].str.replace(r'([A-H])0+(\d+)', r'\1\2', regex=True)
        temp_df = raw_df[['Replicate Name','Total Ion Current Area']].copy()
        temp_df['well'] = temp_df['Replicate Name'].apply(lambda x: x.split('-')[2])
        temp_df['well'] = temp_df['well'].apply(lambda x: x.split('_')[0])
        temp_df['well'] = temp_df['well'].apply(lambda x: x.split('.')[0])
        temp_df = temp_df.merge(metadata, left_on = 'well', right_on = 'well', how ='left').drop_duplicates()
        temp_df["Experimental group"] = temp_df.apply(assign_exp_group, axis=1)
        metadata = temp_df
        
        # Merge raw_df and metadata on the "well" column 
        merged_df = pd.merge(wide_df, metadata, left_on='ReplicateName', right_on = 'Replicate Name').drop('Replicate Name', axis = 1)
        merged_df["Experimental group"] = merged_df.apply(assign_exp_group, axis=1)
        
        # Remove non-injections (as defined by the sipping-time parameter in RapidFire)
        valid_injections = get_wells_below_threshold(EXPERIMENT_DIR / 'results/metabolomics/raw/batch.log', threshold=params['sip_threshold'])
        merged_df = merged_df[merged_df['well'].isin(valid_injections+['MAT1']+['MAT2'])]    
        merged_df.drop(columns=["well"], inplace=True)
    
        # Remove any rows with no "Experimental group" 
        merged_df = merged_df.dropna(subset=['Experimental group'])
        
        # Do the global preprocessing step for outliers and normalization
        X_global, sample_labels_global, replicate_names_global = global_preprocess(merged_df, exp_folder, params)
        print(f'{X_global.shape[1]} metabolites left after curation.')
        
        # Build final processed DataFrame
        final_processed_df = X_global.copy()
        final_processed_df["Experimental group"] = sample_labels_global
    
        # Attach the replicate name to each row
        final_processed_df["ReplicateName"] = replicate_names_global
    
        # Save the final processed_df (with nans)
        final_processed_df.set_index('ReplicateName').to_csv(EXPERIMENT_DIR / 'results/metabolomics/processed/ms_output.tsv', sep = '\t')
        

        imputed_df = do_rf_imputation(final_processed_df)
        imputed_df.set_index('ReplicateName').to_csv(EXPERIMENT_DIR / 'results/metabolomics/processed/ms_output_imputed.tsv', sep='\t')
        
        print(imputed_df.shape)

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="Metabolomics processing")
    parser.add_argument("--output_folder", required=True, type=str, help="Experiment folder")

    args = parser.parse_args()
    main_method(args.output_folder)

