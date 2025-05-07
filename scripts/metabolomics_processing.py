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


def perform_pqn_normalization(
    data,
    sample_labels,
    use_qc_reference: bool = True,
    qc_cv_threshold: float = 0.6,
    min_qc_fraction: int = 0.3
):
    """
    Perform PQN normalization, using QC samples as a reference if available
    and of sufficient quality. QC outliers are first removed by a TIC filter.

    Parameters
    ----------
    data : array-like, shape (n_samples, n_features)
        Raw intensity data.
    sample_labels : array-like of length n_samples
        Labels indicating which rows are "QC", "Blank", or biological samples.
    use_qc_reference : bool, default=True
        Whether to attempt to use QC samples to build the reference spectrum.
    qc_cv_threshold : float, default=0.6
        Maximum acceptable median CV across features in the QC block.
    min_qc_samples : int, default=10
        Minimum number of “good” QC runs required to use QC reference.

    Returns
    -------
    norm_data_filtered : ndarray, shape (n_bio_samples, n_features)
        PQN‑normalized data, excluding QC rows.
    valid_mask : ndarray, shape (n_samples,), dtype=bool
        Boolean mask indicating which original rows (by index) are retained
        (i.e., non-QC samples).
    """
    
    # Convert to numeric array
    data_arr = np.array(data, dtype=np.float64)
    
   

    # Mask to identify biological samples (non-QC, non-Blank)
    bio_mask = ~sample_labels.isin(["QC", "Blank"])

    # Decide on reference spectrum
    if use_qc_reference and ("QC" in sample_labels.values):
        qc_mask = (sample_labels == "QC")
        qc_data = data_arr[qc_mask]

        # 1) Simple (relevant) TIC-based outlier removal on QCs
        total_int = np.nansum(qc_data, axis=1)
        med = np.nanmedian(total_int)
        mad = np.nanmedian(np.abs(total_int - med)) or 1e-9
        keep_qc = np.abs((total_int - med) / mad) < 3
        qc_data = qc_data[keep_qc]
        
        # If you still have the DataFrame version of QC data, do:
        qc_df = data.loc[sample_labels == "QC"].iloc[keep_qc]   # shape (n_good_qc, n_features)
        zero_or_nan_frac = ((qc_df == 0) | qc_df.isna()).mean(axis=0)
        bad_feats = zero_or_nan_frac[zero_or_nan_frac > 0.5].index.tolist()
        
        print(f"{len(bad_feats)} features have >50% zeros or NaNs in QC:")
        for feat in bad_feats:
            print("  –", feat)
        
        # 2) Check minimum number of QCs after filtering
        if qc_data.shape[0] >= min_qc_fraction*len(sample_labels):
            # Compute per-feature CV in the cleaned QC block
            mean_qc = np.nanmean(qc_data, axis=0)
            std_qc  = np.nanstd(qc_data, axis=0)
            cv_qc   = std_qc / (mean_qc + 1e-9)
            median_cv = np.nanmedian(cv_qc)

            print(f"QC median CV: {median_cv:.3f}")
            if median_cv < qc_cv_threshold:
                print("QC quality-check passed. Using cleaned QCs as reference.")
                ref_spec = np.nanmedian(qc_data, axis=0)
            else:
                print("QC CV too high. Falling back to biological samples.")
                ref_spec = np.nanmedian(data_arr[bio_mask], axis=0)
        else:
            print("Too few QC-samples. Falling back to biological samples.")
            ref_spec = np.nanmedian(data_arr[bio_mask], axis=0)
    else:
        print("Using sample spectra as reference.")
        ref_spec = np.nanmedian(data_arr[bio_mask], axis=0)

    # Protect against zeros or NaNs in the reference
    ref_spec = np.where(np.logical_or(np.isnan(ref_spec), ref_spec == 0), 1.0, ref_spec)

    # 3) PQN scaling
    #    Compute sample-wise quotients relative to ref_spec
    quotients = data_arr / ref_spec[np.newaxis, :]
    quotients[np.isnan(data_arr)] = np.nan  # keep NaNs where data was missing

    #    Per-sample scaling factor is the median of those quotients
    scaling_factors = np.nanmedian(quotients, axis=1)
    scaling_factors = np.where(np.isnan(scaling_factors), 1.0, scaling_factors)

    #    Normalize data
    norm_data = data_arr / scaling_factors[:, np.newaxis]
    norm_data[np.isnan(data_arr)] = np.nan

    # 4) Filter out QC rows from the final output
    valid_mask = ~sample_labels.isin(["QC", "Blank"])
    norm_data_filtered = norm_data[valid_mask]

    return norm_data_filtered, valid_mask



def apply_blank_threshold(df, sample_labels, N=3, sn_ratio = 1, remove_by_feature=False):

    blank_mask = sample_labels == "Blank"
    blank_mask[0] = ~blank_mask[0] # Just to remove the first blank of the series
    
    bio_mask = (sample_labels != "Blank") & (sample_labels != "QC")

    if remove_by_feature:
        
        bio_idx = df.index[bio_mask]
        blank_idx = df.index[blank_mask]
        
        # Compute the median intensity of each feature across biological and blank samples
        mean_blank = df.loc[blank_idx].fillna(0).mean(skipna=True)
        mean_bio = df.loc[bio_idx].fillna(0).mean(skipna=True)
        
        # Identify features where the biological mean is lower than the blank mean
        remove_features = mean_bio < mean_blank*sn_ratio
        print(f'{len(df.columns[remove_features])} features removed due to low S/N ratio.')
        df_filtered = df.drop(columns=df.columns[remove_features])
        
    else: 
        
        # Peak-level (per-sample) filtering.
        blank_idx = df[sample_labels == "Blank"].index
        df_filtered = df.copy()
        for i in df.index:
            if sample_labels.loc[i] == "Blank":
                continue
            distances = np.abs(blank_idx - i)
            closest = blank_idx[np.argsort(distances)][:N]
            
            # Calculate the median for each feature from the N closest blanks.
            median_blank = df.loc[closest].fillna(0).median(skipna=True)
            
            # Replace values in the sample that are below this median with NaN.
            df_filtered.loc[i] = df.loc[i].where(df.loc[i] >= median_blank*sn_ratio, np.nan)
    return df_filtered


def detect_pca_outliers(data, labels,
                        n_components=5,
                        conf_level=0.95,
                        residual_threshold=0.95):
    
    from scipy.stats import f

    """
    Identify outliers via two criteria:
      1. Hotelling's T² in PCA score space.
      2. High reconstruction error (residual variance).

    We enforce p <= min(n_features, n_samples-1) to keep the F-statistic valid.
    """
    # 1) Impute missing values (you could swap in a robust or EM approach here)
    #pre_imp_data
    imputed = do_rf_imputation(data)

    N, F = imputed.shape

    # 2) Determine how many components we can actually fit:
    pca = PCA(n_components=n_components)
    scores = pca.fit_transform(imputed)

    # 3) Hotelling T²: sum of (score_i)^2 / eigenvalue_i
    lambdas = pca.explained_variance_
    T2 = np.sum((scores ** 2) / lambdas, axis=1)

    # 4) T² cutoff (F-distribution)
    p = pca.n_components_
    F_val = f.ppf(conf_level, dfn=p, dfd=N - p)
    T2_cutoff = (p * (N - 1) / (N - p)) * F_val

    # 5) Reconstruction residuals
    reconstructed = pca.inverse_transform(scores)
    residuals = np.sum((imputed - reconstructed) ** 2, axis=1)
    residual_cutoff = np.nanquantile(residuals, residual_threshold)

    # 6) Combine criteria
    outlier_mask = (T2 > T2_cutoff) | (residuals > residual_cutoff)

    # 7) Report
    outlier_info = pd.DataFrame({
        'Hotelling_T2': T2,
        'T2_cutoff': T2_cutoff,
        'Residuals': residuals,
        'Residual_cutoff': residual_cutoff,
        'Outlier': outlier_mask
    }, index=data.index)

    cleaned_data = data.loc[~outlier_mask]
    cleaned_labels = labels.loc[~outlier_mask]

    return cleaned_data, cleaned_labels, ~outlier_mask, outlier_info

def drop_high_nan_columns(df, threshold=1/3):
    """
    Drop columns with > threshold fraction of NaNs.
    """
    nan_ratio = df.replace(0.0, np.nan).isna().mean() # This to ensure that ultra low intensities also count as missing in thsi case
    filtered_df = df.loc[:, nan_ratio <= threshold]
    return filtered_df

def get_wells_below_threshold(filepath, threshold=0.600):
    wells_below = []
    current_well = None  # Will store a tuple (row, col) from the process_well line

    # Mapping from row number to well letter (1->A, 2->B, ..., 8->H)
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
                # 1. The well's row is not 1000.
                # 2. The timing value exists and is below the threshold.
                row_num, col_num = current_well
                if row_num != 1000 and t_value is not None and t_value < threshold:
                    well_letter = row_map.get(row_num, f'Row{row_num}')
                    well_name = f"{well_letter}{col_num}"
                    wells_below.append(well_name)

                # Reset current_well for the next block.
                current_well = None

    return wells_below

def do_rf_imputation(df):

    # 8) Do imputation
    imputed = df.copy()
   
    # Select numeric columns for imputation
    numeric_cols = imputed.select_dtypes(include=[np.number]).columns
    
    imputer = IterativeImputer(
        estimator=RandomForestRegressor(n_estimators=100, random_state=0),
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
    """
    Merge columns whose values agree in at least `match_thresh` fraction of all rows.
    Prints the groups of columns that have been merged.

    Parameters:
    - df: DataFrame with samples as rows and features as columns.
    - match_thresh: fraction of total rows where two columns must have identical values to merge.

    Returns:
    - DataFrame with merged columns (mean of group) named by joining originals with '/'.
    """
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

def global_preprocess(merged_df, params):
    """
    Identical to your original pipeline. The only addition:
    - We explicitly preserve 'ReplicateName' as a separate column
      (so it does NOT get treated as numeric).
    """

    # Drop columns that are entirely NaN
    merged_df = merged_df.dropna(axis=1, how='all')

    # Extract numeric features and sample labels
    drop_cols = ["Experimental group", "ReplicateName",'Total Ion Current Area']
    features = merged_df.drop(columns=drop_cols, errors="ignore")
    numeric_features = features.select_dtypes(include=[np.number])
    sample_labels = merged_df["Experimental group"]

    # 1) Select one adduct per molecule using highest SNR
    #numeric_features = select_highest_snr_adduct(numeric_features, sample_labels,
    #                                             min_coverage=params['min_coverage'])
    
    
    # 1) Blank-based filtering
    blank_filtered = apply_blank_threshold(numeric_features, 
                                           sample_labels,
                                           params['n_blanks'], 
                                           params['sn_ratio'],
                                           remove_by_feature=params['remove_by_feature'])
    
    # 2) Remove blank samples from further analysis
    non_blank_mask = sample_labels != "Blank"
    blank_filtered = blank_filtered.loc[non_blank_mask]
    sample_labels = sample_labels.loc[non_blank_mask]
    replicate_series = merged_df.loc[non_blank_mask, "ReplicateName"]

    # 3) Drop columns with too many NaNs
    blank_filtered = drop_high_nan_columns(blank_filtered, 
                                           threshold=params['missing_threshold'])
    
    # 3.5) # impute with LoD/2 for peaks that are detected, but set at 0
    # Important to do this after missingness removal
    blank_filtered = blank_filtered.apply(lambda col: col.replace(0, col[col != 0].min() / 2) if (col != 0).any() else col)

    
    # 4) Select the best adducts
    reduced_features = select_lowest_variation_qc(blank_filtered, 
                                                  sample_labels)
    
    # 4.5) # Merge the metabolites which we cant reliably separate?
    reduced_features = merge_columns_by_value_agreement(reduced_features, match_thresh = 0.5)

    # 4.75) Remove samples with more than % missing
    initial_rows = reduced_features.shape[0]
    reduced_features = reduced_features[reduced_features.isnull().mean(axis=1) <= params['missing_per_sample']]
    print(f'{initial_rows - reduced_features.shape[0]} injections removed due to feature missingness.')

    # 5) PQN normalization with QC quality check
    pqn_arr, valid_mask = perform_pqn_normalization(reduced_features, sample_labels,
                                        use_qc_reference=params['use_qc_reference'],
                                        qc_cv_threshold=params['qc_cv_threshold'],
                                        min_qc_fraction=params['min_qc_fraction'])
    
    # Create a new DataFrame for normalized data using the filtered index
    filtered_index = reduced_features.index[valid_mask]
    pqn_df = pd.DataFrame(pqn_arr, columns=reduced_features.columns, index=filtered_index)
    
    # Remove outliers
    cleaned_features, cleaned_labels, valid_mask, outlier_info = detect_pca_outliers(
        pqn_df, 
        sample_labels.loc[filtered_index],
        n_components=params['pca_num_pcs'],
        conf_level = params['pca_threshold'],
        residual_threshold = params['res_threshold'],
    )
    
    print(f'{pqn_df.shape[0]-cleaned_features.shape[0]} outliers removed.')

    # Update sample_labels and replicate_series accordingly
    filtered_sample_labels = sample_labels.loc[filtered_index].loc[valid_mask]
    replicate_series = merged_df.loc[filtered_index, "ReplicateName"].loc[valid_mask]

    # Return final numeric feature DataFrame, final labels, plus replicate names
    return cleaned_features, filtered_sample_labels, replicate_series

# Update only rows with MAT1 or MAT2 in the material tag
def assign_exp_group(row):
    if row['well'] == "MAT1":
        return "Blank"
    elif row['well'] == "MAT2":
        return "QC"
    else:
        # Leave the existing value unchanged.
        return row["Experimental group"]


def main_method(exp_folder):
    
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
    ]
    for exp_folder in exp_paths:
    
        EXPERIMENT_DIR = Path(exp_folder)
        
        # Remember to move these to global config later
        params = {
            'min_qc': 5, # Just to make sure we have enough QCs to actually do PQN
            'n_blanks': 6, # N closest blanks in the runorder, used to do the blank-filtering
            'use_qc_reference': True, # For PQN. If true, uses the QC spectras as ther reference (if qc passes quality checks), else it uses median spectra
            'pca_num_pcs': 5, # How many PCA dimensions to use for outlier removal
            'pca_threshold': 0.95, # T2 threshold for outlier removal
            'res_threshold': 0.99, # threshold for unexplained residual removal
            'missing_threshold': 1/3, # inverse of the one above, dummy
            'qc_cv_threshold': 0.7, # maximum acceptable median CV among QCs for PQN
            'min_qc_fraction': 0.2, # minimum number of QC samples required for PQN
            'remove_by_feature': True, # If we blanket-remove a feature if it is lower than in blank
            'sip_threshold': 0.600, # Threshold for outlier removal based on rapidfire sensor values
            'sn_ratio': 1,
            'missing_per_sample': 0.5 # missingness threshold per sample. i.e. if one sample has more than X% missing peaks, it is likely not a reliable one
        }
    
        # Paths
        #data_path = EXPERIMENT_DIR / 'results/metabolomics/raw/output_report_3000030lower5.tsv'
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
        X_global, sample_labels_global, replicate_names_global = global_preprocess(merged_df, params)
        #print(f'{X_global.shape[0]} metabolites left after curation.')
        print(f'{X_global.shape[1]} metabolites left after curation.')
        #missing_values = X_global.isna().mean() * 100
        
        # Build final processed DataFrame
        final_processed_df = X_global.copy()
        final_processed_df["Experimental group"] = sample_labels_global
    
        # Attach the replicate name to each row
        final_processed_df["ReplicateName"] = replicate_names_global
    
        # Save the final processed_df (with nans)
        #final_processed_df.set_index('ReplicateName').to_csv(EXPERIMENT_DIR / 'results/metabolomics/processed/ms_output.tsv', sep = '\t')
        final_processed_df.set_index('ReplicateName').to_csv(EXPERIMENT_DIR / 'results/metabolomics/processed/ms_output.tsv', sep = '\t')
        
        #imputed_df = do_rf_imputation(final_processed_df)
        #imputed_df = do_knn_imputation(final_processed_df)
        imputed_df = do_rf_imputation(final_processed_df)
        imputed_df.set_index('ReplicateName').to_csv(EXPERIMENT_DIR / 'results/metabolomics/processed/ms_output_imputed.tsv', sep='\t')
    

if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="Metabolomics processing")
    parser.add_argument("--output_folder", required=True, type=str, help="Experiment folder")

    args = parser.parse_args()
    main_method(args.output_folder)

