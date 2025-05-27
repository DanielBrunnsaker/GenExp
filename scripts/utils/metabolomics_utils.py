#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue May 27 08:58:01 2025

@author: danbru
"""

import re
import numpy as np
from pathlib import Path
import pandas as pd

from sklearn.decomposition import PCA
from sklearn.experimental import enable_iterative_imputer  # noqa
from sklearn.impute import IterativeImputer
from sklearn.ensemble import RandomForestRegressor
from sklearn.preprocessing import StandardScaler
from scipy.stats import f, chi2

from sklearn.linear_model import ElasticNetCV
from sklearn.model_selection import KFold, RepeatedKFold, RepeatedStratifiedKFold
from sklearn.metrics import r2_score
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
from rpy2 import robjects
from rpy2.robjects import StrVector
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, default_converter


def get_wells_below_threshold(filepath, threshold=0.600):
    """
    Reads the Rapidfire logs and removes outliers based on sensor-readings.

    Parameters
    ----------
    filepath : path/string
        Path to the batch-file.
    threshold : float, optional
        The sipping parameter. Typically set to the same
        value as was used for aspiration timing. The default is 0.600.

    Returns
    -------
    wells_below : list
        Wells that passed the sipping-parameter threshold.

    """
    
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

def apply_blank_threshold(df, sample_labels, sn_ratio = 1):
    """
    Removes features (peaks) that are below a set signal-to-noise ratio, 
    defined by sn_ratio.

    Parameters
    ----------
    df : dataframe
        Dataframe containing metabolomics-data in wide format.
    sample_labels : series
        Series containing sample descriptions.
    sn_ratio : float
        the signal-to-noise-ratio. The default is 1.

    Returns
    -------
    df_filtered : dataframe
        Dataframe containing metabolomics data in wide format, with peaks
        non-distinguishable from background removed.

    """

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

def drop_high_nan_columns(df, threshold=1/3):
    """
    Removes columns/features with too many missing values.

    Parameters
    ----------
    df : dataframe
        Metabolomics-data in wide format.
    threshold : float, optional
        What fraction of missing values to allow for 
        the peak to be considered valid. The default is 1/3.

    Returns
    -------
    filtered_df : dataframe
        Dataframe with invalid peaks removed.

    """
    nan_ratio = df.replace(0.0, np.nan).isna().mean() # This to ensure that ultra low intensities also count as missing in thsi case
    filtered_df = df.loc[:, nan_ratio <= threshold]
    return filtered_df


def select_lowest_variation_qc(df, sample_labels):
    """
    
    In case of several adducts per metabolite, selects the one with the lowest
    RSD in included QC-samples and returns them.
    
    Parameters
    ----------
    df : dataframe
        Dataframe containing metabolomics data, in wide format.
    sample_labels : series
        Series containing the sample types, denoting experimental group and/or 
        type of sample.

    Returns
    -------
    df[selected_features]: dataframe
        Returns a dataframe containing the selected adducts.

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

def merge_columns_by_value_agreement(df, match_thresh = 2/3):
    """
    Merges peaks that we cant reliably separate with current resolving power.
    Merges based on column similarity (not correlation). If 2/3 of rows
    are exactly the same, we assume that we cant reliably separate them.

    Parameters
    ----------
    df : pd.DataFrame
        Metabolomics data in wide format.
    match_thresh : float, optional
        The matching threshold, if X rows are exactly matching. The default is 2/3.

    Returns
    -------
    merged : dataframe
        Returns a dataframe with non-distinguishable columns merged. 
        Columns have their names merged as well.

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
    #for group in groups:
    #    if len(group) > 1:
    #        print(f"Merged columns: {', '.join(group)} into { '/'.join(group) }")

    # Build merged DataFrame
    merged = pd.DataFrame(index=df.index)
    for group in groups:
        if len(group) == 1:
            merged[group[0]] = df[group[0]]
        else:
            name = "/".join(group)
            merged[name] = df[group].mean(axis=1)
    return merged


def perform_pqn_normalization(data, sample_labels):
    """
    Performs PQN normalization using the sQCs as the reference spectra.
    Will only return study samples.

    Parameters
    ----------
    data : dataframe
        Wide-format metabolomics data.
    sample_labels : series
        Series containing descriptions of the sample type for all rows in the 
        dataframe.

    Returns
    -------
    final_df : dataframe
        Finalized, normalized, dataframe.
    valid_mask : series
        Mask denoting only study samples.

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
    # Drop features with zero median reference. This one should not matter,
    # as this was handled before. But just to be safe.
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


def detect_pca_outliers(data, labels, var_explained=0.95, conf_level=0.95, res_level=0.95):
    """
    Function that takes the metabolomics-data, and curates it four outliers
    using Hotelling T2-ellipses of PCA loading vectors (and unexplained residuals)

    Parameters
    ----------
    data : dataframe
        Dataframe containing the metabolomics-data.
    labels : series
        Series containing the sample descriptions.
    var_explained : float, optional
        Generate principal components until % variance explained. 
        The default is 0.95 (95%).
    conf_level : float, optional
        The statistic/threshold to use for the T2 ellipse. The default is 0.95.
    res_level : float, optional
        threshold to use for residual filtering. The default is 0.95.

    Returns
    -------
    data.loc[mask_inlier]: dataframe
        Data-output without outliers.
    labels.loc[mask_inlier]: dataframe
        Series denoting which samples were kept.
    mask_inlier : series
        boolean mask for outlier/inliers.
    outlier_info : dataframe
        Outlier report, along with stats.

    """
    
    # 1) Impute & scale
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

    # 4) get something similar to DModX via chi² (QCOmics...))
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

def assign_exp_group(row):
    """
    Update only rows with MAT1 or MAT2 in the material tag
    """
    if row['well'] == "MAT1":
        return "Blank"
    elif row['well'] == "MAT2":
        return "QC"
    else:
        # Leave the existing value unchanged.
        return row["Experimental group"]




def do_rf_imputation(df):
    """
    Imputes missing values with random forest.
    """

    # Do imputation
    imputed = df.copy()
   
    # Select numeric columns for imputation
    numeric_cols = imputed.select_dtypes(include=[np.number]).columns
    
    imputer = IterativeImputer(
        estimator=RandomForestRegressor(n_estimators=100, random_state=42),
        max_iter=25,
        initial_strategy="mean",
        imputation_order="ascending",
        skip_complete=True,
        random_state=0
    )
    
    # Perform imputation on numeric columns
    imputed[numeric_cols] = imputer.fit_transform(imputed[numeric_cols])
    
    return imputed



def global_preprocess(merged_df, exp_path, params):
    """
    Identical to your original pipeline. The only addition:
    - We explicitly preserve 'ReplicateName' as a separate column
      (so it does NOT get treated as numeric).
    """

    # Drop columns that are entirely NaN
    merged_df = merged_df.dropna(axis=1, how='all')
    
    # Subset to QC rows and drop features with more than 40% missing values in QCs
    qc_df = merged_df[merged_df['Experimental group'] == 'QC']
    frac_bad = qc_df[qc_df.replace(0, np.nan).select_dtypes(include=[np.number]).columns].isna().mean()
    cols_to_drop = frac_bad[frac_bad > 0.4].index.tolist()
    merged_df = merged_df.drop(columns=cols_to_drop)

    # Extract numeric features and sample labels
    drop_cols = ["Experimental group", "ReplicateName",'Total Ion Current Area']
    features = merged_df.drop(columns=drop_cols, errors="ignore")
    numeric_features = features.select_dtypes(include=[np.number])
    sample_labels = merged_df["Experimental group"]

   
    # 1) Blank-based filtering
    blank_filtered = apply_blank_threshold(numeric_features, 
                                           sample_labels,
                                           params['sn_ratio'])
    
    # 2) Remove blank samples from further analysis
    non_blank_mask = sample_labels != "Blank"
    blank_filtered = blank_filtered.loc[non_blank_mask]
    sample_labels = sample_labels.loc[non_blank_mask]
    replicate_series = merged_df.loc[non_blank_mask, "ReplicateName"]

    # 3) Drop columns with too many NaNs
    blank_filtered = drop_high_nan_columns(blank_filtered, 
                                           threshold=params['missing_threshold'])
    
    # 3.5) # impute with LoD/2 for peaks that are detected, 
    # but set at 0 (due to skyline pre-processing). 
    # Important to do this AFTER missingness removal.
    blank_filtered = blank_filtered.apply(lambda col: col.replace(0, col[col != 0].min() / 2) if (col != 0).any() else col)

    # 4) Select the best adducts
    reduced_features = select_lowest_variation_qc(blank_filtered, 
                                                  sample_labels)
    
    # 4.5) # Merge the metabolites which we cant reliably separate with our resolution
    reduced_features = merge_columns_by_value_agreement(reduced_features, match_thresh = 2/3)

    # 4.75) Remove samples with more than % missing, 
    # to avoid super low concentration samples biasing our analysis
    mask = reduced_features.isnull().mean(axis=1) <= params['missing_per_sample']
    reduced_features = reduced_features.loc[mask]
    sample_labels = sample_labels.loc[mask]
    
    
    # 5) PQN normalization with QC quality check
    pqn_arr, valid_mask = perform_pqn_normalization(reduced_features, sample_labels)    
    filtered_index = reduced_features.index[valid_mask]
    pqn_df = pd.DataFrame(pqn_arr, columns=reduced_features.columns, index=filtered_index)
    
    # 6) Outlier removal with PCA scores and Hotelling T2 ellipses
    cleaned_features, cleaned_labels, valid_mask, outlier_info = detect_pca_outliers(
                            pqn_df, sample_labels.loc[filtered_index],
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


# Activate pandas<->R conversion
pandas2ri.activate()

# Load R packages
mixOmics = importr('mixOmics')
base     = importr('base')

def load_data(exp_path, ms_filename):
    """
    Loads data and matches it with its growth dynamics.

    Parameters
    ----------
    exp_path : str
        Path to experiment/investigation.
    ms_filename : str
        name of metabolomics data output file.

    Returns
    -------
    df : dataframe
        dataframe with all selected features and metrics.
    feats : TYPE
        list containing the metabolite features.

    """
    met_file = exp_path / f'results/metabolomics/processed/{ms_filename}'
    design_file = exp_path / 'protocol/plate_layout/design_table.tsv'
    #met_file = os.path.join(exp_path, 'results', 'metabolomics', 'processed', ms_filename)
    #design_file = os.path.join(exp_path, 'protocol', 'plate_layout', 'design_table.tsv')
    
    met_df = pd.read_csv(met_file, sep='\t', index_col=0)
    met_df['Well'] = met_df.index.str.split('-').str[2].str.split('.').str[0]
    
    design_df = pd.read_csv(design_file, sep='\t')
    df = met_df.merge(design_df, on='Well', how='left')
    
    growth_file = exp_path / 'results/growth/processed/growth_parameters.tsv'
    #growth_file = os.path.join(exp_path, 'results', 'growth', 'processed', 'growth_parameters.tsv')
    growth_df = pd.read_csv(growth_file, sep='\t', index_col=0).reset_index().rename(columns={'index':'Well'})
    df = df.merge(growth_df, on='Well', how='inner')
    
    feats = [c for c in met_df.columns if c not in ['Experimental group','Well']]
    
    return df, feats

def nested_cv_regression(X, y, inner_splits = 10, n_splits = 3):
    """
    Nested CV regression with ElasticNetCV. Returns mean R2 and coefficient summary DataFrame.
    """
    from sklearn.linear_model import RidgeCV, LassoCV
    outer = KFold(n_splits=n_splits, shuffle=True, random_state=0)
    r2_list, coef_list = [], []
    true_vals, pred_vals = [], []
    
    for tr, te in outer.split(X):
        
        X_tr, X_te = X.iloc[tr,:], X.iloc[te,:]
        y_tr, y_te = y.iloc[tr], y.iloc[te]
        
        scaler = StandardScaler().fit(X_tr)
        Xtr_s, Xte_s = scaler.transform(X_tr), scaler.transform(X_te)
        
        model = ElasticNetCV(cv=inner_splits, random_state=0, l1_ratio=[.1, .3, .5, .7, 1])
        model.fit(Xtr_s, y_tr)
        
        y_pred = model.predict(Xte_s)
        
        r2_list.append(r2_score(y_te, model.predict(Xte_s)))
        coef_list.append(model.coef_)
        true_vals.extend(y_te)
        pred_vals.extend(y_pred)
        
    mean_r2 = float(np.mean(r2_list))
    coef_arr = np.vstack(coef_list)
    coef_df = pd.DataFrame(coef_arr, columns=X.columns)
    summary = coef_df.agg(['mean','std']).T.rename(columns={'mean':'Mean_Coefficient','std':'Std_Coefficient'})
    predictions_df = pd.DataFrame({'True': true_vals, 'Predicted': pred_vals})
    return mean_r2, r2_list, summary, predictions_df


def univariate_stats(df, feats, contrasts):
    """
    Do a mann-whitney U test for the metabolites for all the specified 
    contrasts.

    Parameters
    ----------
    df : dataframe
        annotated dataframe with metabolomics-data.
    feats : list
        metabolite feature names.
    contrasts : dict
        dict of pairwise comparisons.

    Returns
    -------
    out : dataframe 
        dataframe with all stats.

    """
    
    results = []
    for spec1, spec2 in contrasts:
        # filter each subgroup:
        g1 = df.copy()
        for col, val in spec1.items():
            g1 = g1[g1[col].isna() if val is None else g1[col] == val]
        g2 = df.copy()
        for col, val in spec2.items():
            g2 = g2[g2[col].isna() if val is None else g2[col] == val]


        label1 = ", ".join(f"{k}={v}" for k, v in spec1.items())
        label2 = ", ".join(f"{k}={v}" for k, v in spec2.items())

        for f in feats:
            
            # compute U and p
            u_stat, p_raw = mannwhitneyu(
                g1[f], g2[f], alternative="two-sided", nan_policy="omit"
            )
            
            # compute medians
            m1 = g1[f].median()
            m2 = g2[f].median()

            results.append({
                "Group1":      label1,
                "Group2":      label2,
                "Feature":     f,
                "U_stat":      u_stat,
                "p_raw":       p_raw,
                "Median1":     m1,
                "Median2":     m2,
            })

    out = pd.DataFrame(results)
    
    # FDR‐correct across all tests using benjamini-hochberg
    rej, p_corr, _, _ = multipletests(out["p_raw"], alpha=0.05, method="fdr_bh")
    out["p_corrected"] = p_corr
    return out



def train_plsda(X_tr ,y_tr ,n_components = 5):
    """
    Trains a mixOmics PLS-DA model and returns loadings and scores

    Parameters
    ----------
    X_tr : dataframe
    y_tr : series
    n_components : integer, optional
        How many components to use. The default is 5.

    Returns
    -------
    plsda_mod : 
        the trained PLS-DA object.
    scores_df : TYPE
        dataframe with scoring vectors.
    loadings_df : TYPE
        dataframe with loading vectors.

    """
   
    # 1) Convert training data to R objects
    with localconverter(default_converter + pandas2ri.converter):
        rX = pandas2ri.py2rpy(X_tr)
    rY = base.factor(StrVector(y_tr.astype(str).tolist()))
    
    # 2) Fit PLS-DA
    plsda_mod = mixOmics.plsda(X=rX, Y=rY, ncomp=n_components, scale = True)
    
    # 3) Extract variate scores (X‑scores) and loadings from the R model
    scores_df = pd.DataFrame(plsda_mod.rx2('variates').rx2('X'))
    loadings_df = pd.DataFrame(plsda_mod.rx2('loadings').rx2('X'))
    
    comps = [f'PLS{i+1}' for i in range(n_components)]
    scores_df.columns = comps
    scores_df.index=X_tr.index
    
    loadings_df.columns = comps
    loadings_df.index=X_tr.columns
    
    return plsda_mod, scores_df, loadings_df

def predict_plsda(plsda_mod, X_new, n_components):
    """
    Given a trained PLS-DA model, predict new samples.

    Parameters
    ----------
    plsda_mod : 
        the trained PLS-DA object.
    X_new : dataframe
    n_components : integer

    Returns
    -------
    final_preds : array
        Predictions.

    """
    
    # 1) Convert new data to R
    with localconverter(default_converter + pandas2ri.converter):
        rX_new = pandas2ri.py2rpy(X_new)
        
    # 2) Call the S3 generic predict() from mixOmics
    pr = robjects.r['predict'](plsda_mod, rX_new, method='max.dist')
    
    class_vec = pr.rx2('class').rx2('max.dist')
    
    # Convert to a Python list (or ndarray) and reshape:
    flat = list(class_vec)                         
    n_samples = X_new.shape[0]
    ncomp     = n_components       
    arr       = np.array(flat).reshape((n_samples, ncomp), order='F')
    
    # Wrap in a DataFrame:
    cols     = [f'PLS{i+1}' for i in range(ncomp)]
    class_mat = pd.DataFrame(arr, index=X_new.index, columns=cols)
    
    final_preds = class_mat.iloc[:, -1]
    
    return final_preds


def nested_cv_regression_repeated(X, y, inner_splits=10, outer_splits=3,
                                  outer_repeats=10, random_state=0):
    """
    Performs repeated cross-validation to estimate generalization performance
    for treated-only samples. Returns estimated coefficients and std.

    Parameters
    ----------
    X : dataframe
        Dataframe containing the independent variables.
    y : series
        series containing the dependent variable.
    inner_splits : integer, optional
        CV-splits to estimate alpha in the nested loop. The default is 10.
    outer_splits : integer, optional
        K-folds. The default is 3.
    outer_repeats : integer, optional
        How many times to repeat the cv-sampling. The default is 10.
    random_state : integer, optional

    Returns
    -------
    mean_r2 : float
    r2_list : list
    coef_summary : dataframe
        Dataframe containing estimated coefficients across all repeats.
    predictions_df : dataframe
        Dataframe containing all the outputs.

    """
    
    # Set up repeated outer CV
    outer_cv = RepeatedKFold(
        n_splits=outer_splits,
        n_repeats=outer_repeats,
        random_state=random_state
    )

    r2_list = []
    coef_list = []
    true_vals, pred_vals = [], []

    # Outer loop over repeated splits
    for train_idx, test_idx in outer_cv.split(X):
        X_tr, X_te = X.iloc[train_idx, :], X.iloc[test_idx, :]
        y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]

        # Standardize
        scaler = StandardScaler().fit(X_tr)
        Xtr_s = scaler.transform(X_tr)
        Xte_s = scaler.transform(X_te)

        # Inner model selection
        model = ElasticNetCV(
            cv=inner_splits,
            random_state=random_state,
            l1_ratio = np.linspace(0.05, 1.0, 20, endpoint = False),
        )
        model.fit(Xtr_s, y_tr)

        # Predictions
        y_pred = model.predict(Xte_s)
        r2_fold = r2_score(y_te, y_pred)

        # Collect results
        r2_list.append(r2_fold)
        coef_list.append(model.coef_)
        true_vals.extend(y_te)
        pred_vals.extend(y_pred)

    # Collect metrics
    mean_r2 = float(np.mean(r2_list))
    coef_arr = np.vstack(coef_list)
    coef_df = pd.DataFrame(coef_arr, columns=X.columns)
    mean_series = coef_df.mean()
    sd_series = coef_df.std()
    coef_summary = pd.DataFrame({
        'Mean_Coefficient': mean_series,
        'SD_Coefficient': sd_series
    })

    # Predictions DataFrame
    predictions_df = pd.DataFrame({'True': true_vals, 'Predicted': pred_vals})

    return mean_r2, r2_list, coef_summary, predictions_df
