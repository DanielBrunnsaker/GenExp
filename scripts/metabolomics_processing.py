#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 15 14:51:58 2025

@author: danbru
"""

from pathlib import Path
import pandas as pd

from utils.metabolomics_utils import *
import config


def met_processing(exp_folder):
    
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/partially_completed/glutamine_202504251440' # glutamine_acetate
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/partially_completed/lysine_202504291748' # lysine sucrose
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/partially_completed/aminoadipate_202504291411' # aminoadipate
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756' # FA
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539' # Caffeine
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618' # Spermine
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407' # Lactic acid
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655' # LiCl
    #exp_paths = [
    #    '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618',
    #    '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756',
    #    '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539',
    #    '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655',
    #    '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407',
    #    '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamine_202504251440', # glutamine_acetate
    #    '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/lysine_202504291748', # lysine sucrose
    #    '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/aminoadipate_202504291411' # aminoadipate
    #]
    #for exp_folder in exp_paths:
    
    EXPERIMENT_DIR = Path(exp_folder)
    print("\nProcessing metabolomics data...\n")
    
    # Remember to move these to global config later
    params = config.PARAMS

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
    


if __name__ == '__main__':
    
    import argparse
    parser = argparse.ArgumentParser(description="Metabolomics processing")
    parser.add_argument("--output_folder", required=True, type=str, help="Experiment folder")

    args = parser.parse_args()
    met_processing(args.output_folder)

