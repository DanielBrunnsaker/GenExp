import argparse
import numpy as np
import pandas as pd
from pathlib import Path

from sklearn.model_selection import KFold, RepeatedKFold, RepeatedStratifiedKFold
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
from rpy2.robjects.packages import importr
from sklearn.metrics import classification_report
from collections import defaultdict


from utils.metabolomics_utils import *
import config

# Suppress warnings
simplefilter("ignore", category=ConvergenceWarning)

# Activate pandas<->R conversion
pandas2ri.activate()

# Load R packages
mixOmics = importr('mixOmics')
base     = importr('base')

def met_analysis(exp_path, ms_filename):

    # Contrasts for the univariate testing
    contrasts = [({'Treatment': 'None', 'Supplement': 'None'},
      {'Treatment': 'None', 'Supplement': 'PosLow'}),
     ({'Treatment': 'None', 'Supplement': 'None'},
      {'Treatment': 'None', 'Supplement': 'PosHigh'}),
     ({'Treatment': 'None', 'Supplement': 'None'},
      {'Treatment': 'None', 'Supplement': 'NegHigh'}),
     ({'Treatment': 'Yes', 'Supplement': 'None'},
      {'Treatment': 'Yes', 'Supplement': 'PosLow'}),
     ({'Treatment': 'Yes', 'Supplement': 'None'},
      {'Treatment': 'Yes', 'Supplement': 'PosHigh'}),
     ({'Treatment': 'Yes', 'Supplement': 'None'},
      {'Treatment': 'Yes', 'Supplement': 'NegHigh'}),
     ({'Treatment': 'None', 'Supplement': 'None'}, 
      {'Treatment': 'Yes', 'Supplement': 'None'})]


    ms_filename = 'ms_output_imputed.tsv'


    # decide on components, move this to config?
    n_components = 5
        
    EXPERIMENT_DIR = Path(exp_path)
    print('\nAnalysing metabolomics data...\n')
    
    # Load data
    df, feats = load_data(EXPERIMENT_DIR, ms_filename)
    
    # Regression on treated, estimate variables important for resistance
    mask = df['Treatment'].str.lower()=='yes'
    Xr = df.loc[mask, feats]
    yr = df.loc[mask, 'AUC']
   
    mean_r2, r2_list, coef_sum, predictions_df = nested_cv_regression_repeated(
        Xr, yr,
        inner_splits=5,
        outer_splits=3,
        outer_repeats=10,
        random_state=0
    )
    
    print(f'Resistance to treatment, R²: {mean_r2:.2f}')
    coef_sum.to_csv(EXPERIMENT_DIR / 'results/metabolomics/coefficients/enet_coef_summary.csv')
    print('Regression results on resistance saved in results/metabolomics/coefficients/supp_stats.csv...\n')
    
    # Multiclass PLS-DA, to estimate metabolic separability
    Xp, yp = df[feats], df['Experimental group']
    _, scores_df, loadings_df = train_plsda(Xp, yp, n_components=config.PLS_DA_COMPONENTS)
    
    
    # Performance evaluation for PLS-DA, repeated Kfold
    rkf = RepeatedStratifiedKFold(n_splits=3, n_repeats=10, random_state=0)
    all_f1s = []
    for tr, te in rkf.split(Xp, yp):
        X_tr, X_te = Xp.iloc[tr], Xp.iloc[te]
        y_tr, y_te = yp.iloc[tr], yp.iloc[te]
        
        plsda_model, _, _ = train_plsda(X_tr, y_tr, n_components = config.PLS_DA_COMPONENTS)
        predictions = predict_plsda(plsda_model, X_te, n_components = config.PLS_DA_COMPONENTS)

        # Combine with true labels:
        preds_df = pd.DataFrame({
            'TrueClass':      y_te.values,
            'PredictedClass': predictions
        }, index=X_te.index)

        report = classification_report(
            preds_df['TrueClass'],
            preds_df['PredictedClass'],
            output_dict=True, 
            zero_division=0  
        )
        
        f1 = report['macro avg']['f1-score']
        all_f1s.append(f1)
    

    avg_f1_macro = np.mean(all_f1s)
    print(f"Average Macro F₁ (metabolic separability) = {avg_f1_macro:.2f}")
    
    # Save outputs
    loadings_df.to_csv(EXPERIMENT_DIR / 'results/metabolomics/coefficients/plsda_loadings.csv')
    scores_df.to_csv(EXPERIMENT_DIR / 'results/metabolomics/predictions/plsda_scores.tsv', sep='\t')
    
    # Do some univariate stat with mann-whitney
    uni = univariate_stats(df, feats=feats, contrasts=contrasts)
    uni.to_csv(EXPERIMENT_DIR / 'results/metabolomics/coefficients/supp_stats.csv', index=False, sep = '\t')
    
    print('Univariate statistics saved in results/metabolomics/coefficients/supp_stats.csv...')
    
    # Save the data, just in case
    df.to_csv(EXPERIMENT_DIR / 'results/metabolomics/processed/data_with_annotation.tsv', index=False, sep = '\t')
    
    


if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Multi and univariate analysis')
    parser.add_argument('--exp_path', required=True)
    args = parser.parse_args()
    met_analysis(args.exp_path, 'ms_output_imputed.tsv')