

import os
import argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold
from sklearn.metrics import r2_score, f1_score
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
from rpy2 import robjects
from rpy2.robjects import StrVector, IntVector
from rpy2.robjects.packages import importr
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, default_converter

# Suppress warnings
simplefilter("ignore", category=ConvergenceWarning)

# Activate pandas<->R conversion
pandas2ri.activate()

# Load R packages
mixOmics = importr('mixOmics')
base     = importr('base')


def load_data(exp_path: str, ms_filename: str, include_growth: bool = True):
    """
    Load and merge metabolomics, design, (optional) growth data; median-impute missing values.
    Returns merged DataFrame and list of feature columns.
    """
    met_file = os.path.join(exp_path, 'results', 'metabolomics', 'processed', ms_filename)
    design_file = os.path.join(exp_path, 'protocol', 'plate_layout', 'design_table.tsv')
    met_df = pd.read_csv(met_file, sep='\t', index_col=0)
    met_df['Well'] = met_df.index.str.split('-').str[2].str.split('.').str[0]
    design_df = pd.read_csv(design_file, sep='\t')
    df = met_df.merge(design_df, on='Well', how='left')
    if include_growth:
        growth_file = os.path.join(exp_path, 'results', 'growth', 'processed', 'growth_parameters.tsv')
        growth_df = pd.read_csv(growth_file, sep='\t', index_col=0).reset_index().rename(columns={'index':'Well'})
        df = df.merge(growth_df, on='Well', how='inner')
    feats = [c for c in met_df.columns if c not in ['Experimental group','Well']]
    df[feats] = df[feats].fillna(df[feats].median())
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
    Pairwise Mann–Whitney U tests with FDR correction,
    plus group medians and a 'Direction' of effect.
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

        if g1.empty or g2.empty:
            raise ValueError(f"No rows for contrast {spec1} vs {spec2}")

        label1 = ", ".join(f"{k}={v}" for k, v in spec1.items())
        label2 = ", ".join(f"{k}={v}" for k, v in spec2.items())

        for f in feats:
            # compute U and p
            u_stat, p_raw = mannwhitneyu(
                g1[f], g2[f], alternative="two-sided", nan_policy="omit"
            )
            # compute medians and direction
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
    # FDR‐correct across all tests
    rej, p_corr, _, _ = multipletests(out["p_raw"], alpha=0.05, method="fdr_bh")
    out["p_corrected"] = p_corr
    out["Significant"] = rej
    return out



def train_plsda(X_tr ,y_tr ,n_components = 5):
    """
    Train a classic mixOmics PLS-DA model on (X_tr, y_tr).
    Returns:
      - plsda_mod   : the underlying R plsda object
      - scores_tr_df: DataFrame (n_tr × n_components) of X variates
      - loadings_df : DataFrame (n_features × n_components) of loadings
    """
    # 1) Convert training data to R objects
    with localconverter(default_converter + pandas2ri.converter):
        rX = pandas2ri.py2rpy(X_tr)
    rY = base.factor(StrVector(y_tr.astype(str).tolist()))
    
    # 2) Fit PLS-DA
    #plsda_mod = mixOmics.plsda(X=rX, Y=rY, ncomp=n_components, scale = True)
    plsda_mod = mixOmics.plsda(X=rX, Y=rY, ncomp=n_components, scale = True)
    
    # 3) Extract variate scores (X‑scores) and loadings from the R model
    scores_df = pd.DataFrame(plsda_mod.rx2('variates').rx2('X'))
    loadings_df = pd.DataFrame(plsda_mod.rx2('loadings').rx2('X'))
    
    comps = [f'PLS{i+1}' for i in range(n_components)]
    #scores_tr_df = pd.DataFrame(scores_tr, columns=comps)
    scores_df.columns = comps
    scores_df.index=X_tr.index
    
    loadings_df.columns = comps
    loadings_df.index=X_tr.columns
    
    return plsda_mod, scores_df, loadings_df

def predict_plsda(plsda_mod, X_new, n_components):
    """
    Given a trained mixOmics plsda model, predict new samples in X_new.
    
    Returns:
      - scores_df: pandas DataFrame of the X‑variates (n_new × n_components)
      - class_pred: numpy array of predicted class labels (length n_new)
    """
    # 1) Convert new data to R
    with localconverter(default_converter + pandas2ri.converter):
        rX_new = pandas2ri.py2rpy(X_new)
        
    # 2) Call the S3 generic predict() from mixOmics
    pr = robjects.r['predict'](plsda_mod, rX_new, method='max.dist')
    
    # pr = robjects.r['predict'](plsda_mod, rX_new, method='max.dist')
    class_vec = pr.rx2('class').rx2('max.dist')    # flat StrVector of length n_samples * ncomp
    
    # Convert to a Python list (or ndarray) and reshape:
    flat = list(class_vec)                         # length = n_samples * ncomp
    n_samples = X_new.shape[0]
    ncomp     = n_components       # e.g. 5
    arr       = np.array(flat).reshape((n_samples, ncomp), order='F')
    
    # Wrap in a DataFrame so you can inspect all columns if you like:
    cols     = [f'PLS{i+1}' for i in range(ncomp)]
    class_mat = pd.DataFrame(arr, index=X_new.index, columns=cols)
    
    # Now extract the “final” predictions (last column):
    final_preds = class_mat.iloc[:, -1]
    
    return final_preds

def main(exp_path, ms_filename):
        
    from sklearn.metrics import classification_report
    from collections import defaultdict

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

    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756' # FA
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539' # Caffeine
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618' # Spermine
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407' # Lactic acid
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655' # LiCl
    
    # decide on components, move this to config?
    n_components = 5
    
    
    exp_paths = [
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655',
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407',
    ]
    for exp_path in exp_paths:
        
        
        

    
        # Load data
        df, feats = load_data(exp_path, ms_filename, include_growth=True)
        os.makedirs(os.path.join(exp_path,'results','metabolomics','coefficients'), exist_ok=True)
        os.makedirs(os.path.join(exp_path,'results','metabolomics','predictions'), exist_ok=True)
        results = {}
        
        # Regression on treated
        mask = df['Treatment'].str.lower()=='yes'
        Xr = df.loc[mask, feats]
        yr = df.loc[mask, 'AUC']
        mean_r2, r2_list, coef_sum, pred_df = nested_cv_regression(Xr, yr, 10, 3)
        print(f'Resistance R²: {mean_r2:.2f}')
        
        coef_sum.to_csv(os.path.join(exp_path,'results','metabolomics','coefficients','enet_coef_summary.csv'))
        pred_df.to_csv(os.path.join(exp_path,'results','metabolomics','predictions','enet_coef_preds.csv'))
        results['r2'] = mean_r2
        
        # Multiclass PLS-DA
        Xp, yp = df[feats], df['Experimental group']
        _, scores_df, loadings_df = train_plsda(Xp, yp, n_components=n_components)
        
        # Performance evaluation for PLS-DA
        from sklearn.model_selection import RepeatedStratifiedKFold
        rkf = RepeatedStratifiedKFold(n_splits=2, n_repeats=50, random_state=42)
        
        f1_scores = defaultdict(list)
    
        all_f1s = []
        for tr, te in rkf.split(Xp, yp):
            X_tr, X_te = Xp.iloc[tr], Xp.iloc[te]
            y_tr, y_te = yp.iloc[tr], yp.iloc[te]
            
            
            # Train on training data:
            plsda_model, _, _ = train_plsda(X_tr, y_tr, n_components = n_components)
            predictions = predict_plsda(plsda_model, X_te, n_components = n_components)
            #predictions = classify_scores(scores_df, y_te)
    
            # Combine with true labels:
            preds_df = pd.DataFrame({
                'TrueClass':      y_te.values,
                'PredictedClass': predictions
            }, index=X_te.index)
    
            report = classification_report(
                preds_df['TrueClass'],
                preds_df['PredictedClass'],
                output_dict=True,    # so we can turn it into a DataFrame
                zero_division=0      # avoid NaNs if a class never gets predicted
            )
            
            # accumulate per-class F1
            for cls, metrics in report.items():
                # skip the summary rows
                if cls in ('accuracy','macro avg','weighted avg'):
                    continue
                f1_scores[cls].append(metrics['f1-score'])
            
    
            # Compute macro F1:
            #f1 = f1_score(preds_df['TrueClass'], preds_df['PredictedClass'], average='macro')
            f1 = report['macro avg']['f1-score']
            #print(f"Macro F₁ = {f1:.2f}")
    
            all_f1s.append(f1)
        
        
        # Build a summary DataFrame
        summary = pd.DataFrame([
            {
                'class': cls,
                'mean_f1': np.mean(scores),
                'std_f1':  np.std(scores, ddof=1)
            }
            for cls, scores in f1_scores.items()
        ]).set_index('class')
        
        
        avg_f1_macro = np.mean(all_f1s)
        print(f"Average Macro F₁ = {avg_f1_macro:.2f}")
        
        # Save outputs
        loadings_df.to_csv(os.path.join(exp_path,'results','metabolomics','coefficients','plsda_loadings.csv'))
        scores_df.to_csv(os.path.join(exp_path,'results','metabolomics','predictions','plsda_scores.tsv'), sep='\t')    
        
    
        uni = univariate_stats(df, feats=feats, contrasts=contrasts)
        
        #uni = univariate_stats(df, 'Experimental group', feats)
        uni.to_csv(os.path.join(exp_path,'results','metabolomics','coefficients','supp_stats.csv'), index=False, sep = '\t')
        
        df.to_csv(os.path.join(exp_path,'results','metabolomics','processed','data_with_annotation.tsv'), index=False, sep = '\t')
        
        # Save to CSV
        out_csv = os.path.join(exp_path,
            'results','metabolomics','predictions','plsda_classwise_f1_summary.csv'
        )
        summary.to_csv(out_csv)
    
        # Save a heatmap?
   

if __name__=='__main__':
    parser = argparse.ArgumentParser(description='Metabolomics PLS-DA Analysis')
    parser.add_argument('--exp_path', required=True)
    parser.add_argument('--ms_filename', default='ms_output.tsv')
    args = parser.parse_args()
    main(args.exp_path, args.ms_filename)






import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

def plot_heatmap_sns(
    df,
    data_cols,
    group_cols=('Supplement', 'Treatment'),
    scale=None,             # one of {'zscore', 'minmax', None}
    scale_axis='rows',      # 'rows' to scale features, 'columns' to scale samples
    cmap='viridis',
    figsize=(10, 8),
    row_cluster=False,
    col_cluster=False,
    group_order=None,       # dict e.g. {'Treatment': ['A','B'], 'Supplement': ['X','Y','Z']}
    sort_by=None            # list of metadata columns defining sample order
):
    """
    Seaborn clustermap with colored sample-bars, explicit scaling axis, custom ordering,
    legends for groups, narrow colorbar aligned with heatmap, and full feature labels.

    Parameters
    ----------
    df : pd.DataFrame
        One row per sample, with metadata and feature columns.
    data_cols : list[str]
        Names of numeric feature columns to plot.
    group_cols : tuple(str,str)
        Two metadata columns for color bars and legends.
    scale : {'zscore','minmax',None}
        'zscore' or 'minmax' to standardize; None for raw values.
    scale_axis : {'rows','columns'}
        Axis along which to scale: 'rows' to scale each feature across samples,
        'columns' to scale each sample across features.
    cmap : str
        A matplotlib colormap name.
    figsize : tuple
        Figure size.
    row_cluster, col_cluster : bool
        Whether to cluster rows/columns.
    group_order : dict or None
        If provided, maps each group_col to desired category order.
    sort_by : list or None
        List of metadata columns defining the sorting priority of samples.
    """
    # 1) Copy and apply custom ordering/categories
    df_plot = df.copy()
    if group_order:
        for col, order in group_order.items():
            if col in df_plot:
                df_plot[col] = pd.Categorical(df_plot[col], categories=order, ordered=True)
    # 2) Sort samples
    if sort_by is None:
        sort_by = list(group_cols)
    df_plot = df_plot.sort_values(sort_by)
    # 3) Build data matrix: features as rows, samples as columns
    mat = df_plot[data_cols].T
    # 4) Build lut and col_colors for sample annotations
    lut = {}
    for col in group_cols:
        cats = df_plot[col].cat.categories if hasattr(df_plot[col], 'cat') else sorted(df_plot[col].unique())
        palette = sns.color_palette(None, n_colors=len(cats))
        lut[col] = dict(zip(cats, palette))
    color_rows = [
        [lut[col][df_plot.iloc[i][col]] for col in group_cols]
        for i in range(len(df_plot))
    ]
    col_colors = pd.DataFrame(color_rows, index=df_plot.index, columns=group_cols)
    # 5) Determine scaling args
    z_score = None
    standard_scale = None
    if scale == 'zscore':
        z_score = 0 if scale_axis == 'rows' else 1
    elif scale == 'minmax':
        standard_scale = 0 if scale_axis == 'rows' else 1
    # 6) Plot clustermap with narrow, aligned colorbar
    g = sns.clustermap(
        mat,
        cmap=cmap,
        row_cluster=row_cluster,
        col_cluster=col_cluster,
        col_colors=col_colors,
        figsize=figsize,
        z_score=z_score,
        standard_scale=standard_scale,
        cbar_pos=(0.92, 0.2, 0.015, 0.6),  # x, y, width, height
        cbar_kws={'label': f"{scale} ({scale_axis})" if scale else 'Value'}
    )
    # 7) Add legends for group colors
    for idx, col in enumerate(group_cols):
        handles = [mpatches.Patch(color=color, label=cat)
                   for cat, color in lut[col].items()]
        g.ax_heatmap.legend(
            handles=handles,
            title=col,
            bbox_to_anchor=(1.02 + idx*0.15, 1),
            loc='upper left',
            frameon=False
        )
    # 8) Tidy labels and show all features
    g.ax_heatmap.set_xlabel("Samples")
    g.ax_heatmap.set_ylabel("Features")
    plt.setp(g.ax_heatmap.get_xticklabels(), rotation=90)
    plt.setp(g.ax_heatmap.get_yticklabels(), rotation=0, fontsize=6)
    plt.tight_layout()
    plt.show()

plot_heatmap_sns(
    df,
    df.columns[:-10],
    group_cols=('Supplement','Treatment'),
    scale='zscore',            # or 'minmax', or None
    scale_axis='rows',         # 'rows' or 'columns'
    cmap='coolwarm',
    figsize=(12,6),
    row_cluster=True,
    col_cluster=False,
    group_order={
        'Treatment': ['None','Yes'],
        'Supplement': ['None','PosLow','PosHigh','NegHigh']
    },
    sort_by=['Treatment','Supplement']
)



# 1) Example loadings DataFrame
#    index: metabolite names
#    columns: ['Comp1','Comp2',...]
#loadings_df = pd.read_csv("plsda_loadings.csv", index_col=0)

# 2a) Top-N approach: 
N = 20
# for each metabolite, find its largest abs-loading over all components:
max_abs = loadings_df.abs().max(axis=1)
top_feats = max_abs.nlargest(N).index.tolist()

# 2b) Threshold approach:
threshold = 0.5
mask = (loadings_df.abs() > threshold).any(axis=1)
top_feats = loadings_df.index[mask].tolist()

# 3) Subset your data
#    assume 'df' has metadata + all metabolites as columns
#    and your group_cols are ('Supplement','Treatment')
filtered_df = df[['Supplement','Treatment'] + top_feats]

# 4) Plot
plot_heatmap_sns(
    filtered_df,
    data_cols=top_feats,
    group_cols=('Supplement','Treatment'),
    scale='zscore',
    scale_axis='rows',
    cmap='coolwarm',
    figsize=(12,6),
    row_cluster=False,
    col_cluster=False,
    group_order={
        'Treatment': ['None','Yes'],
        'Supplement': ['None','PosLow','PosHigh','NegHigh']
    },
    sort_by=['Treatment','Supplement']
)
