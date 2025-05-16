#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri May 16 11:10:05 2025

@author: danbru
"""

import os
import glob
import argparse
import pandas as pd
import umap
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import numpy as np
from pycombat import Combat



def load_and_concatenate(tsv_dir, join_type):
    """
    Reads all TSV files in `tsv_dir`, aligns their columns based on `join_type`, and concatenates them.
    """
    pattern = os.path.join(tsv_dir, '*.tsv')
    files = sorted(glob.glob(pattern))
    if not files:
        raise FileNotFoundError(f"No TSV files found in {tsv_dir}")

    dfs = []
    for f in files:
        df = pd.read_csv(f, sep='\t')
        df['source_file'] = os.path.basename(f)
        dfs.append(df)

    if join_type == 'inner':
        common = set(dfs[0].columns)
        for df in dfs[1:]:
            common &= set(df.columns)
        common = sorted(common)
        dfs = [df[common] for df in dfs]
        concatenated = pd.concat(dfs, ignore_index=True)
    else:
        concatenated = pd.concat(dfs, ignore_index=True, sort=False)

    return concatenated

def prepare_data(data, impute=True):
    """
    Select numeric columns and optionally impute missing values.
    """
    numeric = data.select_dtypes(include=['number'])
    if numeric.shape[1] == 0:
        raise ValueError("No numeric columns found for UMAP embedding")
    if impute:
        numeric = numeric.fillna(numeric.mean())
    return numeric

def apply_combat(numeric_df, batch_series):
    """
    Uses pycombat to correct batch effects on numeric_df.
    """
    combat = Combat()
    # pycombat expects input as numpy array; batch_series as array-like
    corrected_array = combat.fit_transform(numeric_df.values, batch_series.values)
    return pd.DataFrame(corrected_array, index=numeric_df.index, columns=numeric_df.columns)


def run_umap(data, args):
    """
    Processes data, applies optional batch correction, then computes UMAP embedding.
    """
    # Extract numeric data and impute
    numeric = prepare_data(data)

    # Batch correction
    if args.combat:
        if args.batch_col not in data.columns:
            raise KeyError(f"Batch column '{args.batch_col}' not found in DataFrame.")
        batch = data[args.batch_col].astype(str)
        numeric = apply_combat(numeric, batch)

    reducer = umap.UMAP(
        n_neighbors=args.n_neighbors,
        min_dist=args.min_dist,
        metric=args.metric
    )
    embedding = reducer.fit_transform(numeric)
    return embedding


def plot_umap(data, embedding, output_file='umap_plot.png'):
    """
    Plots a 2D UMAP embedding using 'Supplement' for color and 'Treatment' for marker shape.
    Adds legends outside to the right.
    """
    df_emb = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'], index=data.index)
    df_emb['Supplement'] = data['Supplement'].astype(str)
    df_emb['Treatment'] = data['Treatment'].astype(str)

    supplements = np.unique(df_emb['Supplement'])
    treatments = np.unique(df_emb['Treatment'])

    # Colors for supplements
    cmap = cm.get_cmap('tab10', len(supplements))
    color_dict = {s: cmap(i) for i, s in enumerate(supplements)}

    # Marker styles for treatments
    markers = ['o', 's', '^', 'P', '*', 'X', 'D', 'v', '<', '>']
    if len(treatments) > len(markers):
        raise ValueError(f"Too many Treatment categories: {len(treatments)} > {len(markers)} markers available.")
    marker_dict = {t: markers[i] for i, t in enumerate(treatments)}

    fig, ax = plt.subplots(figsize=(10, 8))
    for s in supplements:
        for t in treatments:
            subset = df_emb[(df_emb['Supplement'] == s) & (df_emb['Treatment'] == t)]
            if subset.empty:
                continue
            ax.scatter(
                subset['UMAP1'], subset['UMAP2'],
                color=color_dict[s], marker=marker_dict[t],
                s=30, alpha=0.8
            )

    ax.set_title('UMAP projection')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')

    from matplotlib.lines import Line2D
    color_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color_dict[s], markersize=8, label=s) for s in supplements]
    marker_handles = [Line2D([0], [0], marker=marker_dict[t], color='k', linestyle='None', markersize=8, label=t) for t in treatments]

    # Make room on the right for legends
    plt.subplots_adjust(right=0.75)
    legend1 = ax.legend(handles=color_handles, title='Supplement', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, fontsize='small')
    ax.add_artist(legend1)
    legend2 = ax.legend(handles=marker_handles, title='Treatment', bbox_to_anchor=(1.02, 0.6), loc='upper left', borderaxespad=0, fontsize='small')

    plt.tight_layout()
    plt.show()
    #plt.savefig(output_file, dpi=300)
    #plt.close()
    
    
    
data = load_and_concatenate('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/unified_ms_data', 'inner')

data.to_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/unified_ms_data/unified_data/unified_data.csv')




data = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/unified_ms_data/unified_data/adjData.tsv', sep ='\t', index_col = 0)

def prepare_data(data, impute=True):
    """
    Select numeric columns and optionally impute missing values.
    """
    numeric = data.select_dtypes(include=['number'])
    if numeric.shape[1] == 0:
        raise ValueError("No numeric columns found for UMAP embedding")
    return numeric.fillna(numeric.mean()) if impute else numeric


def apply_combat(numeric_df, batch_series):
    """
    Uses pycombat to correct batch effects on numeric_df.
    """
    combat = Combat()
    corrected = combat.fit_transform(numeric_df.values, batch_series.values)
    return pd.DataFrame(corrected, index=numeric_df.index, columns=numeric_df.columns)


def compute_umap(numeric_df, n_neighbors, min_dist, metric):
    """
    Given numeric DataFrame, return UMAP embedding.
    """
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric)
    return reducer.fit_transform(numeric_df)


def get_color_dict(categories, cmap_name='tab20'):
    """
    Map each category to a distinct color using a ListedColormap of length len(categories).
    """
    # Use a colormap with enough distinct entries
    cmap = cm.get_cmap(cmap_name, len(categories))
    return {cat: cmap(i) for i, cat in enumerate(categories)}


def get_marker_dict(categories, marker_list=None):
    if marker_list is None:
        marker_list = ['o', 's', '^', 'P', '*', 'X', 'D', 'v', '<', '>']
    if len(categories) > len(marker_list):
        raise ValueError(f"Too many categories: {len(categories)} > {len(marker_list)} markers available.")
    return {cat: marker_list[i] for i, cat in enumerate(categories)}


def plot_umap(data, embedding, args):
    df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'], index=data.index)
    df['Supplement'] = data['Supplement'].astype(str)
    df['Treatment'] = data['Treatment'].astype(str)

    supplements = np.unique(df['Supplement'])
    treatments = np.unique(df['Treatment'])

    color_dict = get_color_dict(supplements, cmap_name=args.cmap)
    marker_dict = get_marker_dict(treatments)

    fig, ax = plt.subplots(figsize=(10, 8))
    for supp in supplements:
        for treat in treatments:
            subset = df[(df['Supplement'] == supp) & (df['Treatment'] == treat)]
            if subset.empty: continue
            ax.scatter(
                subset['UMAP1'], subset['UMAP2'],
                color=color_dict[supp], marker=marker_dict[treat],
                s=30, alpha=0.8
            )

    ax.set_title('UMAP projection')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    _add_legends(ax, color_dict, marker_dict)

    plt.tight_layout()
    #plt.savefig(args.output, dpi=300)
    #plt.close()
    plt.show()

def _add_legends(ax, color_dict, marker_dict):
    """
    Add separate legends for colors and markers to the right of the plot.
    """
    from matplotlib.lines import Line2D
    color_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=cat)
                     for cat, color in color_dict.items()]
    marker_handles = [Line2D([0], [0], marker=marker, color='k', linestyle='None', markersize=8, label=cat)
                      for cat, marker in marker_dict.items()]

    plt.subplots_adjust(right=0.75)
    leg1 = ax.legend(handles=color_handles, title='Supplement', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.add_artist(leg1)
    ax.legend(handles=marker_handles, title='Treatment', bbox_to_anchor=(1.02, 0.6), loc='upper left', borderaxespad=0)


def main():
    # data assumed pre-loaded
    data = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/unified_ms_data/unified_data/adjData.tsv', sep ='\t', index_col = 0)
    data = data.drop(columns = ['AUC','MaxOD', 'mu','Unnamed: 0','Summary_y'])
    data = data[data['Treatment'] == 'None']

    numeric = prepare_data(data)
    batch = data['source_file'].astype(str)
    numeric = apply_combat(numeric, batch)
    
    embedding = compute_umap(numeric, 15, 0.1,'cosine')
    plot_umap(data, embedding, output_file='')    
    

if __name__ == '__main__':
    main()



def prepare_data(data, impute=True):
    numeric = data.select_dtypes(include=['number'])
    if numeric.shape[1] == 0:
        raise ValueError("No numeric columns found for UMAP embedding")
    return numeric.fillna(numeric.mean()) if impute else numeric


def apply_combat(numeric_df, batch_series):
    combat = Combat()
    corrected = combat.fit_transform(numeric_df.values, batch_series.values)
    return pd.DataFrame(corrected, index=numeric_df.index, columns=numeric_df.columns)


def compute_umap(numeric_df, n_neighbors, min_dist, metric):
    reducer = umap.UMAP(n_neighbors=n_neighbors, min_dist=min_dist, metric=metric)
    return reducer.fit_transform(numeric_df)


def get_color_dict(categories, cmap_name='tab20'):
    """
    Map each category to a distinct color using a ListedColormap of length len(categories).
    """
    # Use a colormap with enough distinct entries
    cmap = cm.get_cmap(cmap_name, len(categories))
    return {cat: cmap(i) for i, cat in enumerate(categories)}


def get_marker_dict(categories, marker_list=None):
    if marker_list is None:
        marker_list = ['o', 's', '^', 'P', '*', 'X', 'D', 'v', '<', '>']
    if len(categories) > len(marker_list):
        raise ValueError(f"Too many categories: {len(categories)} > {len(marker_list)} markers available.")
    return {cat: marker_list[i] for i, cat in enumerate(categories)}


def plot_umap(data, embedding):
    df = pd.DataFrame(embedding, columns=['UMAP1', 'UMAP2'], index=data.index)
    df['Supplement'] = data['Supplement'].astype(str)
    df['Treatment'] = data['Treatment'].astype(str)

    supplements = np.unique(df['Supplement'])
    treatments = np.unique(df['Treatment'])

    color_dict = get_color_dict(supplements)
    marker_dict = get_marker_dict(treatments)

    fig, ax = plt.subplots(figsize=(10, 8))
    for supp in supplements:
        for treat in treatments:
            subset = df[(df['Supplement'] == supp) & (df['Treatment'] == treat)]
            if subset.empty: continue
            ax.scatter(
                subset['UMAP1'], subset['UMAP2'],
                color=color_dict[supp], marker=marker_dict[treat],
                s=30, alpha=0.8
            )

    ax.set_title('UMAP projection')
    ax.set_xlabel('UMAP1')
    ax.set_ylabel('UMAP2')
    _add_legends(ax, color_dict, marker_dict)

    plt.tight_layout()
    #plt.savefig(args.output, dpi=300)
    #plt.close()
    plt.show()

def _add_legends(ax, color_dict, marker_dict):
    from matplotlib.lines import Line2D
    color_handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=color, markersize=8, label=cat)
                     for cat, color in color_dict.items()]
    marker_handles = [Line2D([0], [0], marker=marker, color='k', linestyle='None', markersize=8, label=cat)
                      for cat, marker in marker_dict.items()]

    plt.subplots_adjust(right=0.75)
    leg1 = ax.legend(handles=color_handles, title='Supplement', bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
    ax.add_artist(leg1)
    ax.legend(handles=marker_handles, title='Treatment', bbox_to_anchor=(1.02, 0.6), loc='upper left', borderaxespad=0)


def main():
    args = parse_args()
    
    data = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/unified_ms_data/unified_data/adjData.tsv', sep ='\t', index_col = 0)
    data = data.drop(columns = ['AUC','MaxOD', 'Unnamed: 0','Summary_y'])
    data = data[data['Treatment'] == 'None']

    # try to predict growth?
    import numpy as np
    from sklearn.model_selection import KFold, cross_val_score, cross_val_predict
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score
    from sklearn.linear_model import LinearRegression

    import numpy as np
    import pandas as pd
    
    from sklearn.model_selection import KFold
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.linear_model import LinearRegression
    from sklearn.metrics import r2_score
    import shap
    from sklearn.preprocessing import StandardScaler
    
    def residualize_batch(X_train, X_test, batch_train, batch_test):
        """
        Remove batch effects by regressing each feature on batch (train-only),
        then apply the same correction to the test set.
        """
        B_train = pd.get_dummies(batch_train, drop_first=True)
        B_test = pd.get_dummies(batch_test, drop_first=True).reindex(
            columns=B_train.columns, fill_value=0
        )
    
        X_train_corr = np.zeros_like(X_train.values, dtype=float)
        X_test_corr  = np.zeros_like(X_test.values,  dtype=float)
        
        for i in range(X_train.shape[1]):
            y_feat = X_train.iloc[:, i].values
            lm = LinearRegression().fit(B_train.values, y_feat)
            X_train_corr[:, i] = y_feat - lm.predict(B_train.values)
            X_test_corr[:, i]  = X_test.iloc[:, i].values - lm.predict(B_test.values)
        
        return (
            pd.DataFrame(X_train_corr, columns=X_train.columns, index=X_train.index),
            pd.DataFrame(X_test_corr,  columns=X_test.columns,  index=X_test.index)
        )

    
    data = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/unified_ms_data/unified_data/adjData.tsv', sep ='\t', index_col = 0)
    data = data.drop(columns = ['AUC','MaxOD', 'Unnamed: 0','Summary_y'])
    data = data[data['Treatment'] == 'None']
    #data = data[data['Supplement'] == 'None']
    
    y = data['mu']
    
    batch = data['source_file'].astype(str)
    X = prepare_data(data.drop(columns=['mu']))
    X.columns = X.columns.str.split('\[\[').str[0]
    
    model = RandomForestRegressor(n_estimators=500, random_state=42)
    #model = ElasticNetCV(random_state = 0, l1_ratio = np.linspace(0.05, 1, 10))
    cv    = KFold(n_splits=10, shuffle=True, random_state=42)
    
    scores      = []
    y_pred      = np.zeros(len(y))
    importances = np.zeros(X.shape[1])
    
    # 1) CV loop with in-fold residualization
    for train_idx, test_idx in cv.split(X):
        X_tr, X_te = X.iloc[train_idx], X.iloc[test_idx]
        y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]
        b_tr, b_te = batch.iloc[train_idx], batch.iloc[test_idx]
    
        X_tr_corr, X_te_corr = residualize_batch(X_tr, X_te, b_tr, b_te)
        scaler = StandardScaler()
        
        #X_tr_corr = scaler.fit_transform(X_tr_corr)
        #X_te_corr = scaler.transform(X_te_corr)
    
        model.fit(X_tr_corr, y_tr)
        y_pred[test_idx] = model.predict(X_te_corr)
        scores.append(r2_score(y_te, y_pred[test_idx]))
        importances += model.feature_importances_
        #importances += model.coef_
    
    scores      = np.array(scores)
    overall_r2  = r2_score(y, y_pred)
    importances /= cv.get_n_splits()
    
    # 2) Print CV results & mean tree-based importances
    print("R² per fold:", np.round(scores, 3))
    print("Mean R²:    ", np.round(scores.mean(), 3), "±", np.round(scores.std(), 3))
    print("Overall R²:", np.round(overall_r2, 3))
    
    feat_imp_df = pd.DataFrame({
        "feature":    X.columns,
        "importance": importances
    }).sort_values("importance", ascending=False)
    
    print("\nAverage tree-feature importances:")
    print(feat_imp_df.to_string(index=False))
    
    # 3) Train final model on full residualized data
    X_full_corr, _ = residualize_batch(X, X, batch, batch)
    model_full = RandomForestRegressor(n_estimators=500, random_state=42)
    model_full.fit(X_full_corr, y)
    
    explainer = shap.TreeExplainer(model_full, model_output="raw", feature_perturbation="tree_path_dependent")
    shap_values = explainer.shap_values(X_full_corr)
    shap.summary_plot(shap_values, X_full_corr, plot_type="dot", show=True)
    plt.tight_layout()
    plt.show()
    
    
    
    
    
    
    import os
    import re
    import numpy as np
    import pandas as pd
    import shap
    import matplotlib.pyplot as plt
    
    from sklearn.model_selection import StratifiedKFold
    from sklearn.linear_model import LinearRegression
    from sklearn.preprocessing import StandardScaler
    from sklearn.ensemble import RandomForestRegressor
    from sklearn.metrics import r2_score
    
    os.environ["CUDA_VISIBLE_DEVICES"] = ""  # for SHAP on CPU
    
    # --- helper: in‐fold residualization ---
    def residualize_batch(X_train, X_test, batch_train, batch_test):
        B_train = pd.get_dummies(batch_train, drop_first=True)
        B_test = pd.get_dummies(batch_test,  drop_first=True).reindex(
            columns=B_train.columns, fill_value=0
        )
        X_tr_corr = np.zeros_like(X_train.values, dtype=float)
        X_te_corr = np.zeros_like(X_test.values,  dtype=float)
        for i in range(X_train.shape[1]):
            y_feat = X_train.iloc[:, i].values
            lm = LinearRegression().fit(B_train.values, y_feat)
            X_tr_corr[:, i] = y_feat - lm.predict(B_train.values)
            X_te_corr[:, i] = X_test.iloc[:, i].values - lm.predict(B_test.values)
        return (
            pd.DataFrame(X_tr_corr, columns=X_train.columns, index=X_train.index),
            pd.DataFrame(X_te_corr, columns=X_test.columns, index=X_test.index),
        )
    
    # --- load & pre‐filter ---
    data = pd.read_csv(
        '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/unified_ms_data/unified_data/adjData.tsv',
        sep ='\t', index_col=0
    )
    data = data.drop(columns=['AUC','MaxOD','Unnamed: 0','Summary_y'])
    # include all treatments (not just None)
    # data = data[data['Treatment']=='None']
    
    # --- parse Supplement into dose & amino acid ---
    # e.g. "0.5mM Arginine"
    dose_aa = data['Supplement'].str.strip()
    #data['dose'] = dose_aa.str.extract(r'([0-9\.]+)mM').astype(float)
    data['aa']   = dose_aa.str.extract(r'mM\s*([A-Za-z]+)')
    
    # drop the raw Supplement column if you like
    # data = data.drop(columns=['Supplement'])
    
    # --- response & batch & raw X ---
    y     = data['mu']
    batch = data['source_file'].astype(str)
    
    X_raw = prepare_data(data.drop(columns=['mu']))
    X_raw.columns = X_raw.columns.str.split('\[\[').str[0]
    
    # --- add dose + one-hot aa to features ---
    X = X_raw.copy()
    #X['dose'] = data['dose']
    #X = pd.get_dummies(X.join(data['aa']), columns=['aa'], drop_first=True)
    #model = XGBRegressor(n_jobs = 1, n_estimators = 100, learning_rate = 0.1)
    #model = ElasticNetCV(random_state = 0, l1_ratio = [.01 ,.1, .5, .9])
    model = RandomForestRegressor(n_estimators=1000, random_state=42)
    
    # Instead of binning dose, just stratify by amino acid
    strata = data['aa'].fillna('None')
    
    # set up a StratifiedKFold on the amino‐acid groups
    cv = StratifiedKFold(n_splits=20, shuffle=True, random_state=42)
    
    scores      = []
    y_pred      = np.zeros(len(y))
    importances = np.zeros(X.shape[1])
    
    for train_idx, test_idx in cv.split(X, strata):
        X_tr, X_te = X.iloc[train_idx], X.iloc[test_idx]
        y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]
        b_tr, b_te = batch.iloc[train_idx], batch.iloc[test_idx]
    
        # in‐fold batch correction
        X_tr_corr, X_te_corr = residualize_batch(X_tr, X_te, b_tr, b_te)
    
        # (optional) scaling
        #scaler = StandardScaler()
        #X_tr_corr = scaler.fit_transform(X_tr_corr)
        #X_te_corr = scaler.transform(X_te_corr)
    
        #X_tr_corr = X_tr
        #X_te_corr = X_te
    
    
        # fit & predict
        model.fit(X_tr_corr, y_tr)
        y_pred[test_idx] = model.predict(X_te_corr)
        scores.append(r2_score(y_te, y_pred[test_idx]))
        #importances += model.feature_importances_
    
    # aggregate metrics
    scores     = np.array(scores)
    overall_r2 = r2_score(y, y_pred)
    #importances /= cv.get_n_splits()
    
    print("R² per fold:", np.round(scores, 3))
    print("Mean R²:    ", np.round(scores.mean(), 3), "±", np.round(scores.std(), 3))
    print("Overall R²:", np.round(overall_r2, 3))
    
    
    
    
    
    
    
    
    
    
    
    
    # --- your data here ---
    # X: array-like of shape (n_samples, n_features)
    # y: array-like of shape (n_samples,)
    
    y = data['mu']
    
    batch = data['source_file'].astype(str)
    X = prepare_data(data.drop(columns=['mu']))
    X = apply_combat(X, batch)

    # set up model and CV splitter
    model = RandomForestRegressor(n_estimators=500, random_state=42)
    cv = KFold(n_splits=10, shuffle=True, random_state=42)
    
    # 1) per‐fold R² scores
    scores = cross_val_score(model, X, y, cv=cv, scoring='r2')
    print("R² per fold:", np.round(scores, 3))
    print("Mean R²:    ", np.round(scores.mean(), 3),
          "±", np.round(scores.std(), 3))
    
    # 2) overall R² via out‐of‐fold predictions
    y_pred = cross_val_predict(model, X, y, cv=cv)
    overall_r2 = r2_score(y, y_pred)
    print("Overall R²:", np.round(overall_r2, 3))
    
    # 3) compute feature importances averaged over folds
    importances = np.zeros(X.shape[1])
    for train_idx, _ in cv.split(X):
        X_train, y_train = X.iloc[train_idx], y.iloc[train_idx] \
                           if hasattr(X, "iloc") else (X[train_idx], y[train_idx])
        model.fit(X_train, y_train)
        importances += model.feature_importances_
    
    importances /= cv.get_n_splits()
    
    # put into a DataFrame (if X is array, label features as 0,1,2…)
    feat_names = X.columns if hasattr(X, "columns") else [f"X{i}" for i in range(X.shape[1])]
    feat_imp_df = pd.DataFrame({
        "feature": feat_names,
        "importance": importances
    }).sort_values("importance", ascending=False)
    
    print("\nAverage feature importances:")
    print(feat_imp_df.to_string(index=False))


    
    numeric = prepare_data(data)
    batch = data['source_file'].astype(str)
    numeric = apply_combat(numeric, batch)
    
    embedding = compute_umap(numeric, 3, 0.05, 'cosine')
    plot_umap(data, embedding)
    print(f"UMAP plot saved to {args.output} (combat={'on' if args.combat else 'off'})")

if __name__ == '__main__':
    main()