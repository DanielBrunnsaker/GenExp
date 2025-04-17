

import os
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, StratifiedKFold
from sklearn.metrics import r2_score, f1_score, roc_auc_score
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning
simplefilter("ignore", category=ConvergenceWarning)


def average_coefficients(coef_list, feature_names):
    """Average coefficient vectors (one per fold) and return a DataFrame of means and standard deviations."""
    
    coefs = np.array(coef_list)
    mean_coef = np.mean(coefs, axis=0)
    std_coef = np.std(coefs, axis=0)
    
    return pd.DataFrame({
        "Feature": feature_names,
        "Mean_Coefficient": mean_coef,
        "Std_Coefficient": std_coef
    }).set_index("Feature")

def pairwise_comparison_df(merged_df, group_col, metabolite_columns, alpha=0.05, correction_method='fdr_bh'):
    """
    Perform pairwise Mann–Whitney U tests for each metabolite between every pair of groups in group_col,
    and apply multiple testing correction.

    Parameters:
        merged_df (pd.DataFrame): DataFrame containing the data.
        group_col (str): Column name in merged_df indicating group membership.
        metabolite_columns (list of str): List of columns in merged_df representing metabolite measurements.
        alpha (float): Significance level.
        correction_method (str): Method for multiple testing correction (e.g., 'fdr_bh' or 'bonferroni').

    Returns:
        pd.DataFrame: DataFrame containing the pairwise test results with raw and corrected p-values.
    """
    tests = []
    groups = merged_df[group_col].unique()

    # Loop over all unique pairs of groups
    for i in range(len(groups)):
        for j in range(i + 1, len(groups)):
            group1, group2 = groups[i], groups[j]
            data1 = merged_df[merged_df[group_col] == group1]
            data2 = merged_df[merged_df[group_col] == group2]
            for feat in metabolite_columns:
                # Perform Mann–Whitney U test with two-sided alternative
                stat, p = mannwhitneyu(data1[feat], data2[feat], alternative='two-sided', nan_policy='omit')
                tests.append({
                    "Group 1": group1,
                    "Group 2": group2,
                    "Metabolite": feat,
                    "U_stat": stat,
                    "raw_p_value": p
                })
    
    # Extract raw p-values for correction
    raw_p_values = [test["raw_p_value"] for test in tests]
    
    # Apply multiple testing correction
    reject, corrected_p, _, _ = multipletests(raw_p_values, alpha=alpha, method=correction_method)
    
    # Attach corrected p-values and significance flags to results
    for i, test in enumerate(tests):
        test["corrected_p_value"] = corrected_p[i]
        test["Significant"] = reject[i]
    
    return pd.DataFrame(tests)

def merge_data(experiment_path, ms_filename):
    """
    Load metabolomics, design, and growth data; merge on "Well"; impute missing values.
    Used for regression and overall classification analyses.
    """
    
    met_file = os.path.join(experiment_path, "results", "metabolomics", "processed", ms_filename)
    design_file = os.path.join(experiment_path, "protocol", "plate_layout", "design_table.tsv")
    growth_file = os.path.join(experiment_path, "results", "growth", "processed", "growth_parameters.tsv")
    
    met_df = pd.read_csv(met_file, sep="\t", index_col=0)
    
    # Extract Well from index 
    met_df["Well"] = met_df.index.str.split("-").str[2].str.split(".").str[0]
    design_df = pd.read_csv(design_file, sep="\t")
    growth_df = pd.read_csv(growth_file, sep="\t", index_col=0).reset_index().rename(columns={'index': 'Well'})
    
    merged_df = pd.merge(met_df, design_df, on="Well", how="left")
    merged_df = pd.merge(merged_df, growth_df, on="Well", how="inner")
    metabolite_columns = [col for col in met_df.columns if col not in ["Experimental group", "Well"]]
    merged_df[metabolite_columns] = merged_df[metabolite_columns].apply(lambda col: col.fillna(col.median()))
    
    return merged_df, metabolite_columns

def load_metabolomics_data_with_design(experiment_path, ms_filename):
    """
    Load metabolomics data from ms_output.tsv and design data from design_table.tsv;
    merge on "Experimental group" and "Summary"; impute missing values.
    Used for pairwise classification on Supplement groups.
    """
    met_file = os.path.join(experiment_path, "results", "metabolomics", "processed", ms_filename)
    design_file = os.path.join(experiment_path, "protocol", "plate_layout", "design_table.tsv")
    met_df = pd.read_csv(met_file, sep="\t", index_col=0)
    design_df = pd.read_csv(design_file, sep="\t")
    
    merged_df = pd.merge(met_df, design_df, left_on="Experimental group", right_on="Summary", how="left")
    metabolite_columns = [col for col in met_df.columns if col != "Experimental group"]
    merged_df = merged_df.drop_duplicates(subset=metabolite_columns)
    #merged_df[metabolite_columns] = merged_df[metabolite_columns].apply(lambda col: col.fillna(col.median()))
    
    return merged_df, metabolite_columns

def nested_cv_regression(X, y, outer_cv_folds=3):
    """
    Nested CV for regression using ElasticNetCV.
    Outer loop splits data; inner tuning is handled by ElasticNetCV.
    Returns outer R² scores, coefficient vectors, and predictions for each fold.
    """
    outer_cv = KFold(n_splits=outer_cv_folds, shuffle=True, random_state=0)
    r2_scores, coef_list, predictions = [], [], []
    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X)):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train_raw)
        X_test = scaler.transform(X_test_raw)
        
        model = ElasticNetCV(cv=5, random_state=0, l1_ratio = [.1, .5, .7, .9, .95, .99, 1], max_iter = 5000)
        model.fit(X_train, y_train)
        
        y_pred = model.predict(X_test)
        r2_scores.append(r2_score(y_test, y_pred))
        coef_list.append(model.coef_)
        fold_df = pd.DataFrame({
            'fold': fold,
            'y_true': y_test.values,
            'y_pred': y_pred
        }, index=y_test.index)
        predictions.append(fold_df)
    all_predictions = pd.concat(predictions)
    return r2_scores, coef_list, all_predictions

def nested_cv_classification(X, y, outer_cv_folds=3):
    """
    Nested CV for binary classification using LogisticRegressionCV.
    Outer loop splits data; inner tuning is handled by LogisticRegressionCV.
    Returns outer F1 scores, coefficient vectors, and predictions (with probabilities) for each fold.
    """
    outer_cv = StratifiedKFold(n_splits=outer_cv_folds, shuffle=True, random_state=0)
    f1_scores, coef_list, predictions = [], [], []
    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train_raw)
        X_test = scaler.transform(X_test_raw)
        
        #model = LogisticRegressionCV(cv=3, penalty='l2', solver='liblinear', random_state=123)
        model = LogisticRegressionCV(cv=3, penalty='elasticnet', solver='saga', random_state=0, l1_ratios = [.1, .5, .7, .9, .95, .99, 1], max_iter = 5000)
        model.fit(X_train, y_train)
        
        preds = model.predict(X_test)
        probas = model.predict_proba(X_test)[:, 1]  # probability for positive class
        f1_scores.append(f1_score(y_test, preds))
        coef_list.append(model.coef_.flatten())
        fold_df = pd.DataFrame({
            'fold': fold,
            'y_true': y_test.values,
            'predicted_class': preds,
            'predicted_proba': probas
        }, index=y_test.index)
        predictions.append(fold_df)
    all_predictions = pd.concat(predictions)
    return f1_scores, coef_list, all_predictions

def nested_cv_pairwise_classification(X, y, feature_names, outer_cv_folds=3):
    """
    Nested CV for pairwise binary classification using LogisticRegressionCV.
    Returns outer F1 scores, coefficient vectors, and predictions for each fold.
    """
    outer_cv = StratifiedKFold(n_splits=outer_cv_folds, shuffle=True, random_state=0)
    f1_scores, coef_list, predictions = [], [], []
    for fold, (train_idx, test_idx) in enumerate(outer_cv.split(X, y)):
        X_train_raw, X_test_raw = X[train_idx], X[test_idx]
        y_train, y_test = y.iloc[train_idx], y.iloc[test_idx]
        scaler = StandardScaler()
        X_train = scaler.fit_transform(X_train_raw)
        X_test = scaler.transform(X_test_raw)
        
        #model = LogisticRegressionCV(cv=3, penalty='l2', solver='liblinear', random_state=123)
        model = LogisticRegressionCV(cv=3, penalty='elasticnet', solver='saga', random_state=0, l1_ratios = [.1, .5, .7, .9, .95, .99, 1], max_iter = 5000)
        model.fit(X_train, y_train)
        
        preds = model.predict(X_test)
        probas = model.predict_proba(X_test)[:, 1]
        f1_scores.append(f1_score(y_test, preds, average='macro'))
        coef_list.append(model.coef_.flatten())
        fold_df = pd.DataFrame({
            'fold': fold,
            'y_true': y_test.values,
            'predicted_class': preds,
            'predicted_proba': probas
        }, index=y_test.index)
        predictions.append(fold_df)
    all_predictions = pd.concat(predictions)
    return f1_scores, coef_list, all_predictions


def main_method(exp_path):
    
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756' # FA
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539' # Caffeine
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618' # Spermine
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407' # Lactic acid
    # exp_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655' # LiCl
    
    #exp_path = exp_folder  # Ensure exp_folder is defined in your environment
    ms_filename = 'ms_output_knnimputed.tsv'
    
    # --- Regression Analysis on treated samples ---
    merged_df, mets = merge_data(exp_path, ms_filename)
    treated_df = merged_df[merged_df["Treatment"].str.lower() == "yes".lower()]
    X_reg = treated_df[mets].values
    y_reg = treated_df["AUC"]
    
    r2_scores, reg_coef_list, reg_predictions = nested_cv_regression(X_reg, y_reg, outer_cv_folds=5)
    print("Regression Nested CV (R²): {:.2f} ± {:.2f}".format(np.mean(r2_scores), np.std(r2_scores)))
    
    avg_reg_coef = average_coefficients(reg_coef_list, mets)
    reg_coef_out = os.path.join(exp_path, "results", "metabolomics", "coefficients", "elasticnet_averaged_coefficients.csv")
    avg_reg_coef.to_csv(reg_coef_out)
    reg_pred_out = os.path.join(exp_path, "results", "metabolomics", "predictions", "elasticnet_regression_predictions.csv")
    reg_predictions.to_csv(reg_pred_out)
    #print(f"Saved regression coefficients to '{reg_coef_out}' and predictions to '{reg_pred_out}'.")
    
    # --- Classification Analysis (Treated vs. Untreated) ---
    merged_df["Treatment_binary"] = merged_df["Treatment"].str.lower().map({"yes": 1, "none": 0})
    X_clf = merged_df[mets].values
    y_clf = merged_df["Treatment_binary"]
    
    f1_scores_clf, clf_coef_list, clf_predictions = nested_cv_classification(X_clf, y_clf, outer_cv_folds=5)
    print("Logistic Regression Nested CV (F1): {:.2f} ± {:.2f}".format(np.mean(f1_scores_clf), np.std(f1_scores_clf)))
    
    avg_clf_coef = average_coefficients(clf_coef_list, mets)
    clf_coef_out = os.path.join(exp_path, "results", "metabolomics", "coefficients", "logreg_averaged_coefficients.csv")
    avg_clf_coef.to_csv(clf_coef_out)
    clf_pred_out = os.path.join(exp_path, "results", "metabolomics", "predictions", "logreg_classification_predictions.csv")
    clf_predictions.to_csv(clf_pred_out)
    #print(f"Saved classification coefficients to '{clf_coef_out}' and predictions to '{clf_pred_out}'.")
    
    # --- Pairwise Binary Classification on Supplement Groups ---
    merged_pair, mets_pair = load_metabolomics_data_with_design(exp_path, ms_filename)
    nontreated_merged_pair = merged_pair[merged_pair['Treatment'] == 'None']
    
    # Pairwise: "None" vs "PosLow"
    subset_low = nontreated_merged_pair[nontreated_merged_pair["Supplement"].isin(["None", "PosLow"])]
    X_low = subset_low[mets_pair].values
    y_low = subset_low["Supplement"]
    f1_low, coef_low, pair_pred_low = nested_cv_pairwise_classification(X_low, y_low, mets_pair, outer_cv_folds=3)
    print("Pairwise Logistic Regression (None vs PosLow) F1: {:.2f} ± {:.2f}".format(np.mean(f1_low), np.std(f1_low)))
    
    avg_coef_low = average_coefficients(coef_low, mets_pair)
    pair_low_pred_out = os.path.join(exp_path, "results", "metabolomics", "predictions", "logreg_none_vs_poslow_predictions.csv")
    pair_pred_low.to_csv(pair_low_pred_out)
    
    # Pairwise: "None" vs "PosHigh"
    subset_high = nontreated_merged_pair[nontreated_merged_pair["Supplement"].isin(["None", "PosHigh"])]
    X_high = subset_high[mets_pair].values
    y_high = subset_high["Supplement"]
    
    f1_high, coef_high, pair_pred_high = nested_cv_pairwise_classification(X_high, y_high, mets_pair, outer_cv_folds=3)
    print("Pairwise Logistic Regression (None vs PosHigh) F1: {:.2f} ± {:.2f}".format(np.mean(f1_high), np.std(f1_high)))
    
    avg_coef_high = average_coefficients(coef_high, mets_pair)
    pair_high_pred_out = os.path.join(exp_path, "results", "metabolomics", "predictions", "logreg_none_vs_poshigh_predictions.csv")
    pair_pred_high.to_csv(pair_high_pred_out)
    
    # Pairwise: "None" vs "NegHigh"
    subset_neghigh = nontreated_merged_pair[nontreated_merged_pair["Supplement"].isin(["None", "NegHigh"])]
    X_neghigh = subset_neghigh[mets_pair].values
    y_neghigh = subset_neghigh["Supplement"]
    
    f1_neghigh, coef_neghigh, pair_pred_neghigh = nested_cv_pairwise_classification(X_neghigh, y_neghigh, mets_pair, outer_cv_folds=3)
    print("Pairwise Logistic Regression (None vs Negigh) F1: {:.2f} ± {:.2f}".format(np.mean(f1_neghigh), np.std(f1_neghigh)))
    
    avg_coef_neghigh = average_coefficients(coef_neghigh, mets_pair)
    pair_neghigh_pred_out = os.path.join(exp_path, "results", "metabolomics", "predictions", "logreg_none_vs_neghigh_predictions.csv")
    pair_pred_neghigh['predicted_proba'] = 1 - pair_pred_neghigh['predicted_proba'] # Note that due to alphabetical reasons, NegHigh is the zero class.
    pair_pred_neghigh.to_csv(pair_neghigh_pred_out)
    
    # --- Pairwise Statistical Comparisons ---
    comp_df = pairwise_comparison_df(merged_pair, group_col="Supplement", metabolite_columns=mets_pair)
    comp_out = os.path.join(exp_path, "results", "metabolomics", "coefficients", "pairwise_comparisons.csv")
    comp_df.to_csv(comp_out, index=False)
    #print(f"Saved pairwise comparisons to '{comp_out}'.")
    
    # Save the annotated dataframe   
    merged_df.to_csv(os.path.join(exp_path,'results/metabolomics/processed/data_with_annotation.tsv'), sep = '\t')
    


    '''
    # Create the parity plot for reg trhingy
    
    import matplotlib.pyplot as plt
    
    plt.figure(figsize=(4, 4))
    plt.scatter(reg_predictions['y_true'], reg_predictions['y_pred'], color='blue', label='Data points')
    
    # Determine limits for the parity line
    min_val = min(reg_predictions['y_true'].min(), reg_predictions['y_pred'].min())
    max_val = max(reg_predictions['y_true'].max(), reg_predictions['y_pred'].max())
    plt.plot([min_val, max_val], [min_val, max_val], 'r--', label='Parity line')
    
    # Set labels and title
    plt.xlabel('True')
    plt.ylabel('Predicted')
    plt.title('Parity Plot')
    #plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()
    
    # Compute absolute values and select top 10 features by absolute mean coefficient
    df_temporary = avg_reg_coef
    df_temporary['Abs_Mean'] = df_temporary['Mean_Coefficient'].abs()
    top_df = df_temporary.sort_values(by='Abs_Mean', ascending=False).head(10)
    
    # Create barplot with error bars
    plt.figure(figsize=(6, 5))
    plt.bar(top_df.index, top_df['Mean_Coefficient'], yerr=top_df['Std_Coefficient'], capsize=5)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Coefficient")
    plt.title("Top 10 Metabolites for prediction of formic acid resistance")
    plt.tight_layout()
    plt.show()
    
    
    
    
    '''

if __name__ == "__main__":
    
    import argparse
    parser = argparse.ArgumentParser(description="Multivariate and univariate analysis for metabolomics")
    parser.add_argument("--output_folder", required=True, type=str, help="Experiment folder")

    args = parser.parse_args()
    main_method(args.output_folder)
    
    main_method()



import os
import argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV, LogisticRegressionCV
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.metrics import r2_score, roc_auc_score, matthews_corrcoef
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning

# Suppress convergence warnings
simplefilter("ignore", category=ConvergenceWarning)


def average_coefficients(coef_list, feature_names):
    """Average coefficient vectors and return a DataFrame of means and std deviations."""
    coefs = np.vstack(coef_list)
    df = pd.DataFrame(coefs, columns=feature_names)
    summary = df.agg(['mean', 'std']).T.rename(columns={'mean': 'Mean_Coefficient', 'std': 'Std_Coefficient'})
    return summary


def pairwise_comparison_df(df, group_col, features, alpha=0.05, method='fdr_bh'):
    """Perform pairwise Mann–Whitney U tests with FDR correction for given groups in df."""
    records = []
    groups = sorted(df[group_col].unique())
    for g1 in groups:
        for g2 in groups:
            if g1 >= g2:
                continue
            for feat in features:
                d1 = df[df[group_col] == g1][feat]
                d2 = df[df[group_col] == g2][feat]
                stat, p = mannwhitneyu(d1, d2, alternative='two-sided', nan_policy='omit')
                records.append((g1, g2, feat, stat, p))
    res = pd.DataFrame(records, columns=['Group1','Group2','Feature','U_stat','p_raw'])
    reject, p_corr, _, _ = multipletests(res['p_raw'], alpha=alpha, method=method)
    res['p_corrected'] = p_corr
    res['significant'] = reject
    return res


def merge_data(exp_path, ms_filename, include_growth=False):
    """Load metabolomics (and optionally growth) data, merge with design, impute features."""
    met_path = os.path.join(exp_path, 'results', 'metabolomics', 'processed', ms_filename)
    design_path = os.path.join(exp_path, 'protocol', 'plate_layout', 'design_table.tsv')

    met_df = pd.read_csv(met_path, sep='\t', index_col=0)
    met_df['Well'] = met_df.index.str.split('-').str[2].str.split('.').str[0]
    design_df = pd.read_csv(design_path, sep='\t')
    df = met_df.merge(design_df, on='Well', how='left')

    if include_growth:
        growth_path = os.path.join(exp_path, 'results', 'growth', 'processed', 'growth_parameters.tsv')
        growth_df = pd.read_csv(growth_path, sep='\t', index_col=0).reset_index().rename(columns={'index':'Well'})
        df = df.merge(growth_df, on='Well', how='inner')

    features = [c for c in met_df.columns if c not in ['Experimental group','Well']]
    df[features] = df[features].fillna(df[features].median())
    return df, features


def nested_cv_regression(X, y, n_splits=5):
    """Nested CV regression with ElasticNetCV. Returns R2 scores, coefs, and preds."""
    outer = KFold(n_splits=n_splits, shuffle=True, random_state=0)
    scores, coefs, preds = [], [], []
    for fold, (tr, te) in enumerate(outer.split(X)):
        X_tr, X_te = X[tr], X[te]
        y_tr, y_te = y.iloc[tr], y.iloc[te]
        scaler = StandardScaler().fit(X_tr)
        X_tr_s, X_te_s = scaler.transform(X_tr), scaler.transform(X_te)
        model = ElasticNetCV(cv=5, random_state=0,
                             l1_ratio=[.1, .5, .7, .9, 1], max_iter=5000)
        model.fit(X_tr_s, y_tr)
        y_pred = model.predict(X_te_s)
        scores.append(r2_score(y_te, y_pred))
        coefs.append(model.coef_)
        preds.append(pd.DataFrame({'Well': y_te.index, 'y_true': y_te.values, 'y_pred': y_pred}))
    print(f"Treatment resistance R²: {np.mean(scores):.2f} ± {np.std(scores):.2f}")
    return scores, coefs, pd.concat(preds, ignore_index=True)


def leave_one_out_classification_metrics(X, y, threshold=0.5):
    """LOO CV classification: returns ROC-AUC, MCC, coefs, and preds."""
    loo = LeaveOneOut()
    y_true_all, probas_all, coefs = [], [], []
    preds = []
    for idx_train, idx_test in loo.split(X):
        X_tr, X_te = X[idx_train], X[idx_test]
        y_tr, y_te = y.iloc[idx_train], y.iloc[idx_test]
        scaler = StandardScaler().fit(X_tr)
        X_tr_s, X_te_s = scaler.transform(X_tr), scaler.transform(X_te)
        model = LogisticRegressionCV(cv=3, penalty='elasticnet', solver='saga', random_state=0,
                                     l1_ratios=[.1, .5, .7, .9, 1], max_iter=5000)
        model.fit(X_tr_s, y_tr)
        probas = model.predict_proba(X_te_s)[0,1]
        y_true_all.append(y_te.values[0])
        probas_all.append(probas)
        coefs.append(model.coef_.flatten())
        preds.append(pd.DataFrame({'Well': y_te.index, 'y_true': y_te.values, 'predicted_proba': probas}))

    # Compute metrics
    auc = roc_auc_score(y_true_all, probas_all)
    y_pred_all = [1 if p >= threshold else 0 for p in probas_all]
    mcc = matthews_corrcoef(y_true_all, y_pred_all)
    print(f"Classification ROC-AUC (LOO): {auc:.2f}")
    print(f"Classification MCC (threshold={threshold}) (LOO): {mcc:.2f}")
    return auc, mcc, coefs, pd.concat(preds, ignore_index=True)


def main(exp_path, ms_filename='ms_output_knnimputed.tsv'):
    
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756' # FA
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539' # Caffeine
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618' # Spermine
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407' # Lactic acid
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655' # LiCl
    exp_path = exp_folder
    
    df, features = merge_data(exp_path, ms_filename, include_growth=True)
    score_dict = {}
    pred_list = []

    # Regression
    treat_df = df[df['Treatment'].str.lower() == 'yes']
    r2_scores, reg_coefs, reg_preds = nested_cv_regression(treat_df[features].values, treat_df['AUC'])
    reg_preds.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'predictions', 'elasticnet_regression_preds.csv'), index=False)

    # Classification treated vs untreated (LOO)
    df['Treatment_binary'] = df['Treatment'].str.lower().map({'yes':1,'none':0})
    key = 'treatment_vs_none'
    print(key)
    auc, mcc, clf_coefs, clf_preds = leave_one_out_classification_metrics(df[features].values, df['Treatment_binary'])
    score_dict[f"{key}_auc"] = auc
    score_dict[f"{key}_mcc"] = mcc
    pred = clf_preds.copy()
    pred['classification'] = key
    pred_list.append(pred)

    # Pairwise supplement classifications
    pairs = [
        ('None', 'None', 'PosLow'),
        ('None', 'None', 'PosHigh'),
        ('None', 'None', 'NegHigh'),
        ('Yes', 'None', 'PosLow'),
        ('Yes', 'None', 'PosHigh'),
        ('Yes', 'None', 'NegHigh'),
    ]
    
    pairs = [
        ('Yes', 'None', 'PosLow'),
    ]
    for treatment, g1, g2 in pairs:
        subset = df[(df['Treatment'] == treatment) & df['Supplement'].isin([g1, g2])]
        y = subset['Supplement'].map({g1:0, g2:1})
        key = f"{treatment}_{g1}_vs_{g2}"
        print(f"Pairwise classification: {key}")
        auc, mcc, coefs, preds = leave_one_out_classification_metrics(subset[features].values, y)
        score_dict[f"{key}_auc"] = auc
        score_dict[f"{key}_mcc"] = mcc
        pred = preds.copy()
        pred['classification'] = key
        pred_list.append(pred)

    # Save scores (one row with metrics)
    scores_df = pd.DataFrame([score_dict])
    scores_df.index.name = 'LOO'
    scores_df.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'predictions', 'classification_scores.tsv'), sep='\t')

    # Save predictions (for ROC/MCC analysis)
    all_preds = pd.concat(pred_list, ignore_index=True)
    all_preds.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'predictions', 'all_classification_predictions.tsv'), sep='\t', index=False)

    # Univariate stats
    comp_df = pairwise_comparison_df(df, 'Supplement', features)
    comp_df.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'coefficients', 'pairwise_comparisons.csv'), index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Metabolomics analysis')
    parser.add_argument('--exp_path', required=True, help='Experiment folder')
    args = parser.parse_args()
    main(args.exp_path)










import os
import argparse
import numpy as np
import pandas as pd
from sklearn.linear_model import ElasticNetCV
from sklearn.svm import SVC
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import KFold, LeaveOneOut
from sklearn.metrics import r2_score, roc_auc_score, matthews_corrcoef
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests
from warnings import simplefilter
from sklearn.exceptions import ConvergenceWarning

# Suppress convergence warnings
simplefilter("ignore", category=ConvergenceWarning)


def average_coefficients(coef_list, feature_names):
    """Average coefficient vectors and return a DataFrame of means and std deviations."""
    coefs = np.vstack(coef_list)
    df = pd.DataFrame(coefs, columns=feature_names)
    summary = df.agg(['mean', 'std']).T.rename(columns={'mean': 'Mean_Coefficient', 'std': 'Std_Coefficient'})
    return summary


def pairwise_comparison_df(df, group_col, features, alpha=0.05, method='fdr_bh'):
    """Perform pairwise Mann–Whitney U tests with FDR correction for given groups in df."""
    records = []
    groups = sorted(df[group_col].unique())
    for g1 in groups:
        for g2 in groups:
            if g1 >= g2:
                continue
            for feat in features:
                d1 = df[df[group_col] == g1][feat]
                d2 = df[df[group_col] == g2][feat]
                stat, p = mannwhitneyu(d1, d2, alternative='two-sided', nan_policy='omit')
                records.append((g1, g2, feat, stat, p))
    res = pd.DataFrame(records, columns=['Group1','Group2','Feature','U_stat','p_raw'])
    reject, p_corr, _, _ = multipletests(res['p_raw'], alpha=alpha, method=method)
    res['p_corrected'] = p_corr
    res['significant'] = reject
    return res


def merge_data(exp_path, ms_filename, include_growth=False):
    """Load metabolomics (and optionally growth) data, merge with design, impute features."""
    met_path = os.path.join(exp_path, 'results', 'metabolomics', 'processed', ms_filename)
    design_path = os.path.join(exp_path, 'protocol', 'plate_layout', 'design_table.tsv')

    met_df = pd.read_csv(met_path, sep='\t', index_col=0)
    met_df['Well'] = met_df.index.str.split('-').str[2].str.split('.').str[0]
    design_df = pd.read_csv(design_path, sep='\t')
    df = met_df.merge(design_df, on='Well', how='left')

    if include_growth:
        growth_path = os.path.join(exp_path, 'results', 'growth', 'processed', 'growth_parameters.tsv')
        growth_df = pd.read_csv(growth_path, sep='\t', index_col=0).reset_index().rename(columns={'index':'Well'})
        df = df.merge(growth_df, on='Well', how='inner')

    features = [c for c in met_df.columns if c not in ['Experimental group','Well']]
    df[features] = df[features].fillna(df[features].median())
    return df, features


def nested_cv_regression(X, y, n_splits=5):
    """Nested CV regression with ElasticNetCV. Returns R2 scores, coefs, and preds."""
    outer = KFold(n_splits=n_splits, shuffle=True, random_state=0)
    scores, coefs, preds = [], [], []
    for fold, (tr, te) in enumerate(outer.split(X)):
        X_tr, X_te = X[tr], X[te]
        y_tr, y_te = y.iloc[tr], y.iloc[te]
        scaler = StandardScaler().fit(X_tr)
        X_tr_s, X_te_s = scaler.transform(X_tr), scaler.transform(X_te)
        
        # model
        model = ElasticNetCV(cv=5, random_state=0,
                             l1_ratio=[.1, .5, .7, .9, 1], max_iter=5000)
        model.fit(X_tr_s, y_tr)
        
        # interpretation
        y_pred = model.predict(X_te_s)
        scores.append(r2_score(y_te, y_pred))
        coefs.append(model.coef_)
        preds.append(pd.DataFrame({'Well': y_te.index, 'y_true': y_te.values, 'y_pred': y_pred}))
        
    print(f"Treatment resistance R²: {np.mean(scores):.2f} ± {np.std(scores):.2f}")
    return scores, coefs, pd.concat(preds, ignore_index=True)


def leave_one_out_svm(X, y):
    """LOO CV classification with linear SVM; returns ROC-AUC, MCC, coefs, and preds."""
    loo = LeaveOneOut()
    y_true_all, decision_scores, y_pred_all, coefs = [], [], [], []
    preds = []
    for train_idx, test_idx in loo.split(X):
        X_tr, X_te = X[train_idx], X[test_idx]
        y_tr, y_te = y.iloc[train_idx], y.iloc[test_idx]
        scaler = StandardScaler().fit(X_tr)
        X_tr_s, X_te_s = scaler.transform(X_tr), scaler.transform(X_te)
        
        # model
        model = SVC(kernel='linear', C=1.0, class_weight='balanced')
        model.fit(X_tr_s, y_tr)
        
        # interpretation
        score = model.decision_function(X_te_s)[0]
        pred_label = int(score >= 0)
        y_true_all.append(y_te.values[0])
        decision_scores.append(score)
        y_pred_all.append(pred_label)
        coefs.append(model.coef_.flatten())
        preds.append(pd.DataFrame({'Well': y_te.index,
                                   'y_true': y_te.values,
                                   'decision_score': score,
                                   'predicted_class': pred_label}))
    # Compute metrics
    auc = roc_auc_score(y_true_all, decision_scores)
    #mcc = matthews_corrcoef(y_true_all, y_pred_all)
    print(f"Linear SVM ROC-AUC (LOO): {auc:.2f}")
    #print(f"Linear SVM MCC (LOO): {mcc:.2f}")
    return auc, coefs, pd.concat(preds, ignore_index=True)


def main(exp_path, ms_filename='ms_output_knnimputed.tsv'):
    
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202503141756' # FA
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503131539' # Caffeine
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/glutamate_202501281618' # Spermine
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/proline_202503051407' # Lactic acid
    # exp_folder = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/completed_experiments/arginine_202503141655' # LiCl
    
    exp_path = exp_folder
    
    df, features = merge_data(exp_path, ms_filename, include_growth=True)
    score_dict = {}
    pred_list = []

    # Regression (ElasticNet)
    treat_df = df[df['Treatment'].str.lower() == 'yes']
    _, _, reg_preds = nested_cv_regression(treat_df[features].values, treat_df['AUC'])
    reg_preds.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'predictions', 'elasticnet_regression_preds.csv'), index=False)

    # Linear SVM Classification treated vs untreated (LOO)
    df['Treatment_binary'] = df['Treatment'].str.lower().map({'yes':1,'none':0})
    key = 'treatment_vs_none_svm'
    auc, svm_coefs, preds = leave_one_out_svm(df[features].values, df['Treatment_binary'])
    score_dict[f"{key}"] = auc
    avg_coefs = average_coefficients(svm_coefs, features)
    avg_coefs.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'coefficients', f'svm_{key}_avg_coefs.csv'))
    preds['classification'] = key
    pred_list.append(preds)
    preds.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'predictions', f'svm_{key}_preds.csv'), index=False)

    # Pairwise supplement classifications with SVM
    pairs = [
        ('None', 'None', 'PosLow'),
        ('None', 'None', 'PosHigh'),
        ('None', 'None', 'NegHigh'),
        ('Yes',  'None', 'PosLow'),
        ('Yes',  'None', 'PosHigh'),
        ('Yes',  'None', 'NegHigh'),
    ]
    for treatment, g1, g2 in pairs:
        subset = df[(df['Treatment'] == treatment) & df['Supplement'].isin([g1, g2])]
        y = subset['Supplement'].map({g1:0, g2:1})
        key = f"{treatment}_{g1}_vs_{g2}_svm"
        print(f"Pairwise SVM: {key}")
        auc, coefs, preds = leave_one_out_svm(subset[features].values, y)
        score_dict[f"{key}"] = auc
        #score_dict[f"{key}_mcc"] = mcc
        avg_coefs = average_coefficients(coefs, features)
        avg_coefs.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'coefficients', f'svm_{key}_avg_coefs.csv'))
        preds['classification'] = key
        pred_list.append(preds)
        preds.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'predictions', f'svm_{key}_preds.csv'), index=False)

    # Save combined scores
    scores_df = pd.DataFrame([score_dict]).transpose()
    scores_df.columns = ['AUC']
    scores_df.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'predictions', 'classification_scores.tsv'), sep='\t')

    # Save combined predictions
    all_preds = pd.concat(pred_list, ignore_index=True)
    all_preds.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'predictions', 'all_classification_predictions.tsv'), sep='\t', index=False)

    # Univariate stats
    comp_df = pairwise_comparison_df(df, 'Supplement', features)
    comp_df.to_csv(os.path.join(exp_path, 'results', 'metabolomics', 'coefficients', 'pairwise_comparisons_svm.csv'), index=False)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Metabolomics analysis with SVM')
    parser.add_argument('--exp_path', required=True, help='Experiment folder')
    args = parser.parse_args()
    main(args.exp_path)








