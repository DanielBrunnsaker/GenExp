import pandas as pd
import numpy as np
import os
from sklearn.metrics import root_mean_squared_error, mean_absolute_error, r2_score
from sklearn.model_selection import KFold
from sklearn.linear_model import ElasticNet, ElasticNetCV
from skopt import BayesSearchCV
from skopt.space import Real
import json
from pathlib import Path

def return_scoring(y_true, y_pred):
    
    r2 = r2_score(y_true, y_pred)
    mae = mean_absolute_error(y_true, y_pred)
    rmse = root_mean_squared_error(y_true, y_pred)
    
    return r2, mae, rmse

def weight_patterns(home_folder):
    
    import warnings
    warnings.simplefilter("ignore")
    os.environ["PYTHONWARNINGS"] = "ignore" # Also affect subprocesses
    
    os.chdir(home_folder) # or set to whatever your scripts folder is
    
    aap_set = pd.read_excel('../data/AA.xls', 
                            sheet_name='intracellular_concentration_mM', index_col = 0).iloc[:,1:]
    aap_set = aap_set.reset_index().groupby('ORF').mean()
    kf = KFold(n_splits=10, shuffle=True, random_state=0)
    
    train_dict = {}
    test_dict = {}
    
    result_df_linear = pd.DataFrame(columns=['Amino acid', 'Fold', 'R2', 'MAE', 'RMSE'])
    
    # Loop through amino acids first!
    for kk, aa in enumerate(aap_set.columns[:19]):
    
        print(aa)
        # Save scores and predictions per amino acid
        full_frequent_dataset = pd.read_feather('../results/patterns/datasets/frequent_20241114.feather') #CHANGE THIS LINE
        full_frequent_dataset.set_index('ORF', inplace = True)
        coefficients = pd.DataFrame(index = full_frequent_dataset.columns)
    
        for fold_idx, (train_index, test_index) in enumerate(kf.split(list(aap_set.index))):
            
        
            fold = str(fold_idx+1)
            print('Fold:',fold)
            
            xtrain = aap_set.iloc[train_index,:].merge(full_frequent_dataset, left_index = True, right_index = True, how = 'left')
            xtest_full = aap_set.iloc[test_index,:].merge(full_frequent_dataset, left_index = True, right_index = True, how = 'left')
            
            ilp_indices = np.array([xtrain.columns.get_loc(col) for col in [col for col in xtrain.columns if col.startswith('Y')]]) # Define the pattern indices
            filtered_dx = xtrain.iloc[:,ilp_indices]#.loc[:, (xtrain.iloc[:,ilp_indices] == 1).sum(axis=0) >= 3] #remove with less than three positive examples
            
            y = xtrain[aa]
            y_test = xtest_full[aa]
            
            elasticnet_cv = BayesSearchCV(
                ElasticNet(), scoring="r2", search_spaces={"alpha": Real(5e-5, 1, 'log-uniform')}, 
                cv=10, n_jobs=-1, verbose = 0, random_state = 0, n_iter = 50 
            ).fit(filtered_dx, xtrain[aa])
            
            elasticnet = ElasticNet(alpha = elasticnet_cv.best_params_['alpha']).fit(filtered_dx, y)
            
            linear_predictions = elasticnet.predict(xtest_full[filtered_dx.columns])
            
            # Save scores
            r2_linear, mae_linear, rmse_linear = return_scoring(y_test, linear_predictions)
            
            new_row_df = pd.DataFrame([{
                'Amino acid': aa, 
                'Fold': fold, 
                'R2': r2_linear, 
                'MAE': mae_linear, 
                'RMSE': rmse_linear
            }])
            print('R2:',r2_linear)
            
            # Append the new row to the DataFrame
            result_df_linear = pd.concat([result_df_linear, new_row_df], ignore_index=True)
            result_df_linear.to_csv('../results/metrics/eCV_results_20241114.csv')
    
            # Define the coefficients
            coefficients = coefficients.merge(pd.DataFrame(elasticnet.coef_, index = filtered_dx.columns, columns = [str(fold)]), how = 'left', left_index = True, right_index = True)
            
            # Save the indices and coefficients for downstream analysis
            coefficients.to_csv('../results/coefficients/20241114/'+aa+'_eCV_coefficients.csv')
    
    avg_results = result_df_linear.groupby('Amino acid').mean()        

if __name__ == '__main__':
    
    BASE_DIR = Path(__file__).resolve().parent
    weight_patterns(BASE_DIR)

