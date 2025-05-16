

import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
import argparse

import config

from utils.processing_utils import *
from growth_stat_testing import growth_testing#, two_way_anova_testing



#EXPERIMENT_DIR = Path('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/histidine_202504270959')
#EXPERIMENT_DIR = Path('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/lysine_202505071046')


#@task
#def save_plots(EXPERIMENT_DIR, growth_curves_filtered_smoothed,unfiltered_growth_curves_smoothed, violinplot_df, mu_per_well):
def save_plots(EXPERIMENT_DIR, growth_curves_filtered_smoothed, raw_growth_curves, violinplot_df, mu_per_well, layout):
    
    raw_growth_curves.index = raw_growth_curves.index.map(lambda x: unify_well_format(str(x)))
    layout = layout.copy()
    layout['well'] = layout['well'].astype(str).apply(unify_well_format)
    
    # Merge layout information (e.g. Summary) into the data.
    raw_growth_curves = raw_growth_curves.merge(layout[['well', 'Summary']], left_index=True, right_on='well').set_index('well')

    # Use a non-interactive backend
    plt.switch_backend('Agg')  # Ensures figures are not displayed, only saved

    plot_group_averages_in_hours(growth_curves_filtered_smoothed, group_col = "Summary")
    plt.savefig(EXPERIMENT_DIR / 'results/plots/averaged_growth_curves.png')   # save the figure to file
    plt.close('all')    # close the figure window

    # Plot a grid of them? along with the windows extracted for the growth rate calculation
    plot_growth_curves(raw_growth_curves, annotations_df=mu_per_well)
    plt.savefig(EXPERIMENT_DIR / 'results/plots/growth_curves.png')   # save the figure to file
    plt.close('all')    # close the figure window

    # Prep boxplot/violinplots?
    
    plot_boxplots(violinplot_df.drop('MaxOD', axis = 1))
    #plot_violinplots(violinplot_df.drop('MaxOD', axis = 1))
    plt.savefig(EXPERIMENT_DIR / 'results/plots/growth_metrics.png')   # save the figure to file
    plt.close('all')    # close the figure window

#@task
def save_growth_properties(EXPERIMENT_DIR, mu_per_well, auc_per_well,
                           finalOD_per_well, growth_curves_filtered_smoothed):
    
    # Save file with properties
    violinplot_df = pd.merge(mu_per_well['mu'], auc_per_well, 
                             left_index = True, right_index = True).merge(finalOD_per_well, 
                                                      left_index = True, right_index = True).merge(growth_curves_filtered_smoothed['Summary'], 
                                                                          left_index = True, right_index = True)
    violinplot_df.to_csv(EXPERIMENT_DIR / 'results/growth/processed/growth_parameters.tsv', sep = '\t')
    
    return violinplot_df

#@task
def extract_growth_properties(growth_curves_filtered_smoothed, layout):
    
    auc_per_well = compute_auc(growth_curves_filtered_smoothed.iloc[:,:-1])
    mu_per_well = extract_growth_rates(layout, growth_curves_filtered_smoothed)
    finalOD_per_well = extract_finalOD(layout, growth_curves_filtered_smoothed)
    
    return auc_per_well, mu_per_well, finalOD_per_well

#@task
#def save_interim_data(EXPERIMENT_DIR, growth_curves, 
#                      growth_curves_filtered, 
#                      unfiltered_growth_curves_smoothed, 
#                      growth_curves_filtered_smoothed):
#    # Save the processed curves and unprocessed curves?
#    growth_curves.to_csv(EXPERIMENT_DIR / 'results/growth/processed/growth_curves.tsv', sep = '\t')
#    growth_curves_filtered.to_csv(EXPERIMENT_DIR / 'results/growth/processed/curated_growth_curves.tsv', sep = '\t')
#    unfiltered_growth_curves_smoothed.to_csv(EXPERIMENT_DIR / 'results/growth/processed/smoothed_growth_curves.tsv', sep = '\t')
#    growth_curves_filtered_smoothed.to_csv(EXPERIMENT_DIR / 'results/growth/processed/curated_smoothed_growth_curves.tsv', sep = '\t')
    
    
    
def save_interim_data(EXPERIMENT_DIR, raw_growth_curves, 
                      curated_growth_curves, 
                      blanked_growth_curves, 
                      smoothed_growth_curves):
    # Save the processed curves and unprocessed curves?
    raw_growth_curves.to_csv(EXPERIMENT_DIR / 'results/growth/processed/raw_growth_curves.tsv', sep = '\t')
    curated_growth_curves.to_csv(EXPERIMENT_DIR / 'results/growth/processed/curated_growth_curves.tsv', sep = '\t')
    blanked_growth_curves.to_csv(EXPERIMENT_DIR / 'results/growth/processed/blanked_growth_curves.tsv', sep = '\t')
    smoothed_growth_curves.to_csv(EXPERIMENT_DIR / 'results/growth/processed/curated_smoothed_growth_curves.tsv', sep = '\t')

#@task
def smooth(growth_curves_filtered, growth_curves, window_size):
    growth_curves_filtered_smoothed = smooth_growth_curves(growth_curves_filtered, window=window_size, method="mirror")
    #unfiltered_growth_curves_smoothed = smooth_growth_curves(growth_curves.drop('growth_summary', axis = 1), window=window_size, method="mirror")
    unfiltered_growth_curves_smoothed = smooth_growth_curves(growth_curves, window=window_size, method="mirror")
    return growth_curves_filtered_smoothed, unfiltered_growth_curves_smoothed

#@task
def filter_curves(layout, growth_curves, group_col, method, threshold=3):
    # Filter for outliers curves using MAD (default),
    growth_curves_filtered = filter_outlier_growth_curves(layout, growth_curves, group_col="Summary", method=method, threshold=threshold)
    return growth_curves_filtered

#@task   
def read_growth_data(EXPERIMENT_DIR):
    # Read in the data and metadata
    layout = pd.read_excel(EXPERIMENT_DIR / 'protocol/hamilton/pipetting_layout.xlsx')
    raw_growth_curves = process_measurement_data(EXPERIMENT_DIR / 'results/growth/raw')
    return layout, raw_growth_curves

#@task   
def blank_processing(layout, raw_growth_curves, n, fillin_value, blank_bool):
    growth_curves = subtract_and_impute_blanks(layout, raw_growth_curves, n, fillin_value, blank_bool)
    return growth_curves

#@flow
def run_growth_processing(EXPERIMENT_DIR, testing_bool):
    
    # Load config parameters
    EXPERIMENT_DIR = Path(EXPERIMENT_DIR)
    n = config.N_CLOSEST_BLANKS
    fillin_value = config.BLANK_FILLIN_VALUE
    window_size = config.ROLLING_MEAN_WINDOWSIZE
    threshold = config.THRESHOLD
    outlier_method = config.OUTLIER_METHOD
    
    # Read in the growth curves from the raw polarstar output
    layout, raw_growth_curves = read_growth_data(EXPERIMENT_DIR)
    
    # Filter low-quality growth curves using MAD df = raw_growth_curves
    filtered_growth_curves = filter_curves(layout, raw_growth_curves, 'summary', outlier_method, threshold)
    
    #growth_curves = blank_processing(layout,raw_growth_curves, n, fillin_value, config.BLANK_SUBTRACTION)
    filtered_blanked_growth_curves = blank_processing(layout,filtered_growth_curves, n, fillin_value, config.BLANK_SUBTRACTION)
    
    # Smooth curves using a rolling median
    #growth_curves_filtered_smoothed, unfiltered_growth_curves_smoothed = smooth(growth_curves_filtered, growth_curves, window_size)
    growth_curves_filtered_smoothed, unfiltered_growth_curves_smoothed = smooth(filtered_blanked_growth_curves, filtered_growth_curves, window_size)
    
    # Save growth curves in different stages of processing
    save_interim_data(EXPERIMENT_DIR, raw_growth_curves, filtered_growth_curves, filtered_blanked_growth_curves, growth_curves_filtered_smoothed)
    
    # Extract growth parameters
    auc_per_well, mu_per_well, finalOD_per_well = extract_growth_properties(growth_curves_filtered_smoothed, layout)
    
    # Save growth parameters
    violinplot_df = save_growth_properties(EXPERIMENT_DIR, mu_per_well, auc_per_well, finalOD_per_well, growth_curves_filtered_smoothed)
    
    # Save partial report
    save_plots(EXPERIMENT_DIR, growth_curves_filtered_smoothed, raw_growth_curves, violinplot_df, mu_per_well, layout)
    
    # if we want to run some basic sign testing
    if testing_bool == 'yes':
        growth_testing(EXPERIMENT_DIR)
        #two_way_anova_testing(EXPERIMENT_DIR)
        
    
if __name__ == "__main__":
    
    parser = argparse.ArgumentParser(description="Run growth processing pipeline.")
    parser.add_argument('--output_folder', required=True,
                        help="folder")
    parser.add_argument('--testing', required=True, help="do testing?")
    args = parser.parse_args()
    
    run_growth_processing(args.output_folder, args.testing)


    # Example:
    # python growth_processing.py --output_folder "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_202501281618" --testing yes


    # EXPERIMENT_DIR = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_202504291411'
    '''
     df = raw_growth_curves.copy()
     df.index = df.index.map(lambda x: unify_well_format(str(x)))
     layout_rev = layout.copy()
     layout_rev['well'] = layout_rev['well'].astype(str).apply(unify_well_format)
     
     # Merge layout information (e.g. Summary) into the data.
     df = df.merge(layout_rev[['well', 'Summary']], left_index=True, right_on='well').set_index('well')
    '''