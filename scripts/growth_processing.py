#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 15:45:22 2024

@author: danbru
"""

def generate_list(start_letter, end_letter, start_num, end_num):
    
    import string

    result = []
    
    for letter in string.ascii_uppercase[string.ascii_uppercase.index(start_letter):string.ascii_uppercase.index(end_letter) + 1]:
        for num in range(start_num, end_num + 1):
            result.append(f"{letter}{num:02}")
    
    return result

def well_to_position(well):
    row = well[0]
    col = int(well[1:])
    return row, col

def get_well_index(well):
    row, col = well_to_position(well)
    return ord(row) - ord('A'), col - 1

def well_distance(well1, well2):
    row1, col1 = get_well_index(well1)
    row2, col2 = get_well_index(well2)
    return np.sqrt((row1 - row2) ** 2 + (col1 - col2) ** 2)

def subtract_closest_blanks(df, blank_wells, n):
    result = df.copy()
    for well in df.index:
        
        # Calculate distances from the current well to each blank well. Maybe change how i calculate this distance? 
        # Maybe Manhattah makes more sense?
        distances = [(blank, well_distance(well, blank)) for blank in blank_wells]
        
        # Sort by distance and take the closest n blanks
        closest_blanks = sorted(distances, key=lambda x: x[1])[:n]
        
        # Get the values of the closest blank wells for each timepoint
        closest_blank_values = df.loc[[blank for blank, _ in closest_blanks]].mean()
        
        result.loc[well] -= closest_blank_values
    
    return result

def remove_timepoint_outliers_mad(df, blank_wells,threshold=1):
  # Calculate the median and MAD across all wells for each time point
  #df = df[~df.index.isin(blank_wells)]
  time_point_medians = df[~df.index.isin(blank_wells)].median(axis=0)
  mad = (df - time_point_medians).abs().median(axis=0)

  # Set a threshold based on MAD
  mad_threshold = mad.median() + threshold * mad.median()

  # Identify outlier time points where the deviation exceeds the MAD threshold
  outlier_timepoints = mad > mad_threshold

  # Replace outlier time points with NaN
  filtered_df = df.copy()
  filtered_df.loc[:, outlier_timepoints] = np.nan

  return filtered_df, outlier_timepoints


def remove_timepoint_outliers_iqr(data_matrix, iqr_threshold=1.5):
    # Calculate the 1st and 3rd quartiles (25th and 75th percentiles) at each time point
    Q1 = np.percentile(data_matrix, 25, axis=0)
    Q3 = np.percentile(data_matrix, 75, axis=0)
    IQR = Q3 - Q1

    # Define outlier thresholds for each time point
    lower_bound = Q1 - iqr_threshold * IQR
    upper_bound = Q3 + iqr_threshold * IQR

    # Flag time points where most or all wells are outside the expected range
    outlier_timepoints = np.any((data_matrix < lower_bound) | (data_matrix > upper_bound), axis=0)

    # Replace outlier time points with NaN across all wells
    filtered_data = data_matrix.copy()
    filtered_data[:, outlier_timepoints] = np.nan

    return filtered_data, outlier_timepoints

def remove_timepoint_outliers_derivative(df, max_change=0.5):
    # Calculate the difference between consecutive time points (rate of change)
    rate_of_change = df.diff(axis=1).abs()

    # Identify time points where most wells have a large rate of change
    outlier_timepoints = (rate_of_change > max_change).any(axis=0)

    # Replace the values in outlier time points with NaN
    filtered_df = df.copy()
    filtered_df.loc[:, outlier_timepoints] = np.nan

    return filtered_df, outlier_timepoints

def remove_outliers_per_well_mad(df, threshold=3):
    # Create a copy of the DataFrame to store the filtered data
    filtered_df = df.copy()
    
    # Iterate over each well (each row in the DataFrame)
    for well in df.index:
        # Get the data for the current well
        well_data = df.loc[well]
        
        # Calculate the median for the well
        median_value = np.median(well_data)
        
        # Calculate the absolute deviation from the median
        mad = np.median(np.abs(well_data - median_value))
        
        # Set a threshold based on MAD
        mad_threshold = threshold * mad
        
        # Replace values that deviate more than the threshold with NaN
        filtered_df.loc[well] = np.where(np.abs(well_data - median_value) > mad_threshold, np.nan, well_data)
    
    return filtered_df

def remove_outliers_based_on_growth_rate(df, max_growth_rate=0.5, max_drop_rate=-0.5):
    """
    Removes outliers based on a sudden unrealistic growth rate followed by a sharp drop.
    
    Parameters:
    df: pandas DataFrame
        DataFrame where each row is a well and each column is a time point.
    max_growth_rate: float
        Maximum allowed growth rate between consecutive time points.
    max_drop_rate: float
        Maximum allowed drop rate between consecutive time points (optional).
        
    Returns:
    filtered_df: pandas DataFrame
        DataFrame with outliers replaced with NaN.
    """
    
    # Create a copy to avoid modifying the original DataFrame
    filtered_df = df.copy()

    # Loop through each well
    for well in df.index:
        # Calculate growth rate (difference between consecutive time points)
        growth_rate = df.loc[well].diff()

        # Flag time points where growth rate exceeds the max_growth_rate
        growth_spikes = growth_rate > max_growth_rate

        # Optionally, flag time points where a sharp drop happens
        drop_spikes = growth_rate < max_drop_rate

        # Replace values where growth spikes or drop spikes occur with NaN
        filtered_df.loc[well][growth_spikes | drop_spikes] = np.nan
    
    return filtered_df

def combined_outlier_detection(df, max_growth_rate=0.5, max_drop_rate=-0.5, mad_threshold=3):
    """
    Combines growth rate-based outlier detection with MAD-based outlier detection for time series.
    
    Parameters:
    df: pandas DataFrame
        DataFrame where each row is a well and each column is a time point.
    max_growth_rate: float
        Maximum allowed growth rate between consecutive time points.
    max_drop_rate: float
        Maximum allowed drop rate between consecutive time points.
    mad_threshold: float
        MAD-based outlier threshold (multiplier for the median absolute deviation).
    
    Returns:
    filtered_df: pandas DataFrame
        DataFrame with outliers replaced with NaN.
    """
    
    # Step 1: Growth rate-based outlier detection
    filtered_df = df.copy()

    for well in df.index:
        # Calculate growth rate (difference between consecutive time points)
        growth_rate = df.loc[well].diff()

        # Flag time points where growth rate exceeds the max_growth_rate
        growth_spikes = growth_rate > max_growth_rate

        # Flag time points where a sharp drop happens
        drop_spikes = growth_rate < max_drop_rate

        # Replace values where growth spikes or drop spikes occur with NaN
        filtered_df.loc[well][growth_spikes | drop_spikes] = np.nan

    # Step 2: MAD-based outlier detection (applied after growth rate filtering)
    for well in df.index:
        well_data = filtered_df.loc[well]

        # Calculate the median and MAD for each well after growth rate filtering
        median_value = np.nanmedian(well_data)
        mad = np.nanmedian(np.abs(well_data - median_value))

        # Set a threshold based on MAD
        mad_threshold_value = mad_threshold * mad

        # Flag sustained outliers (values far from the median)
        outliers = np.abs(well_data - median_value) > mad_threshold_value

        # Replace sustained outliers with NaN
        filtered_df.loc[well][outliers] = np.nan
    
    return filtered_df

def mad_sliding_window(df, window_size=5, mad_threshold=3):
    filtered_df = df.copy()

    # Iterate over each well
    for well in df.index:
        well_data = df.loc[well]
        
        # Iterate over sliding windows
        for start in range(len(well_data) - window_size + 1):
            window_data = well_data[start:start + window_size]

            # Calculate median and MAD within the window
            median_value = np.nanmedian(window_data)
            mad = np.nanmedian(np.abs(window_data - median_value))

            # Flag values as outliers if they exceed the MAD threshold
            mad_threshold_value = mad_threshold * mad
            outliers = np.abs(window_data - median_value) > mad_threshold_value

            # Replace outliers in the filtered DataFrame
            filtered_df.loc[well].iloc[start:start + window_size][outliers] = np.nan

    return filtered_df
def combined_outlier_detection(df, max_growth_rate=0.45, max_drop_rate=-0.45, mad_threshold=2):
    """
    Combines growth rate-based outlier detection with MAD-based outlier detection for time series.
    
    Parameters:
    df: pandas DataFrame
        DataFrame where each row is a well and each column is a time point.
    max_growth_rate: float
        Maximum allowed growth rate between consecutive time points.
    max_drop_rate: float
        Maximum allowed drop rate between consecutive time points.
    mad_threshold: float
        MAD-based outlier threshold (multiplier for the median absolute deviation).
    
    Returns:
    filtered_df: pandas DataFrame
        DataFrame with outliers replaced with NaN.
    """
    
    # Step 1: Growth rate-based outlier detection
    filtered_df = df.copy()

    for well in df.index:
        # Calculate growth rate (difference between consecutive time points)
        growth_rate = df.loc[well].diff()

        # Flag time points where growth rate exceeds the max_growth_rate (sudden increases)
        growth_spikes = growth_rate > max_growth_rate
        
        # For drops: flag the time point just before a sharp decrease
        drop_spikes = (growth_rate < max_drop_rate).shift(-1)  # Shift by -1 to flag the previous time point

        # Replace values where growth spikes or the time point just before a sharp drop occur with NaN
        filtered_df.loc[well][growth_spikes | drop_spikes] = np.nan

    # Step 2: MAD-based outlier detection (applied after growth rate filtering)
    for well in df.index:
        well_data = filtered_df.loc[well]

        # Calculate the median and MAD for each well after growth rate filtering
        median_value = np.nanmedian(well_data)
        mad = np.nanmedian(np.abs(well_data - median_value))

        # Set a threshold based on MAD
        mad_threshold_value = mad_threshold * mad

        # Flag sustained outliers (values far from the median)
        outliers = np.abs(well_data - median_value) > mad_threshold_value

        # Replace sustained outliers with NaN
        filtered_df.loc[well][outliers] = np.nan
    
    return filtered_df

def apply_rolling_median_filter(df, window_size=3):
    """
    Applies a rolling median filter to smooth the data for each well.
    
    Parameters:
    df: pandas DataFrame
        DataFrame where each row is a well and each column is a time point.
    window_size: int
        The size of the rolling window to calculate the median.
    
    Returns:
    filtered_df: pandas DataFrame
        DataFrame with the rolling median filter applied.
    """
    # Apply rolling median filter to each well (each row)
    filtered_df = df.apply(lambda row: row.rolling(window=window_size, center=True).median(), axis=1)
    
    return filtered_df

def hampel_filter(df, window_size=5, n_sigma=3):
    """
    Applies a Hampel filter to smooth the data for each well.
    
    Parameters:
    df: pandas DataFrame
        DataFrame where each row is a well and each column is a time point.
    window_size: int
        The size of the rolling window to calculate the median and MAD.
    n_sigma: int
        The number of standard deviations (scaled MAD) to use as a threshold for outlier detection.
    
    Returns:
    filtered_df: pandas DataFrame
        DataFrame with outliers replaced by the median of the window.
    """
    # Create a copy of the DataFrame to store the filtered data
    filtered_df = df.copy()

    # Apply the Hampel filter to each well
    for well in df.index:
        well_data = df.loc[well]

        # Apply rolling window to calculate median and MAD
        rolling_median = well_data.rolling(window=window_size, center=True).median()
        mad = lambda x: np.median(np.abs(x - np.median(x)))  # MAD function
        rolling_mad = well_data.rolling(window=window_size, center=True).apply(mad)

        # Hampel filter: replace values where deviation from median exceeds threshold
        threshold = n_sigma * rolling_mad
        outliers = np.abs(well_data - rolling_median) > threshold
        filtered_df.loc[well][outliers] = rolling_median[outliers]

    return filtered_df

def main():
    
    
    '''
    
    default values
    python growth_processing.py --path "/Volumes/EVE/20241010 HPLC" --output /Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/growth/20241010_hplc_growth_experiment.txt
    python growth_processing.py --path "/Volumes/EVE/20241010 HPLC" --output "../experiments/growth/20241010_hplc_growth_experiment.txt"
    python growth_processing.py --path "/Volumes/EVE/20241010 HPLC" --output "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/amiga_results/data/20241010_hplc_growth_experiment.txt " --blank distance
    
    python growth_processing.py --path "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/growth/20241010 HPLC" --output "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/amiga_results/data/20241010_hplc_growth_experiment.txt" --correction "yes" --n 3 --filtering "yes"
    
    /Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/amiga_results/data/20241010_hplc_growth_experiment.txt 
    
    Then run:
    
    python "../amiga-master/amiga.py" summarize -i ../experiments/amiga_results/data/20241010_hplc_growth_experiment.txt
    
    '''
    
    parser = argparse.ArgumentParser(description='Script with a command-line argument.')
    parser.add_argument('--path', type=str, help='Value for the "path" variable.')
    parser.add_argument('--output', type=str, help='Value for the "output" variable.')
    parser.add_argument('--correction', type=str, help='Value for the "correction" variable.')
    parser.add_argument('--filtering', type=str, help='Value for the "filtering" variable.')
    parser.add_argument('--n', type=str, help='Value for the "n" variable.')
   
    args = parser.parse_args()
         
    if args.path is not None:
        path = args.path
        #print(f'Variable "target" set to: {target}')
    else:
        print('Variable "path" not provided.')
        
        
    if args.output is not None:
        output = args.output
        #print(f'Variable "target" set to: {target}')
    else:
        print('Variable "output" not provided.')
        
    if args.correction is not None:
        correction = args.correction
        #print(f'Variable "target" set to: {target}')
    else:
        print('Variable "correction" not provided.')
    
    if args.n is not None:
        n = int(args.n)
        #print(f'Variable "target" set to: {target}')
    else:
        print('Variable "n" not provided.')
 
    if args.filtering is not None:
        filtering = args.filtering
        #print(f'Variable "target" set to: {target}')
    else:
        print('Variable "filtering" not provided.')
    
    # Define the directory where your files are stored
    #directory_path = '/Volumes/EVE/20241010 HPLC'
    directory_path = path
    
    # Get a list of all files starting with "AUTOMATED" in the directory
    files = glob.glob(os.path.join(directory_path, "AUTOMATED*"))
    
    # Sort the files by modification date (earliest first)
    files.sort(key=lambda x: os.path.getmtime(x))
    
    # Initialize an empty dictionary to hold the data
    data_dict = {}
    
    # List of valid keys from A01 to H12
    valid_keys = [f"{row}{col:02}" for row in 'ABCDEFGH' for col in range(1, 13)]
    
    # Read each file and extract the values
    for file_path in files:
        print(f"Processing file: {file_path}")
        try:
            with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
                
                inside_data_section = False
                for line in file:
                    
                    # Check if we're in the measurement data section
                    if "Measurement Data" in line:
                        inside_data_section = True
                        continue
                    
                    # Skip irrelevant lines from the OD file
                    if not inside_data_section or line.strip() == "" or "=" in line:
                        continue
    
                    # Extract key-value pairs
                    parts = line.split(":")
                    if len(parts) == 2:
                        key = parts[0].strip()  # e.g. "A01"
                        value = float(parts[1].strip())  # e.g. 0.1817
    
                        # Only add valid keys (A01-H12) to the dictionary
                        if key in valid_keys:
                            if key not in data_dict:
                                data_dict[key] = []
                            data_dict[key].append(value)
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
    
    
    data = pd.DataFrame.from_dict(data_dict, orient = 'index', columns = [i * 1200 for i in range(len(files))])
    #data.to_csv('../experiments/growth/20241010_hplc_growth_experiment.txt', sep='\t')
    blank_wells = ['A01','B01','C01','D01','E01','F01','G01','H01', 
                   'A07','B07','C07','D07','E07','F07','G07','H07',
                   'A01','B12','C12','D12','E12','F12','G12','H12',]  # Replace with your list of blank wells
    
    if correction == 'yes':
    
        
        # This needs to be loaded in somehow when I generate the experimental design. Should not be too hard
        
        
        # Apply the normalization function
        data = subtract_closest_blanks(data, blank_wells, n)
        #print(data.iloc[0,:])
        #data.to_csv(output, sep='\t')
   
    if filtering == 'yes':
       #data, outliers = remove_timepoint_outliers_mad(data, blank_wells, threshold = 2.0)
       
       #data = remove_outliers_per_well_mad(data, threshold=2)
       
       
       datatemp = remove_outliers_based_on_growth_rate(data, max_growth_rate=0.5*0.34, max_drop_rate=-0.5*0.34)
       #data = combined_outlier_detection(data, max_growth_rate=0.45/0.34, max_drop_rate=-0.45/0.34, mad_threshold=1.5)
       #data = combined_outlier_detection(data, max_growth_rate=0.45, max_drop_rate=-0.45, mad_threshold=1.5)
       
       #data = pd.DataFrame.from_dict(data_dict, orient = 'index', columns = [i * 1200 for i in range(len(files))])
       #data = apply_mad_on_residuals(df, mad_threshold=1.5)
       #data = mad_sliding_window(data, window_size=10, mad_threshold=2)
       
       #data = apply_rolling_median_filter(data, window_size=5)
       #tempdata = hampel_filter(data, window_size=3, n_sigma=2)
       #tempdata = kalman_smoothing(data)
       #tempdata = iterative_mad_outlier_removal(data, mad_threshold=2, max_iterations=5)
       
       
       tempdata = apply_rolling_median_filter(data, window_size=3)
       tempdata = remove_outliers_based_on_growth_rate(tempdata, max_growth_rate=0.5*0.34, max_drop_rate=-0.5*0.34)
       
       #data, outliers = remove_timepoint_outliers_derivative(data[data.index.isin(blank_wells)], max_change=0.5)
       #max_change=0.5
    
    data.to_csv(output, sep='\t')
    #data, outliers = remove_timepoint_outliers_derivative(data, max_change=0.5)

def apply_rolling_median_filter(df, window_size=3):
    return df.apply(lambda row: row.rolling(window=window_size, center=True).median(), axis=1)

def remove_outliers_based_on_growth_rate(df, max_growth_rate=0.5, max_drop_rate=-0.5):
    filtered_df = df.copy()
    for well in df.index:
        growth_rate = df.loc[well].diff()
        growth_spikes = growth_rate > max_growth_rate
        drop_spikes = (growth_rate < max_drop_rate).shift(-1)  # Flag the time point before the drop
        filtered_df.loc[well][growth_spikes | drop_spikes] = np.nan
    return filtered_df


import argparse
import os
import glob
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import os
os.chdir('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/scripts')
main()


'''
import scipy
from scipy.stats import bootstrap

from scipy import stats
stats.bootstrap([[1,2,3],[5,6,7]], np.std, confidence_level=0.9,random_state=rng)

res = bootstrap([[1,2,3],[5,6,7]], np.std, confidence_level=0.9,random_state=rng)



def OD_fitter(x):
    return -0.0635*(x**3)+0.430*(x**2)+0.552*x+0.04  # Example: a linear transformation

# Create time points for the x-axis
num_files = len(files)
time_points = [i * (1/3) for i in range(num_files)]  # Start at 0, increment by 1/3 for each file

# Dictionary for averaging specific keys into groups
average_groups = {
    'Blank': ['A01','B01','C01','D01','E01','F01','G01','H01','A12','B12','C12','D12','E12','F12','G12','H12'],
    'his3-del': generate_list('A', 'A', 2, 11)+generate_list('B', 'B', 2, 11)+
    generate_list('C', 'C', 2, 11)+generate_list('D', 'D', 2, 11)+
    generate_list('E', 'E', 2, 11)+generate_list('F', 'F', 2, 11)+
    generate_list('G', 'G', 2, 11)+generate_list('H', 'H', 2, 11)
}

# Create dictionaries to store averaged data and standard deviations
averaged_data = {}
std_devs = {}

# Compute the averages and standard deviations across time for each group
for group, keys in average_groups.items():
    avg_values = []
    std_values = []
    for i in range(num_files):
        # Collect values for the given time point across the specified keys
        #key_values = [data_dict[key][i] for key in keys if key in data_dict]
        key_values = [OD_fitter(data_dict[key][i]) for key in keys if key in data_dict]
        
        if key_values:
            avg_values.append(np.mean(key_values))
            std_values.append(np.std(key_values))  # Standard deviation
    averaged_data[group] = avg_values
    std_devs[group] = std_values

# Plotting the data
fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(10, 8))  # Create subplots

# --- Subplot 1: Original Data ---
for key, values in data_dict.items():
    ax1.plot(time_points, values, label=key)

ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('OD600')
ax1.set_title('Growth Profiling (all wells)')
ax1.grid(True)

# --- Subplot 2: Averaged Data for Groups with Standard Deviation Fill ---
for group, avg_values in averaged_data.items():
    # Get the corresponding standard deviations
    std_values = std_devs[group]
    
    # Plot the average values
    ax2.plot(time_points, avg_values, label=group)
    
    # Plot the fill between (mean - std) and (mean + std)
    ax2.fill_between(time_points, 
                     np.array(avg_values) - np.array(std_values), 
                     np.array(avg_values) + np.array(std_values), 
                     alpha=0.3)  # Standard deviation shading

ax2.set_xlabel('Time (hours)')
ax2.set_ylabel('Average OD600')
ax2.set_title('Growth Profiling (his3-del vs. blank)')
ax2.legend(title="Group", bbox_to_anchor=(1.05, 1), loc='upper left')
ax2.grid(True)

# Adjust layout and show the plots
plt.tight_layout()
plt.show()


import matplotlib.gridspec as gridspec

# Create a gridspec layout
fig = plt.figure(figsize=(12, 6))  # Define overall figure size (wider than tall)
gs = gridspec.GridSpec(2, 2, width_ratios=[1.5, 1])  # Two rows, two columns, with left plot wider

# --- Subplot 1: Original Data (on the left side, taking both rows) ---
ax1 = plt.subplot(gs[:, 0])  # Use all rows for ax1 in the first column
for key, values in data_dict.items():
    ax1.plot(time_points, values, label=key)

ax1.set_xlabel('Time (hours)')
ax1.set_ylabel('OD600')
ax1.set_title('Growth Profiling (all wells)')
ax1.grid(True)

# --- Subplot 2: Averaged Data for Groups with Standard Deviation Fill (top-right) ---
ax2 = plt.subplot(gs[0, 1])  # Top-right plot
for group, avg_values in averaged_data.items():
    # Get the corresponding standard deviations
    std_values = std_devs[group]
    
    # Plot the average values
    ax2.plot(time_points, avg_values, label=group)
    
    # Plot the fill between (mean - std) and (mean + std)
    ax2.fill_between(time_points, 
                     np.array(avg_values) - np.array(std_values), 
                     np.array(avg_values) + np.array(std_values), 
                     alpha=0.3)  # Standard deviation shading

ax2.set_xlabel('Time (hours)')
ax2.set_ylabel('Average OD600')
ax2.set_title('Growth Profiling (his3-del vs. blank)')
ax2.legend(title="Group", bbox_to_anchor=(1.05, 1), loc='upper left')
ax2.grid(True)

# --- Subplot 3: Log-Scale Averaged Data (bottom-right) ---
ax3 = plt.subplot(gs[1, 1])  # Bottom-right plot
for group, avg_values in averaged_data.items():
    # Get the corresponding standard deviations
    std_values = std_devs[group]
    
    # Plot the average values in log scale
    ax3.plot(time_points, avg_values, label=group)
    
    # Plot the fill between (mean - std) and (mean + std) in log scale
    ax3.fill_between(time_points, 
                     np.array(avg_values) - np.array(std_values), 
                     np.array(avg_values) + np.array(std_values), 
                     alpha=0.3)  # Standard deviation shading

ax3.set_xlabel('Time (hours)')
ax3.set_ylabel('Average OD600 (Log Scale)')
ax3.set_title('Growth Profiling (his3-del vs. blank, Log Scale)')
ax3.set_yscale('log')  # Set y-axis to log scale
ax3.grid(True)

# Adjust layout and show the plots
plt.tight_layout()
plt.show()



data = pd.DataFrame.from_dict(data_dict, orient = 'index', columns = [i * 1200 for i in range(58)])
data.to_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/relStack/experiments/experiment/data/hplc_growth_experiment.txt', sep='\t')





def extract_17th_item(dict_of_lists):
    result = {}
    for key, value_list in dict_of_lists.items():
        if isinstance(value_list, list) and len(value_list) >= 17:
            result[key] = value_list[16]  # Indexing starts at 0, so 16 is the 17th item
        else:
            result[key] = None  # If the list has fewer than 17 items
    return result

# Example dictionary

import pandas as pd
data = pd.DataFrame.from_dict(extract_17th_item(data_dict), orient = 'index', columns = ['Values'])



import pandas as pd
import matplotlib.pyplot as plt

# Generate the indices A01 to H12
rows = [chr(i) for i in range(ord('A'), ord('H') + 1)]  # Letters from A to H
cols = [f'{i:02d}' for i in range(1, 13)]  # Numbers from 01 to 12

# Create the DataFrame with random values for demonstration
index_labels = [f'{row}{col}' for row in rows for col in cols]
data = pd.DataFrame(index=index_labels, data={'Values': range(len(index_labels))})

# Reshape the data for plotting on a grid (8 rows, 12 columns)
data_grid = data['Values'].values.reshape((8, 12))

# Plot the DataFrame values as a heatmap
plt.figure(figsize=(10, 6))
plt.imshow(data_grid, cmap='viridis', aspect='auto')
plt.colorbar(label='Values')

# Set labels for the axes
plt.xticks(ticks=range(12), labels=cols)
plt.yticks(ticks=range(8), labels=rows)
plt.xlabel('Column')
plt.ylabel('Row')
plt.title('Grid Plot of OD600 values')

# Show the plot
plt.show()
'''