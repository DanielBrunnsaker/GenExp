#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  9 10:26:05 2025

@author: danbru
"""


import numpy as np
import os
import glob
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd

def unify_well_format(well_name: str) -> str:
    """
    Convert strings like 'A01' -> 'A1', 'C08' -> 'C8', etc.
    If the numeric part is all zeros, revert to '0'.
    """
    if not well_name:
        return well_name  # edge case handling
    
    letter = well_name[0]               # e.g. 'A'
    number_str = well_name[1:].lstrip('0')  # e.g. '01' -> '1'
    
    if number_str == "":               # if it was all zeros
        number_str = "0"
    
    return f"{letter}{number_str}"      # e.g. 'A1'

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

def remove_outliers_by_group(
    df_growth: pd.DataFrame, 
    df_meta: pd.DataFrame,
    well_col: str = "well", 
    group_col: str = "Summary", 
    z_thresh: float = 3.0
) -> pd.DataFrame:
    """
    1) Standardize wells in both dataframes so they match ('A01' -> 'A1'),
    2) Join growth data to metadata (so each well has a group),
    3) For each group, remove outliers based on row-wise mean and std z-scores.
    """

    # ----------------------------------------------------------------
    # STEP A: Unify well formats in BOTH DataFrames
    # ----------------------------------------------------------------
    
    # If df_growth uses the well as an index:
    df_growth = df_growth.copy()
    df_growth.index = df_growth.index.map(lambda x: unify_well_format(str(x)))
    
    # If df_meta has well in a column named well_col (e.g. "well"):
    df_meta = df_meta.copy()
    df_meta[well_col] = df_meta[well_col].astype(str).apply(unify_well_format)

    # ----------------------------------------------------------------
    # STEP B: Merge (Join) on the well name
    # ----------------------------------------------------------------
    # Make df_meta indexed by well so we can join easily
    df_meta.set_index(well_col, inplace=True)

    # Join growth data (index = wells) with meta (index = wells)
    df_combined = df_growth.join(df_meta, how="left")  
    # Now df_combined has numeric columns (time points) + the group_col (e.g. 'Summary').

    # ----------------------------------------------------------------
    # STEP C: Define a simple outlier filter for a single group's data
    # ----------------------------------------------------------------
    def filter_one_group(group_df: pd.DataFrame) -> pd.DataFrame:
        """
        Within the group_df, remove rows that deviate beyond z_thresh
        in row-wise mean or row-wise std.
        """
        # Numeric-only columns (exclude 'Summary', etc.)
        numeric_cols = group_df.select_dtypes(include=[np.number])
        
        # Row-wise stats
        row_means = numeric_cols.mean(axis=1)
        row_stds  = numeric_cols.std(axis=1)
        
        # Group-level stats (of row_means and row_stds)
        mean_of_means = row_means.mean()
        std_of_means  = row_means.std() + 1e-9
        mean_of_stds  = row_stds.mean()
        std_of_stds   = row_stds.std() + 1e-9
        
        # Keep only rows whose mean and std are within z_thresh standard deviations
        mean_mask = (np.abs(row_means - mean_of_means) / std_of_means) < z_thresh
        std_mask  = (np.abs(row_stds  - mean_of_stds ) / std_of_stds ) < z_thresh

        return group_df[mean_mask & std_mask]

    # ----------------------------------------------------------------
    # STEP D: Group by group_col and filter outliers in each group
    # ----------------------------------------------------------------
    filtered_subdfs = []
    for group_name, subdf in df_combined.groupby(group_col):
        subdf_inliers = filter_one_group(subdf)
        filtered_subdfs.append(subdf_inliers)

    # Re-combine
    df_filtered = pd.concat(filtered_subdfs).sort_index()

    return df_filtered

def plot_group_averages_in_hours(df: pd.DataFrame, group_col: str = "Summary"):
    """
    Plot publication-ready average growth curves by group, converting time from seconds to hours.
    
    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing:
         - numeric time columns (e.g. 0, 1200, 2400, ...)
         - a column named `group_col` specifying each well's group
    group_col : str
        The column in df that indicates the group each row belongs to.
    """
    # ----------------------------------------------------------------
    # STEP 1: Identify numeric columns (time points)
    # ----------------------------------------------------------------
    numeric_cols = df.select_dtypes(include=["float", "int"]).columns
    numeric_cols = numeric_cols.drop(group_col, errors="ignore")  # in case group_col is numeric dtype

    # ----------------------------------------------------------------
    # STEP 2: Compute the mean of each numeric column, grouped by group_col
    #         (Each row is now a group, columns are time points.)
    # ----------------------------------------------------------------
    df_avg = df.groupby(group_col)[numeric_cols].mean().reset_index()

    # ----------------------------------------------------------------
    # STEP 3: Reshape (melt) so that each row is (group, time_in_sec, value)
    #         Then convert to hours.
    # ----------------------------------------------------------------
    df_melt = df_avg.melt(id_vars=group_col, var_name="Time_Seconds", value_name="Mean Growth")

    # Convert melted time column from string to numeric
    df_melt["Time_Seconds"] = pd.to_numeric(df_melt["Time_Seconds"], errors="coerce")

    # Sort by time so lines plot left to right
    df_melt.sort_values(by="Time_Seconds", inplace=True)

    # --- CONVERT SECONDS TO HOURS ---
    df_melt["Time_Hours"] = df_melt["Time_Seconds"] / 3600.0

    # ----------------------------------------------------------------
    # STEP 4: Plot with Seaborn
    # ----------------------------------------------------------------
    sns.set_style("whitegrid")
    plt.figure(figsize=(8, 5))

    # Draw one line per group, using Time_Hours as x-axis
    ax = sns.lineplot(
        data=df_melt,
        x="Time_Hours",
        y="Mean Growth",
        hue=group_col,
        #marker="o",
        dashes=False
    )

    # Decorate plot
    ax.set_title("Mean Growth Curves by Group (Hours)", fontsize=14, weight="bold")
    ax.set_xlabel("Time (hours)", fontsize=12)
    ax.set_ylabel("Mean Growth (OD)", fontsize=12)

    # Legend inside the plot
    ax.legend(title=group_col, loc="best", fontsize=10, title_fontsize=11)

    # Remove top and right spines to get a clean look
    sns.despine()

    plt.tight_layout()
    plt.show()


import numpy as np
import pandas as pd
from scipy.stats import zscore

def filter_outlier_growth_curves(df: pd.DataFrame, group_col: str = "Summary", method: str = "mad", threshold: float = 3.0) -> pd.DataFrame:
    """
    Remove entire growth curves that are outliers within their experimental group.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame where rows = samples, columns = time points + experimental conditions.
    group_col : str
        The column representing experimental groups (e.g., "Summary").
    method : str
        "zscore" for standard deviation-based filtering.
        "iqr" for interquartile range-based filtering.
        "mad" for median absolute deviation filtering.
    threshold : float
        The threshold for detecting outliers (default: 3 for z-score, 1.5 for IQR).

    Returns
    -------
    pd.DataFrame
        A filtered DataFrame with outlier growth curves removed.
    """
    # Identify numeric (time-series) columns
    numeric_cols = df.select_dtypes(include=["float", "int"]).columns

    # Compute summary metric per curve (e.g., AUC, total OD)
    df['growth_summary'] = df[numeric_cols].sum(axis=1)  # AUC-like metric

    def detect_outliers(group):
        if len(group) < 3:  # Not enough data to detect outliers
            group['outlier'] = False
            return group

        if method == "zscore":
            group['outlier'] = np.abs(zscore(group['growth_summary'])) > threshold
        elif method == "iqr":
            Q1, Q3 = group['growth_summary'].quantile([0.25, 0.75])
            IQR = Q3 - Q1
            lower_bound, upper_bound = Q1 - threshold * IQR, Q3 + threshold * IQR
            group['outlier'] = (group['growth_summary'] < lower_bound) | (group['growth_summary'] > upper_bound)
        elif method == "mad":
            median_val = group['growth_summary'].median()
            mad = np.median(np.abs(group['growth_summary'] - median_val))
            modified_z = 0.6745 * (group['growth_summary'] - median_val) / (mad + 1e-9)
            group['outlier'] = np.abs(modified_z) > threshold
        return group

    # Apply outlier detection within each experimental group
    df = df.groupby(group_col, group_keys=False).apply(detect_outliers)

    # Keep only non-outliers
    df_filtered = df[df['outlier'] == False].drop(columns=['growth_summary', 'outlier'])

    return df_filtered


import pandas as pd


import pandas as pd
import numpy as np

def smooth_growth_curves(
    df: pd.DataFrame, 
    window: int = 3, 
    center: bool = True, 
    method: str = "mirror"
) -> pd.DataFrame:
    """
    Apply a rolling average smoothing to each row's time series data,
    while minimizing unrealistic drop-offs at the curve ends.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame where rows = wells, columns = time points.
    window : int
        Size of the moving window for the rolling average.
    center : bool
        If True, the window is centered around each point. If False, the window is trailing.
    method : str
        Method to handle edges ("mirror", "repeat", "none").
        - "mirror": Reflects the edge values to prevent drop-off.
        - "repeat": Repeats the last value to smooth the end.
        - "none": No special handling (default Pandas behavior).

    Returns
    -------
    pd.DataFrame
        A copy of df with numeric (time) columns smoothed.
    """
    numeric_cols = df.select_dtypes(include=["float", "int"]).columns
    df_smoothed = df.copy()

    def smooth_row_values(row):
        values = row.values
        n = len(values)

        # Edge handling - extend the series before smoothing
        if method == "mirror":
            extended = np.concatenate((values[:window][::-1], values, values[-window:][::-1]))
        elif method == "repeat":
            extended = np.concatenate(([values[0]] * (window // 2), values, [values[-1]] * (window // 2)))
        else:
            extended = values  # No edge handling

        smoothed = pd.Series(extended).rolling(window=window, center=center, min_periods=1).mean().values
        
        # Ensure smoothed row matches original length
        start_idx = (len(smoothed) - n) // 2
        smoothed = smoothed[start_idx:start_idx + n]

        return pd.Series(smoothed, index=row.index)  # Ensure Pandas aligns correctly

    # Apply smoothing
    df_smoothed[numeric_cols] = df_smoothed[numeric_cols].apply(smooth_row_values, axis=1)

    return df_smoothed

def rolling_mad_outlier_removal(df: pd.DataFrame, window: int = 5, threshold: float = 3.0) -> pd.DataFrame:
    """
    Detect and remove outliers using a rolling median and Median Absolute Deviation (MAD).

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame where rows = wells, columns = time points.
    window : int
        Size of the rolling window for median and MAD calculation.
    threshold : float
        The threshold for MAD-based outlier detection.

    Returns
    -------
    pd.DataFrame
        A DataFrame with outliers replaced by NaN.
    """
    numeric_cols = df.select_dtypes(include=["float", "int"]).columns
    df_filtered = df.copy()

    def mad_based_outlier_detection(series):
        rolling_median = series.rolling(window, center=True, min_periods=1).median()
        mad = (series - rolling_median).abs().rolling(window, center=True, min_periods=1).median()
        modified_z_score = 0.6745 * (series - rolling_median) / (mad + 1e-9)  # Avoid division by zero
        
        # Convert np.where output to a pandas Series to maintain the correct shape
        return pd.Series(np.where(np.abs(modified_z_score) > threshold, np.nan, series), index=series.index)

    # Apply function row-wise
    df_filtered[numeric_cols] = df_filtered[numeric_cols].apply(mad_based_outlier_detection, axis=1)

    return df_filtered





# Define the directory where your files are stored
#directory_path = '/Volumes/EVE/20241010 HPLC'
directory_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/lactic acid/run/proline_0.5_5_20241219_1102/results/growth/20241219 LLM'
directory_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/spermine3/glutamate_0.0_1_20250128_1618/results/raw'
directory_path = '/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/caffeine/selected/arginine_202502131014/results/growth/raw'


# Get a list of all files starting with "AUTOMATED" in the directory
files = glob.glob(os.path.join(directory_path, "AUTOMATED*"))

# Sort the files by modification date (earliest first)
files.sort(key=lambda x: os.path.getmtime(x))

# Initialize an empty dictionary to hold the data
data_dict = {}

# List of valid keys from A01 to H12
valid_keys = [f"{row}{col:02}" for row in 'ABCDEFGH' for col in range(1, 13)]

import os
import datetime

data_dict = {}
file_times = []  # List to store file times for sorting

for file_path in files:
    print(f"Processing file: {file_path}")
    try:
        with open(file_path, 'r', encoding='utf-8', errors='ignore') as file:
            inside_data_section = False
            file_time = None  # Initialize file time variable

            for line in file:
                # Extract the date and time from the file
                if "Date:" in line and "Time:" in line:
                    parts = line.split("Time:")
                    date_part = parts[0].split("Date:")[1].strip()
                    time_part = parts[1].strip()

                    # Convert to a datetime object for sorting
                    file_time = datetime.datetime.strptime(f"{date_part} {time_part}", "%d/%m/%Y %H:%M:%S")

                # Check if we're in the measurement data section
                if "Measurement Data" in line:
                    inside_data_section = True
                    continue
                
                # Skip irrelevant lines
                if not inside_data_section or line.strip() == "" or "=" in line:
                    continue

                # Extract key-value pairs
                parts = line.split(":")
                if len(parts) == 2:
                    key = parts[0].strip()
                    value = float(parts[1].strip())

                    # Only add valid keys (A01-H12) to the dictionary
                    if key in valid_keys:
                        if key not in data_dict:
                            data_dict[key] = []
                        data_dict[key].append((file_time, value))

            # Store file time for sorting later
            if file_time:
                file_times.append(file_time)

    except Exception as e:
        print(f"Error reading file {file_path}: {e}")

# Sort the dictionary based on time
for key in data_dict:
    data_dict[key].sort(key=lambda x: x[0])  # Sort by datetime

# Convert values back to just floats after sorting
for key in data_dict:
    data_dict[key] = [value for _, value in data_dict[key]]

print("Data sorted by time.")




data = pd.DataFrame.from_dict(data_dict, orient = 'index', columns = [i * 1200 for i in range(len(files))])

layout = pd.read_excel('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/caffeine/selected/arginine_202502131014/protocol/hamilton/pipetting_layout.xlsx')

blank_wells = layout['well'][layout['Summary'] == 'Media control']
blank_wells = blank_wells.apply(lambda x: f"{x[0]}{int(x[1:]):02d}")

# Subtract blanks
blanked_data = subtract_closest_blanks(data, blank_wells, 5)
blanked_data[blanked_data < 0] = 0.01

blanked_data = blanked_data.copy()
blanked_data.index = blanked_data.index.map(lambda x: unify_well_format(str(x)))

# If df_meta has well in a column named well_col (e.g. "well"):
layout = layout.copy()
layout['well'] = layout['well'].astype(str).apply(unify_well_format)
blanked_data = blanked_data.merge(layout[['well','Summary']], left_index = True, right_on = 'well').set_index('well')

df_filtered = filter_outlier_growth_curves(blanked_data, group_col="Summary", method="mad", threshold=3.0)
df_smoothed = smooth_growth_curves(df_filtered, window=5, method="mirror")  # Best for avoiding drop-off
plot_group_averages_in_hours(df_smoothed, group_col="Summary")


# Save df_smoothed and use with amiga?
df_smoothed = df_smoothed.reset_index()
df_smoothed.rename(columns={ df_smoothed.columns[0]: "Well" }, inplace = True)
df_smoothed.set_index('Well', inplace = True)

# Save the blanked one so we can use it for the summary
blanked_data = blanked_data.reset_index()
blanked_data.rename(columns={ blanked_data.columns[0]: "Well" }, inplace = True)
blanked_data.set_index('Well', inplace = True)

#df_smoothed.drop('Summary', axis = 1).to_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/spermine3/glutamate_0.0_1_20250128_1618/results/data/growth_data_corrected.txt',sep = '\t')

blanked_data.reset_index().rename(columns={'Well':'Time'}).set_index('Time').to_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/caffeine/selected/arginine_202502131014/results/growth/processed/growth_data_corrected_forsummary.txt',sep = '\t')
df_smoothed.reset_index().rename(columns={'Well':'Time'}).set_index('Time').drop('Summary', axis = 1).to_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/caffeine/selected/arginine_202502131014/results/growth/processed/growth_data_corrected.txt',sep = '\t')


'''

    python "amiga-master/amiga.py" summarize -i "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/lactic acid/run/proline_0.5_5_20241219_1102/results/data/lactic_acid_smoothed.txt"
    python "amiga-master/amiga.py" fit -i "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/lactic acid/run/proline_0.5_5_20241219_1102/results/data/lactic_acid_smoothed.txt" --skip-first-n 12

'''


