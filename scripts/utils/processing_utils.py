#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 11:41:43 2025

@author: danbru
"""

import os
import glob
import datetime
import numpy as np
import pandas as pd
from scipy.stats import zscore
from croissance import process_curve
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D

from sklearn.isotonic import IsotonicRegression
from scipy.interpolate import PchipInterpolator, UnivariateSpline
import statsmodels.api as sm

def extract_growth_rates(layout, data):
    
    mu_dict = {}
    timing_dict = {}
    for well in data.index:
        
        if well in list(layout['well'][layout['Summary'] == 'Media control']):
            mu = 0
            start = 0
            end = 0
        else:
            try:
                temp_curve = data.iloc[:,18:-1].loc[well]
                temp_curve.index = np.array(data.loc[well].index[18:-1])/3600
                results = process_curve(temp_curve, n0 = 0.1)
                mu = results.growth_phases[0][2]
                start = results.growth_phases[0][0]
                end = results.growth_phases[0][1]
            except:
                print(f'No growth rate could be calculated for well {well}')
                start = 0
                end = 0
                mu = 0
        
        mu_dict[well] = [mu, start, end]
        
    mu_df = pd.DataFrame.from_dict(mu_dict, orient = 'index', columns = ['mu', 'start', 'end'])
    return mu_df


def extract_finalOD(layout, data):
    
    od_dict = {}
    for well in data.index:
        
        if well in list(layout['well'][layout['Summary'] == 'Media control']):
            max_od = 0
        else:
            max_od = data.loc[well].drop('Summary').max()
            
        od_dict[well] = max_od
        
    od_df = pd.DataFrame.from_dict(od_dict, orient = 'index', columns = ['MaxOD'])
    return od_df



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


def process_measurement_data(directory_path):
    """
    Process measurement data from files starting with 'AUTOMATED' in the given folder.

    Parameters:
        directory_path (str): Path to the folder containing the files.

    Returns:
        dict: A dictionary where each valid key (A01-H12) maps to a list of measurement values
              sorted by the file's date and time.
    """
    # Get a list of all files starting with "AUTOMATED" in the directory
    files = glob.glob(os.path.join(directory_path, "AUTOMATED*"))
    
    # Sort the files by modification date (earliest first)
    files.sort(key=lambda x: os.path.getmtime(x))
    
    # List of valid keys from A01 to H12
    valid_keys = [f"{row}{col:02}" for row in 'ABCDEFGH' for col in range(1, 13)]
    
    data_dict = {}
    
    for file_path in files:
        #print(f"Processing file: {file_path}")
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
                        try:
                            value = float(parts[1].strip())
                        except ValueError:
                            continue  # Skip if conversion to float fails
                        
                        # Only add valid keys (A01-H12) to the dictionary
                        if key in valid_keys:
                            if key not in data_dict:
                                data_dict[key] = []
                            data_dict[key].append((file_time, value))
        except Exception as e:
            print(f"Error reading file {file_path}: {e}")
    
    # Sort the dictionary based on time and convert values back to just floats
    for key in data_dict:
        data_dict[key].sort(key=lambda x: x[0])  # Sort by datetime
        data_dict[key] = [value for _, value in data_dict[key]]
    
    #print("Data sorted by time.")
    
    data_df = pd.DataFrame.from_dict(data_dict, orient = 'index', columns = [i * 1200 for i in range(len(files))])
    
    return data_df


def filter_outlier_growth_curves(layout, df, group_col="Summary", method = 'mad', threshold=3, eps=1e-6):
    """
    Filters out wells whose log-AUC or final OD is beyond threshold×MAD
    within each group (e.g. experimental replicate). Uses only the MAD rule.
    
    Parameters
    ----------
    layout : pd.DataFrame
        Must contain columns 'well' and the grouping column (e.g. 'Summary').
    df : pd.DataFrame
        Time-course data, indexed by well, with numeric columns = time points.
    group_col : str
        Column in layout that defines grouping for outlier detection.
    threshold : float
        Multiplier of MAD for outlier detection (default 3.5).
    eps : float
        Small constant to avoid log(0).
    """
    # 1) Copy & unify well names
    df = df.copy()
    df.index = df.index.map(lambda x: unify_well_format(str(x)))
    layout = layout.copy()
    layout['well'] = layout['well'].astype(str).apply(unify_well_format)
    
    # 2) Merge layout info
    df = df.merge(layout[['well', group_col]],
                  left_index=True, right_on='well').set_index('well')
    
    # 3) Identify time-series columns
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    if not numeric_cols:
        raise ValueError("No numeric columns found for time series!")
    
    # temporary, mild, smoothing
    df[numeric_cols] = df[numeric_cols].rolling(axis=1, window=3, center=True, min_periods=1).median()
    
    # 4) Compute AUC (sum approximation) and final OD
    df['AUC']         = df[numeric_cols].sum(axis=1)
    df['log_auc']     = np.log(df['AUC'] + eps)
    final_col         = numeric_cols[-1]
    df['final_value'] = df[final_col]
    df['TV'] = df[numeric_cols].diff(axis=1).abs().sum(axis=1)
    
    def detect_mad_outliers(group):
        # Too few curves → no flags
        if len(group) < 3:
            group['outlier'] = False
            return group
        
        # MAD on log-AUC
        med_auc = np.median(group['log_auc'])
        mad_auc = np.median(np.abs(group['log_auc'] - med_auc))
        out_auc = np.abs(group['log_auc'] - med_auc) > (threshold * mad_auc + eps)
        
        # MAD on final OD
        med_f = np.median(group['final_value'])
        mad_f = np.median(np.abs(group['final_value'] - med_f))
        out_f   = np.abs(group['final_value'] - med_f) > (threshold * mad_f + eps)
        
        # MAD on total variation (spikiness)
        med_tv = np.median(group['TV'])
        mad_tv = np.median(np.abs(group['TV'] - med_tv))
        out_tv = np.abs(group['TV'] - med_tv) > (threshold * mad_tv + eps)
    
        
        #group['outlier'] = out_auc | out_f
        group['outlier'] = out_auc | out_f | out_tv
        return group
    
    # 5) Apply per-group
    df = df.groupby(group_col, group_keys=False).apply(detect_mad_outliers)
    
    # 6) Filter out and drop helper columns
    kept = df.loc[~df['outlier']]
    return kept.drop(columns=['AUC', 'log_auc', 'final_value', 'outlier', 'TV'])


def subtract_and_impute_blanks(layout, data, n_blanks = 3, fillin_value = 0.01, blank_bool = True):
    
    blanked_data = data.copy()
    blank_wells = blanked_data[blanked_data['Summary'] == 'Media control'].index
    blanked_data.iloc[:,:-1] = subtract_closest_blanks(blanked_data.iloc[:,:-1], blank_wells, n_blanks, blank_bool)

    # Subtract blanks
    
    subtracted = blanked_data.iloc[:,:-1].copy()
    subtracted[subtracted < 0] = fillin_value
    blanked_data.iloc[:,:-1] = subtracted
        

    
    return blanked_data


def subtract_closest_blanks(df, blank_wells, n, blank_subtraction):
    
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
    
    result = df.copy()
    for well in df.index:
        
        # Calculate distances from the current well to each blank well. Maybe change how i calculate this distance? 
        # Maybe Manhattah makes more sense?
        distances = [(blank, well_distance(well, blank)) for blank in blank_wells]
        
        # Sort by distance and take the closest n blanks
        closest_blanks = sorted(distances, key=lambda x: x[1])[:n]
        
        
        
        # Get the values of the closest blank wells for each timepoint
        closest_blank_values = df.loc[[blank for blank, _ in closest_blanks]].mean()
        
        if blank_subtraction:
            result.loc[well] -= closest_blank_values
        else:
            result 
    
    return result


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
    
    def smooth_row_values(
        row,
        method="loess",
        window=5,
        loess_frac=0.2,
        spline_s=None,
        enforce_monotonic=False
    ):
        """
        row             : pd.Series of values to smooth
        method          : "mean", "loess", "pchip", or "spline"
        window          : window length for moving average (only if method=="mean")
        loess_frac      : fraction of data used for LOESS (only if method=="loess")
        """
        x = np.arange(len(row))
        y = row.values
    
        loess_sm = sm.nonparametric.lowess(
                y, x,
                frac=loess_frac,
                it=1,               # one robustness pass
                return_sorted=False
            )
        result = loess_sm

    
        return pd.Series(result, index=row.index)

    # Apply smoothing
    df_smoothed[numeric_cols] = df_smoothed[numeric_cols].apply(smooth_row_values, axis=1)

    return df_smoothed

def compute_auc(df):
    # Exclude the last column (experiment descriptor)
    time_series_data = df
    time_points = np.array(time_series_data.columns, dtype=float)  # Convert column names to time points

    auc_results = []
    for well_id, row in time_series_data.iterrows():
        growth_values = row.astype(float).values  # Growth data
        
        # Compute AUC using Trapezoidal Rule
        auc_trapz = np.trapz(growth_values, time_points)
        auc_results.append([well_id, auc_trapz])

    # Create DataFrame with AUC results
    auc_df = pd.DataFrame(auc_results, columns=["Well", "AUC"])
    auc_df.set_index("Well", inplace=True)
    
    return auc_df




# Plotting scripts

def plot_group_averages_in_hours(df, group_col = "Summary"):

    numeric_cols = df.select_dtypes(include=["float", "int"]).columns
    numeric_cols = numeric_cols.drop(group_col, errors="ignore")  # in case group_col is numeric dtype
    df_avg = df.groupby(group_col)[numeric_cols].mean().reset_index()

    df_melt = df_avg.melt(id_vars=group_col, var_name="Time_Seconds", value_name="Mean Growth")

    # Convert melted time column from string to numeric
    df_melt["Time_Seconds"] = pd.to_numeric(df_melt["Time_Seconds"], errors="coerce")

    # Sort by time so lines plot left to right
    df_melt.sort_values(by="Time_Seconds", inplace=True)

    df_melt["Time_Hours"] = df_melt["Time_Seconds"] / 3600.0

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

    #ax.set_ylim([0,2.5])
    ax.set_ylim([0, df.drop(columns='Summary').max().max()+0.25])

    # Remove top and right spines to get a clean look
    sns.despine()

    plt.tight_layout()
    #plt.show()


def get_summary_color_map(df, palette_name="tab10"):
    """
    Generates a color mapping from the unique values in the "Summary" column.
    The groups are sorted alphabetically for consistency.
    
    Parameters:
        df (pd.DataFrame): DataFrame containing a "Summary" column.
        palette_name (str): Name of the seaborn palette to use.
    
    Returns:
        dict: A dictionary mapping each unique group to a color.
    """
    groups = sorted(df["Summary"].unique())
    palette = sns.color_palette(palette_name, len(groups))
    return {group: palette[i] for i, group in enumerate(groups)}

def plot_growth_curves(growth_df, color_map=None, annotations_df=None):
    """
    Plots a grid of growth curves arranged as a 96-well plate.
    
    Parameters:
        growth_df (pd.DataFrame): DataFrame with wells as the index (e.g., "A1", "B3", etc.)
                                  and each column representing a timepoint, except for the "Summary"
                                  column which is used for color coding.
        color_map (dict, optional): A predefined dictionary mapping experimental groups (from
                                  the "Summary" column) to colors. If not provided, one is generated.
        annotations_df (pd.DataFrame, optional): DataFrame with the well as the index and columns:
                                  "start", "end", and "mu". A shaded area is drawn between start
                                  and end times for wells where mu != 0.
    
    The bottom left subplot ("H1") shows numeric tick labels (with "1\nTime (hours)" and "H\nOD600"),
    and a legend based on the color_map is added beneath the grid.
    """
    # Generate a consistent color mapping if one is not provided.
    if color_map is None:
        color_map = get_summary_color_map(growth_df, palette_name="tab10")
    
    # Define expected rows and columns for a 96-well plate.
    rows = list("ABCDEFGH")
    cols = list(range(1, 13))  # 1 to 12
    
    # Identify timepoint columns (all except "Summary").
    time_cols = [col for col in growth_df.columns if col != "Summary"]
    # Convert timepoints from seconds to hours.
    x_vals = [float(tp) / 3600 for tp in time_cols]
    
    # Create a grid of subplots.
    fig, axes = plt.subplots(nrows=8, ncols=12, figsize=(14, 9), sharex=True, sharey=True)
    
    # Loop through each well position.
    for i, row in enumerate(rows):
        for j, col in enumerate(cols):
            well = f"{row}{col}"  # e.g., "A1", "B5", etc.
            ax = axes[i, j]
            
            # Determine the curve color from the "Summary" column.
            color = 'black'
            if well in growth_df.index and "Summary" in growth_df.columns:
                summary_val = growth_df.loc[well, "Summary"]
                color = color_map.get(summary_val, 'black')
            
            # Plot the growth curve using only the timepoint columns.
            if well in growth_df.index:
                ax.plot(x_vals, growth_df.loc[well, time_cols], linestyle='-', color='black', linewidth=2)
            
            # Draw a shaded area for annotated wells if mu != 0.
            if annotations_df is not None and well in annotations_df.index:
                if annotations_df.loc[well]['mu'] != 0:
                    start_time = float(annotations_df.loc[well, "start"])
                    end_time = float(annotations_df.loc[well, "end"])
                    ax.axvspan(start_time, end_time, color=color, alpha=0.4)
            
            # Hide tick labels on all subplots by default.
            ax.tick_params(labelbottom=False, labelleft=False)
            ax.set_ylim([0, growth_df.drop(columns='Summary').max().max()+0.25])
            
            # For the bottom row, add the column label.
            if i == 7:
                ax.set_xlabel(str(col), fontsize=6)
            # For the left column, add the row label.
            if j == 0:
                ax.set_ylabel(row, fontsize=6)
    
    # Enable numeric tick labels only for the bottom left subplot ("H1").
    bottom_left_ax = axes[7, 0]
    bottom_left_ax.tick_params(labelbottom=True, labelleft=True, labelsize=10)
    bottom_left_ax.set_xlabel(f"{cols[0]}\nTime (hours)", fontsize=8)
    bottom_left_ax.set_ylabel(f"{rows[-1]}\nOD600", fontsize=8)
    
    # Add a legend below the grid.
    legend_handles = [Line2D([0], [0], color=color, lw=2.5, label=str(summary))
                      for summary, color in sorted(color_map.items())]
    fig.legend(handles=legend_handles, loc='lower center', ncol=len(color_map),
               fontsize=8, bbox_to_anchor=(0.5, -0.1))
    fig.subplots_adjust(bottom=0.2)
    
    plt.tight_layout()
    #plt.show()

def plot_boxplots(df, color_map=None):
    
    # Generate a consistent color mapping if not provided.
    if color_map is None:
        color_map = get_summary_color_map(df, palette_name="tab10")
    
    # Create a figure with two subplots.
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # Left subplot: AUC boxplot.
    sns.boxplot(x="Summary", y="AUC", data=df, ax=axes[0], palette=color_map)
    axes[0].set_title("AUC")
    axes[0].set_xlabel("")
    axes[0].set_xticklabels([])  # Remove x-axis tick labels.
    axes[0].set_ylabel("AUC")
    
    # Right subplot: Growth rate (mu) boxplot.
    sns.boxplot(x="Summary", y="mu", data=df, ax=axes[1], palette=color_map)
    axes[1].set_title("Growth Rate (mu)")
    axes[1].set_xlabel("")
    axes[1].set_xticklabels([])  # Remove x-axis tick labels.
    axes[1].set_ylabel("mu")
    
    # Create legend handles from the color_map.
    legend_handles = [
        Line2D([0], [0], marker='o', color='w', label=str(group),
               markerfacecolor=color, markersize=10)
        for group, color in sorted(color_map.items())
    ]
    
    # Add the legend beneath the plots.
    fig.legend(handles=legend_handles, loc='lower center',
               ncol=min(len(color_map), 3), fontsize=9, frameon=False)
    fig.subplots_adjust(bottom=0.25)
    
    # plt.show() can be called outside the function if desired.


def plot_violinplots(df, color_map=None):
    """
    Creates side-by-side violin plots for "AUC" and "mu" from a DataFrame,
    using a color mapping based on the "Summary" column. A boxplot is drawn inside each violin.
    
    Parameters:
        df (pd.DataFrame): DataFrame with columns "AUC", "mu", and "Summary" where "Summary"
                           defines the experimental groups.
        color_map (dict, optional): A predefined dictionary mapping groups to colors.
                                    If not provided, one is generated.
    
    A legend based on the color mapping is placed beneath the plots.
    """
    # Generate a consistent color mapping if not provided.
    if color_map is None:
        color_map = get_summary_color_map(df, palette_name="tab10")
    
    # Create a figure with two subplots.
    fig, axes = plt.subplots(1, 2, figsize=(12, 6))
    
    # Left subplot: AUC violin plot with a boxplot inside.
    sns.boxplot(x="Summary", y="AUC", data=df, ax=axes[0],
                   inner="box", palette=color_map)
    axes[0].set_title("AUC")
    axes[0].set_xlabel("")
    axes[0].set_xticklabels([])  # Remove x-axis tick labels.
    axes[0].set_ylabel("AUC")
    
    # Right subplot: Growth rate (mu) violin plot with a boxplot inside.
    sns.boxplot(x="Summary", y="mu", data=df, ax=axes[1],
                   inner="box", palette=color_map)
    axes[1].set_title("Growth Rate (mu)")
    axes[1].set_xlabel("")
    axes[1].set_xticklabels([])  # Remove x-axis tick labels.
    axes[1].set_ylabel("mu")
    
    # Create legend handles from the color_map.
    legend_handles = [Line2D([0], [0], marker='o', color='w', label=str(group),
                             markerfacecolor=color, markersize=10)
                      for group, color in sorted(color_map.items())]
    
    # Add the legend naturally beneath the plots.
    fig.legend(handles=legend_handles, loc='lower center',
               ncol=min(len(color_map), 3), fontsize=9, frameon=False)
    fig.subplots_adjust(bottom=0.25)
    
    #plt.show()
