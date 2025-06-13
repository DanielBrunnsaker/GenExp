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
from matplotlib.patches import Patch
import math


from sklearn.isotonic import IsotonicRegression
from scipy.interpolate import PchipInterpolator, UnivariateSpline
import statsmodels.api as sm

def extract_growth_rates(layout, data):
    """
    Function to extract growth rate from smoothed curves

    Parameters
    ----------
    layout : dataframe
        Dataframe containing the plate layout.
    data : dataframe
        dataframe of growth data (over time).

    Returns
    -------
    mu_df : dataframe
        Dataframe containing mu, and start/stop of the calculation interval.

    """
    
    mu_dict = {}
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
    """ 
    Extracts the final OD
    """
    
    od_dict = {}
    for well in data.index:
        
        if well in list(layout['well'][layout['Summary'] == 'Media control']):
            max_od = 0
        else:
            max_od = data.loc[well].drop('Summary').max()
            
        od_dict[well] = max_od
        
    od_df = pd.DataFrame.from_dict(od_dict, orient = 'index', columns = ['FinalOD'])
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
    process measurement files from the polarstar omega, and reshape it into a 
    dict, and then a dataframe.

    Parameters
    ----------
    directory_path : path
        path to the polarstar raw files.

    Returns
    -------
    data_df : dataframe
        dataframe with temporal growth data for all wells.

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
                file_time = None  
                
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


def filter_outlier_growth_curves(layout, df, group_col="Summary", threshold=3):
    """
    Filters wells based on AUC, TV and final OD. Uses MAD. 

    Parameters
    ----------
    layout : 
    df : 
    group_col : string, optional
        The column used for grouping. The default is "Summary".
    threshold : float, optional
        how many MADs to treat as outliers. The default is 3.

    Returns
    -------
    dataframe
        dataframe with the remaining wells.

    """
    
    
    # Copy & unify well names
    df = df.copy()
    df.index = df.index.map(lambda x: unify_well_format(str(x)))
    layout = layout.copy()
    layout['well'] = layout['well'].astype(str).apply(unify_well_format)
    df = df.merge(layout[['well', group_col]],
                  left_index=True, right_on='well').set_index('well')
    
    # Identify time-series columns and do interim smoothing
    numeric_cols = df.select_dtypes(include=[np.number]).columns.tolist()
    df[numeric_cols] = df[numeric_cols].rolling(axis=1, window=3, center=True, min_periods=1).median()

    df['AUC'] = compute_auc(df.drop(columns=['Summary']))['AUC']
    df['log_auc'] = np.log(df['AUC'] + 1e-6)
    df['final_value'] = df[numeric_cols[-1]]
    df['TV'] = df[numeric_cols].diff(axis=1).abs().sum(axis=1)
    
    def detect_mad_outliers(group):
       
        # MAD on log-AUC
        med_auc = np.median(group['log_auc'])
        mad_auc = np.median(np.abs(group['log_auc'] - med_auc))
        out_auc = np.abs(group['log_auc'] - med_auc) > (threshold * mad_auc + 1e-6)
        
        # MAD on final OD
        med_f = np.median(group['final_value'])
        mad_f = np.median(np.abs(group['final_value'] - med_f))
        out_f   = np.abs(group['final_value'] - med_f) > (threshold * mad_f + 1e-6)
        
        # MAD on total variation (spikiness)
        med_tv = np.median(group['TV'])
        mad_tv = np.median(np.abs(group['TV'] - med_tv))
        out_tv = np.abs(group['TV'] - med_tv) > (threshold * mad_tv + 1e-6)
    
        group['outlier'] = out_auc | out_f | out_tv
        
        return group
    
    # Apply per experimental group
    df = df.groupby(group_col, group_keys=False).apply(detect_mad_outliers)
    
    # Filter out and drop helper columns
    kept = df.loc[~df['outlier']]
    return kept.drop(columns=['AUC', 'log_auc', 'final_value', 'outlier', 'TV']), df[['AUC','log_auc','final_value','TV','outlier']]


def subtract_and_impute_blanks(layout, data, n_blanks = 3, fillin_value = 0.01, blank_bool = True):
    """
    Subtracts the N closest blanks to estimate a baseline OD and 
    fills in a baseline OD if negative.

    """
    
    blanked_data = data.copy()
    blank_wells = blanked_data[blanked_data['Summary'] == 'Media control'].index
    blanked_data.iloc[:,:-1] = subtract_closest_blanks(blanked_data.iloc[:,:-1], blank_wells, n_blanks, blank_bool)

    # Subtract blanks
    subtracted = blanked_data.iloc[:,:-1].copy()
    subtracted[subtracted < 0] = fillin_value
    blanked_data.iloc[:,:-1] = subtracted
        
    return blanked_data


def subtract_closest_blanks(df, blank_wells, n, blank_subtraction):
    """ 
    FInds the closest blanks and subtracts them
    """
    
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
        
        # Calculate distances from the current well to each blank well. 
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


def smooth_growth_curves(df, loess_frac):
    """
    Smooths growth curves using LOESS.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame where rows = wells, columns = time points..
    loess_frac : float, optional
        the loess fraction of samples to estimate 
        the smoothing from. The default is 0.2.

    Returns
    -------
    dataframe
        smoothed dataframe.

    """
    numeric_cols = df.select_dtypes(include=["float", "int"]).columns
    df_smoothed = df.copy()
    
    def smooth_row_values(row, loess_frac=0.2):
        
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
    """
    Cleans the dataframe and calculates an AUC.
    """
    # Exclude the last column (experiment descriptor)
    time_series_data = df
    time_points = np.array(time_series_data.columns, dtype=float)  # Convert column names to time points

    auc_results = []
    for well_id, row in time_series_data.iterrows():
        growth_values = row.astype(float).values  # Growth data
        
        # Compute AUC using np.trapz
        auc_trapz = np.trapz(growth_values, time_points)
        auc_results.append([well_id, auc_trapz])

    # Create DataFrame with AUC results
    auc_df = pd.DataFrame(auc_results, columns=["Well", "AUC"])
    auc_df.set_index("Well", inplace=True)
    
    return auc_df


# Plotting scripts

def plot_group_averages_in_hours(df, group_col = "Summary", out_pdf = None):

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
    ax.set_title("Averaged growth curves", fontsize=14, weight="bold")
    ax.set_xlabel("Time (hours)", fontsize=12)
    ax.set_ylabel("OD600", fontsize=12)

    # Legend inside the plot
    ax.legend(title=group_col, loc="best", fontsize=10, title_fontsize=11)

    #ax.set_ylim([0,2.5])
    ax.set_ylim([0, df.drop(columns='Summary').max().max()+0.25])

    # Remove top and right spines to get a clean look
    sns.despine()

    plt.tight_layout()
    
    if out_pdf:
        plt.savefig(out_pdf, bbox_inches="tight")
    
    #plt.show()

def get_summary_color_map(df, palette_name="tab10"):
    """
    Generates a color mapping from the unique values in the "Summary" column.
    """
    groups = sorted(df["Summary"].unique())
    palette = sns.color_palette(palette_name, len(groups))
    return {group: palette[i] for i, group in enumerate(groups)}

def plot_growth_curves(growth_df, color_map=None, annotations_df=None, out_pdf = None):
    """
    Plots a grid of growth curves arranged as a 96-well plate,
    with each subplot’s border colored by its 'Summary' group.
    Legend is placed just below the grid in 3 rows.
    """
    if color_map is None:
        color_map = get_summary_color_map(growth_df, palette_name="tab10")

    rows = list("ABCDEFGH")
    cols = list(range(1, 13))
    time_cols = [c for c in growth_df.columns if c != "Summary"]
    x_vals = [float(tp) / 3600 for tp in time_cols]

    # Precompute global y-axis limit
    global_max = growth_df[time_cols].to_numpy().max() + 0.25

    fig, axes = plt.subplots(8, 12, figsize=(14, 9), sharex=True, sharey=True)

    for i, row in enumerate(rows):
        for j, col in enumerate(cols):
            well = f"{row}{col}"
            ax = axes[i, j]

            # Determine border color
            color = 'lightgray'
            if well in growth_df.index:
                summary = growth_df.loc[well, "Summary"]
                color = color_map.get(summary, 'lightgray')
                # Plot the curve
                ax.plot(x_vals,
                        growth_df.loc[well, time_cols],
                        linestyle='-',
                        color='black',
                        linewidth=1.5)

            # Annotated window shading
            if annotations_df is not None and well in annotations_df.index:
                mu = annotations_df.loc[well, "mu"]
                if mu != 0:
                    start = float(annotations_df.loc[well, "start"])
                    end   = float(annotations_df.loc[well, "end"])
                    ax.axvspan(start, end, color=color, alpha=0.3)

            # Color the spines
            for spine in ax.spines.values():
                spine.set_edgecolor(color)
                spine.set_linewidth(1.5)

            # Clean ticks and limits
            ax.tick_params(labelbottom=False, labelleft=False)
            ax.set_ylim(0, global_max)

            # Row/column labels
            if i == 7:
                ax.set_xlabel(str(col), fontsize=6)
            if j == 0:
                ax.set_ylabel(row, fontsize=6)

    # Only show axes labels on bottom‐left
    bl = axes[7, 0]
    bl.tick_params(labelbottom=True, labelleft=True, labelsize=8)
    bl.set_xlabel(f"{cols[0]}\nTime (h)", fontsize=8)
    bl.set_ylabel(f"{rows[-1]}\nOD600", fontsize=8)

    # Build legend handles
    legend_handles = [
        Line2D([0], [0], color=col, lw=2, label=str(summary))
        for summary, col in sorted(color_map.items())
    ]

    # Determine columns so legend forms 3 rows
    n_items = len(legend_handles)
    ncol = math.ceil(n_items / 3)

    # Place legend just below the grid
    fig.legend(
        handles=legend_handles,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.05),
        ncol=ncol,
        title="Summary",
        frameon=False,
        fontsize=8,
        title_fontsize=9
    )

    # Adjust margins to bring legend closer
    fig.subplots_adjust(bottom=0.1, top=0.98, left=0.02, right=0.98,
                        hspace=0.2, wspace=0.1)
    plt.tight_layout()

    if out_pdf:
        fig.savefig(out_pdf, bbox_inches="tight")


def plot_boxplots(df, color_map=None, out_pdf=None):
    if color_map is None:
        color_map = get_summary_color_map(df, palette_name="tab10")

    fig, axes = plt.subplots(1, 3, figsize=(15, 5))

    # 1) AUC
    sns.boxplot(x="Summary", y="AUC", data=df, ax=axes[0], palette=color_map, hue = "Summary", legend = False)
    axes[0].set(xlabel="", ylabel="AUC")
    axes[0].set_xticklabels([])
    axes[0].tick_params(bottom=False)

    # 2) Growth rate (mu)
    sns.boxplot(x="Summary", y="mu", data=df, ax=axes[1], palette=color_map, hue = "Summary", legend = False)
    axes[1].set(xlabel="", ylabel="mu")
    axes[1].set_xticklabels([])
    axes[1].tick_params(bottom=False)

    # 3) Final OD600
    sns.boxplot(x="Summary", y="FinalOD", data=df, ax=axes[2], palette=color_map, hue = "Summary", legend = False)
    axes[2].set(xlabel="", ylabel="OD600")
    axes[2].set_xticklabels([])
    axes[2].tick_params(bottom=False)

    sns.despine(fig=fig, bottom=True)

    # Legend handles
    legend_handles = [
        Patch(facecolor=c, label=str(g))
        for g, c in sorted(color_map.items())
    ]

    # Place legend just below the plots, in 3 columns → 3 rows for 8 items
    fig.legend(
        handles=legend_handles,
        loc='lower center',
        bbox_to_anchor=(0.5, -0.08),  # slightly below the axes
        ncol=3,
        title="Summary",
        frameon=False,
        fontsize=8,
        title_fontsize=9
    )

    # Make room at the bottom for the legend
    fig.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.18)

    if out_pdf:
        fig.savefig(out_pdf, bbox_inches="tight")

