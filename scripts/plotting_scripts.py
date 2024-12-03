#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:02:15 2024

@author: danbru
"""
import os
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

gdata = pd.read_csv(os.environ["GEN_EXP_ROOT_DIR"] + '/experiments/amiga_results/data/20241010_hplc_growth_experiment.txt', sep = '\t', index_col = 0)


# Plot wells that were sampled at 10h and those that were sampled at 12h


ten_samples = ['A02','B04','B10','B11','C03','C08','C10','D06','D09','E04','E09','E11','F02','F03','F10','G04','G06']
twelve_samples = ['A03','A04','A05','A06','A08','A09','A10','A11',
                  'B02','B03','B05','B06','B08','B09',
                  'C02','C04','C05','C06','C09','C11',
                  'D02','D03','D04','D05','D08','D10','D11',
                  'E02','E03','E05','E06','E08','E10',
                  'F04','F05','F06','F08','F09','F11',
                  'G02','G03','G05','G08','G09','G10','G11',
                  'H02','H03','H04','H05','H06','H08','H09','H10','H11']
tenhour = gdata.loc[ten_samples]
twelvehour = gdata.loc[twelve_samples]


# Function to plot mean with standard deviation fill
def plot_with_error(ax, df, label, color):
    time_points = df.columns.astype(float)  # Assuming column headers are time points
    mean_values = df.mean(axis=0)  # Mean of each column (time point)
    std_values = df.std(axis=0)    # Std dev of each column (time point)
    
    # Plot the mean values
    ax.plot(time_points, mean_values, label=label, color=color)
    
    # Fill between mean - std and mean + std
    ax.fill_between(time_points, 
                    mean_values - std_values, 
                    mean_values + std_values, 
                    color=color, alpha=0.2)


# Create the plot
fig, ax = plt.subplots(figsize=(8, 6))

# Plot the two dataframes with their mean and standard deviation
plot_with_error(ax, tenhour, label='Sampled at 10h', color='blue')
plot_with_error(ax, twelvehour, label='Sampled at 12h', color='green')

# Customize the plot
ax.set_xlabel('Time (s)')
ax.set_ylabel('OD600')
ax.set_title('OD600 over time (Filtered but uncorrected)')
ax.legend()

# Show the plot
plt.show()