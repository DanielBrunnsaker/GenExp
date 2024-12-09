#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 22 11:02:15 2024

@author: danbru
"""
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

gdata = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_0.25_5_20241121_1256 copy/results/growth/FA_corrected.txt', sep = '\t', index_col = 0)
layout = pd.read_excel('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/glutamate_0.25_5_20241121_1256 copy/protocol/hamilton/pipetting_layout.xlsx')
layout['well'] = layout['well'].apply(lambda x: f"{x[0]}{int(x[1:]):02d}")





# Merge gdata with layout's 'summary' column
#gdata = gdata.T  # Transpose to match well as columns
gdata = gdata.merge(layout[['well', 'Summary']], left_index=True, right_on='well', how='inner')
gdata.set_index('well', inplace=True)


#gdata = gdata[gdata['Summary'].isin(['High glutamate supplementation without formic acid treatment.',
#                            'Low glutamate supplementation without formic acid treatment.',
#                            'Control condition with no supplementation and no formic acid treatment.'])]

# Group data by 'summary'
grouped = gdata.groupby('Summary')


# Function to plot mean with standard deviation fill
def plot_with_error(ax, df, label, color, error_bool):
    time_points = df.columns.astype(float)  # Assuming column headers are time points
    mean_values = df.mean(axis=0)  # Mean of each column (time point)
    std_values = df.std(axis=0)    # Std dev of each column (time point)
    
    # Plot the mean values
    ax.plot(time_points, mean_values, label=label, color=color)
    
    if error_bool == True:
        # Fill between mean - std and mean + std
        ax.fill_between(time_points, 
                        mean_values - std_values, 
                        mean_values + std_values, 
                        color=color, alpha=0.2)



# Create the plot
fig, ax = plt.subplots(figsize=(10, 8))

# Plot each group with unique colors
colors = plt.cm.tab10(np.linspace(0, 1, len(grouped)))  # Generate a list of unique colors
for (summary, group), color in zip(grouped, colors):
    # Drop the 'summary' column to extract only the growth data
    growth_data = group.drop(columns='Summary')
    plot_with_error(ax, growth_data, label=summary, color=color, error_bool = False)

# Customize the plot
ax.set_xlabel('Time (s)')
ax.set_ylabel('OD600')
ax.set_title('OD600 over time by experimental groups')
ax.legend(title='Experimental groups')

# Show the plot
plt.show()



# Create the figure and axes
fig, axes = plt.subplots(1, 2, figsize=(16, 8), sharey=True)

# First plot: Average with error
ax_avg = axes[0]
colors = plt.cm.tab10(np.linspace(0, 1, len(grouped)))  # Generate a list of unique colors
for (summary, group), color in zip(grouped, colors):
    # Drop the 'summary' column to extract only the growth data
    growth_data = group.drop(columns='Summary')
    plot_with_error(ax_avg, growth_data, label=summary, color=color, error_bool=True)

# Customize the average plot
ax_avg.set_xlabel('Time (s)')
ax_avg.set_ylabel('OD600')
ax_avg.set_title('Average OD600 over time by experimental groups')
ax_avg.legend(title='Experimental groups')

# Second plot: Individual traces grouped by summary
ax_individual = axes[1]
for (summary, group), color in zip(grouped, colors):
    # Drop the 'summary' column to extract only the growth data
    growth_data = group.drop(columns='Summary')
    time_points = growth_data.columns.astype(float)  # Assuming column headers are time points
    for row in growth_data.iterrows():  # Iterate through each trace
        ax_individual.plot(
            time_points,
            row[1],
            color=color,
            alpha=0.7,
            label=summary if summary not in ax_individual.get_legend_handles_labels()[1] else "",  # Avoid duplicate labels
        )

# Customize the individual traces plot
ax_individual.set_xlabel('Time (s)')
ax_individual.set_ylabel('OD600')
ax_individual.set_title('Individual OD600 traces by experimental groups')
ax_individual.legend(title='Experimental groups')

# Adjust layout and show the plot
plt.tight_layout()
plt.show()
