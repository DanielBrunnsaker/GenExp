#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 10 11:24:20 2024

@author: danbru
"""





import pandas as pd
import matplotlib.pyplot as plt
import math
import seaborn as sns
import numpy as np

df = pd.read_excel('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/AA.xls', sheet_name = 'intracellular_concentration_mM', index_col = 0).iloc[:,1:]

# Assume 'df' is your DataFrame
# Replace 'df' with the actual variable name of your DataFrame

# Extract the values of the row labeled 'YOR202W'
highlight_values = df.loc["YDL227C"]

# Number of columns
num_columns = len(df.columns)
num_rows = math.ceil(num_columns / 4)  # Adjust number of rows to fit 4 columns per row

# Create subplots
fig, axes = plt.subplots(num_rows, 4, figsize=(20, 5 * num_rows))
axes = axes.flatten()

for idx, column in enumerate(df.columns):
    ax = axes[idx]
    
    # Plot histogram
    sns.histplot(df[column], kde=True, bins=30, ax=ax, color='blue', alpha=0.7)
    
    # Highlight the specific value in the histogram
    if column in highlight_values:
        ax.axvline(highlight_values[column], color='red', linestyle='--', linewidth=2)
    
    ax.set_title(f"{column}")
    ax.set_xlabel(column)
    ax.set_ylabel("Frequency")
    ax.legend()

# Hide any extra subplots if the number of columns is not a perfect multiple of 4
for idx in range(num_columns, len(axes)):
    fig.delaxes(axes[idx])

plt.tight_layout()
plt.show()




# Compute the mean and standard deviation for each column
column_means = df.mean()
column_std = df.std()

# Extract the row values
yor202w_values = df.loc["YDL227C"]

# Compute z-scores for the row values relative to column statistics
z_scores = (yor202w_values - column_means) / column_std

# Define significance threshold (e.g., 2 or 3)
threshold = 2  # Adjust as needed

# Check if the absolute z-score exceeds the threshold
significant_differences = z_scores.abs() > threshold

# Display results
results = pd.DataFrame({
    "Column": df.columns,
    "YOR202W Value": yor202w_values,
    "Column Mean": column_means,
    "Column Std Dev": column_std,
    "Z-Score": z_scores,
    "Significantly Different": significant_differences
}).reset_index(drop=True)



# Calculate the mean for each column
column_means = df.mean()

# Calculate the Euclidean distance of each row to the column means
distances = df.apply(lambda row: np.linalg.norm(row - column_means), axis=1)

# Find the row with the smallest distance
closest_row_index = distances.idxmin()
closest_row = df.loc[closest_row_index]

# Display the result
print(f"The row closest to the mean is: {closest_row_index}")
print(closest_row)




# Assume `df` is your DataFrame and `row_names` is the list of row names to check
row_names = ["YEL024W","YDR497C","YGL256W","YPL058C","YBL039C","YBR084W","YDR035W","YNL030W","YNL117W","YJL070C","YDL227C"]  # Replace with your actual list of row names

# Standardize the DataFrame (z-score normalization)
standardized_df = (df - df.mean()) / df.std()

# Calculate the standardized mean for each column (will be 0, but included for clarity)
standardized_mean = standardized_df.mean()

# Filter the standardized DataFrame to only include the specified rows
filtered_standardized_df = standardized_df.loc[row_names]

# Calculate the Euclidean distance of each specified row to the standardized mean
distances = filtered_standardized_df.apply(lambda row: np.linalg.norm(row - standardized_mean), axis=1)

# Find the row with the smallest distance
most_representative_row = distances.idxmin()
closest_row = df.loc[most_representative_row]

# Display the result
print(f"The row most representative of the mean is: {most_representative_row}")
print(closest_row)

# Optionally, display the distance values for the specified rows
results = pd.DataFrame({
    "Row Name": filtered_standardized_df.index,
    "Distance to Mean": distances
}).sort_values(by="Distance to Mean")

import ace_tools as tools; tools.display_dataframe_to_user(name="Most Representative Row from Specified List", dataframe=results)


