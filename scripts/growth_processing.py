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

def main():
    
    
    '''
    
    default values
    python growth_processing.py --path "/Volumes/EVE/20241010 HPLC" --output /Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/growth/20241010_hplc_growth_experiment.txt
    python growth_processing.py --path "/Volumes/EVE/20241010 HPLC" --output "../experiments/growth/20241010_hplc_growth_experiment.txt"
    python growth_processing.py --path "/Volumes/EVE/20241010 HPLC" --output "/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/amiga_results/data/20241010_hplc_growth_experiment.txt " --blank distance
    
    
    
    /Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/experiments/amiga_results/data/20241010_hplc_growth_experiment.txt 
    
    Then run:
    
    python "../amiga-master/amiga.py" summarize -i ../experiments/amiga_results/data/20241010_hplc_growth_experiment.txt
    
    '''
    
    parser = argparse.ArgumentParser(description='Script with a command-line argument.')
    parser.add_argument('--path', type=str, help='Value for the "path" variable.')
    parser.add_argument('--output', type=str, help='Value for the "output" variable.')
    parser.add_argument('--blank', type=str, help='Value for the "blank" variable.')
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
        
    if args.blank is not None:
        blank = args.blank
        #print(f'Variable "target" set to: {target}')
    else:
        print('Variable "blank" not provided.')
    
    if args.n is not None:
        n = args.n
        #print(f'Variable "target" set to: {target}')
    else:
        print('Variable "n" not provided.')
 

    
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
    
    
    if blank == 'distance':
    
        
        # This needs to be loaded in somehow when I generate the experimental design. Should not be too hard
        blank_wells = ['A01','B01','C01','D01','E01','F01','G01','H01', 
                       'A07','B07','C07','D07','E07','F07','G07','H07',
                       'A01','B12','C12','D12','E12','F12','G12','H12',]  # Replace with your list of blank wells
        
        # Apply the normalization function
        data = subtract_closest_blanks(data, blank_wells, n)
        #print(data.iloc[0,:])
        data.to_csv(output, sep='\t')
    else:
        data.to_csv(output, sep='\t')
        #print(data.iloc[0,:])

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