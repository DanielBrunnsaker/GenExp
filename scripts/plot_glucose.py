#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 29 09:51:42 2024

@author: danbru
"""

import os
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize, curve_fit
import pandas as pd

# Example data - replace with actual measurements
time_data = np.array([0, 10, 12])  # Time in hours
glucose_data = np.array([(19.760+18.392+20.170)/3, (13.909+13.921+14.140)/3, (10.921+9.496+11.047)/3])  # Glucose concentration (arbitrary units)

biomass_data = pd.read_csv(os.environ["GEN_EXP_ROOT_DIR"] + '/experiments/amiga_results/data/20241002_24_growth_experiment.txt', sep = '\t', index_col = 0)
biomass_data = biomass_data[['0','36000','43200']]
biomass_data = biomass_data.loc[['A03','A11','D04','F09','H09']].mean(axis = 0).to_numpy()#*0.34

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint
from scipy.optimize import minimize, curve_fit

# Example data - replace with actual measurements
time_data = np.array([0, 10, 12])  # Time in hours
#glucose_data = np.array([10.0, 5.0, 3.0])  # Glucose concentration (arbitrary units)
#biomass_data = np.array([0.1, 0.8, 1.2])  # Biomass concentration (arbitrary units)
glucose_data = np.array([(19.760+18.392+20.170)/3, (13.909+13.921+14.140)/3, (10.921+9.496+11.047)/3])  # Glucose concentration (arbitrary units)

# Extend the time data for growth measurements
time_data_growth = np.array([0, 10, 12, 14, 16, 18, 20])  # Time in hours
#biomass_data_growth = np.array([0.1, 0.8, 1.2, 1.6, 2.0, 2.4])  # Extended biomass concentrations
biomass_data = pd.read_csv(os.environ["GEN_EXP_ROOT_DIR"] + '/experiments/amiga_results/data/20241002_24_growth_experiment.txt', sep = '\t', index_col = 0)
biomass_data = biomass_data[['0','36000','43200','50400','57600','64800','72000']]
biomass_data_growth = biomass_data.loc[['A02','A03','A04','E03','D03','F03']].mean(axis = 0).to_numpy()*0.34

# Define Monod growth model for biomass
def monod_growth(X, t, mu_max, Ks, Y):
    G = X[1]  # Glucose concentration
    dXdt = mu_max * (G / (Ks + G)) * X[0]
    dGdt = -dXdt / Y  # Glucose consumption proportional to growth
    return [dXdt, dGdt]

# Function to simulate biomass and glucose over time
def simulate_biomass_glucose(params, time_points):
    mu_max, Ks, Y = params
    X0 = [biomass_data_growth[0], glucose_data[0]]  # Initial biomass and glucose concentrations
    solution = odeint(monod_growth, X0, time_points, args=(mu_max, Ks, Y))
    return solution[:, 0], solution[:, 1]  # Biomass and Glucose

# Define a loss function to fit biomass data
def loss(params):
    biomass_simulated, _ = simulate_biomass_glucose(params, time_data_growth)
    return np.sum((biomass_simulated - biomass_data_growth) ** 2)

# Initial guess for mu_max, Ks, and Y
initial_guess = [0.1, 0.1, 0.5]
result = minimize(loss, initial_guess, method='Powell', options={'maxfev': 2000})

# Extract fitted parameters
mu_max, Ks, Y = result.x
print("Fitted parameters for biomass (Monod model):", mu_max, Ks, Y)

# Predict biomass and glucose over a longer time period
time_future = np.linspace(0, 20, 100)
biomass_pred, glucose_pred = simulate_biomass_glucose(result.x, time_future)

# Plotting the results with glucose prediction
plt.figure(figsize=(8, 6))

# Plot biomass
plt.plot(time_future, biomass_pred, 'r-', label='Predicted Biomass', alpha=0.7)
plt.scatter(time_data_growth, biomass_data_growth, color='red', marker='o', label='Observed Biomass')
plt.xlabel('Time (h)')
plt.ylabel('Biomass Concentration (g/L)', color='r')
plt.tick_params(axis='y', labelcolor='r')

# Create a second y-axis for glucose
ax2 = plt.gca().twinx()
ax2.plot(time_future, glucose_pred, 'b-', label='Predicted Glucose', alpha=0.7)
ax2.scatter(time_data, glucose_data, color='blue', marker='o', label='Observed Glucose')
ax2.set_ylabel('Glucose Concentration (g/L)', color='b')
ax2.tick_params(axis='y', labelcolor='b')

# Combine legends
#lines1, labels1 = plt.gca().get_legend_handles_labels()  # Main axis
#lines2, labels2 = ax2.get_legend_handles_labels()         # Secondary axis

# Create a single legend with unique entries
#unique_lines = lines1 + [line for line in lines2 if line not in lines1]
#unique_labels = labels1 + [label for label in labels2 if label not in labels1]

#plt.legend(unique_lines, unique_labels, loc='upper left')

plt.title('Biomass Growth and Predicted Glucose Concentration Over Time')
plt.grid()
plt.tight_layout()
plt.show()

# Print fitted parameters for reference
print("Fitted parameters for biomass (Monod model):", mu_max, Ks, Y)
