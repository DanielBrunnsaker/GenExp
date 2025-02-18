import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.optimize import fsolve
from sklearn.metrics import r2_score

# Step 1: Define True OD and Corresponding Measured OD Readings
true_od = np.array([2.015, 0.999, 0.760, 0.4334, 0.304, .1521])  # Known true OD values (actual concentration)
measured_od = np.array([1.8332, 1.0075, 0.4995, 0.380, 0.2167, 0.1521])  # Spectrophotometer OD readings
blank_od = 0.00  # Blank reference (optional)

# Step 2: Define Various Calibration Models
def linear_model(x, a, b):
    return a * x + b

def exponential_model(x, a, b):
    return a * np.exp(-b * x)

def power_law_model(x, a, b):
    return a * (x ** b)

def michaelis_menten_model(x, a, b):
    return a * x / (b + x)

# Step 3: Fit Each Model to the Data
models = {
    "Linear": linear_model,
    "Exponential": exponential_model,
    "Power Law": power_law_model,
    "Michaelis-Menten": michaelis_menten_model
}

best_model = None
best_r2 = -np.inf
best_params = None
fitted_curves = {}

for name, model in models.items():
    try:
        params, _ = curve_fit(model, true_od, measured_od, maxfev=10000)
        fitted_od = model(true_od, *params)
        r2 = r2_score(measured_od, fitted_od)  # Compute R²
        fitted_curves[name] = (params, fitted_od, r2)
        
        if r2 > best_r2:  # Choose the best model based on R²
            best_r2 = r2
            best_model = model
            best_params = params

        print(f"{name} Model: R² = {r2:.5f}")

    except Exception as e:
        print(f"Could not fit {name} model: {e}")

# Step 4: Plot All Models and the Best Fit
plt.figure(figsize=(8, 5))
plt.scatter(true_od, measured_od, color='black', label="Measured OD", zorder=3)

for name, (params, fitted_od, r2) in fitted_curves.items():
    plt.plot(true_od, fitted_od, linestyle="--", label=f"{name} Fit (R²={r2:.3f})")

plt.axhline(y=blank_od, color='gray', linestyle=':', label="Blank OD")
plt.xlabel("True OD (Actual Concentration)")
plt.ylabel("Measured OD (Spectrophotometer)")
plt.title("OD Calibration - Model Comparison")
plt.legend()
plt.grid(True, linestyle="--", linewidth=0.5)
plt.show()

# Step 5: Define a Function to Estimate True OD from Measured OD
def estimate_true_od(measured_value):
    """ Inverts the best calibration model to estimate true OD from measured OD. """
    if best_model is None:
        raise ValueError("No valid model found for calibration!")
    
    true_od_guess = 1.0  # Initial guess for true OD
    true_od_solution = fsolve(lambda true_od: best_model(true_od, *best_params) - measured_value, true_od_guess)
    return true_od_solution[0]

# Example: Converting a Measured OD to True OD
measured_od_example = 0.3  # Replace with any new OD reading
true_od_estimated = estimate_true_od(measured_od_example)
print(f"\nBest Model: {best_model.__name__} (R² = {best_r2:.5f})")
print(f"Measured OD: {measured_od_example}, Estimated True OD: {true_od_estimated:.4f}")











import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit, fsolve
from sklearn.metrics import r2_score


def get_calibration_function(true_od, measured_od):
    """
    Fits multiple calibration models to OD data and selects the best model.

    Parameters:
        true_od (array-like): Known true OD values (actual concentration).
        measured_od (array-like): Corresponding OD readings from the spectrophotometer.

    Returns:
        function: A function that converts measured OD to true OD.
    """

    # Step 1: Define Calibration Models
    def linear_model(x, a, b): return a * x + b
    def exponential_model(x, a, b): return a * np.exp(-b * x)
    def power_law_model(x, a, b): return a * (x ** b)
    def michaelis_menten_model(x, a, b): return a * x / (b + x)

    models = {
        "Linear": linear_model,
        "Exponential": exponential_model,
        "Power Law": power_law_model,
        "Michaelis-Menten": michaelis_menten_model
    }

    # Step 2: Fit Each Model and Select the Best One
    best_model = None
    best_r2 = -np.inf
    best_params = None

    for name, model in models.items():
        try:
            params, _ = curve_fit(model, true_od, measured_od, maxfev=10000)
            fitted_od = model(true_od, *params)
            r2 = r2_score(measured_od, fitted_od)

            if r2 > best_r2:
                best_r2 = r2
                best_model = model
                best_params = params

        except Exception as e:
            print(f"Could not fit {name} model: {e}")

    if best_model is None:
        raise ValueError("No valid model found for calibration!")

    print(f"Best Calibration Model: {best_model.__name__} (R² = {best_r2:.5f})")

    # Step 3: Return a Function to Convert Measured OD to True OD
    def calibrated_function(measured_value):
        """ Applies the best calibration model to estimate true OD. """
        true_od_guess = 1.0
        true_od_solution = fsolve(lambda true_od: best_model(true_od, *best_params) - measured_value, true_od_guess)
        return true_od_solution[0]

    return np.vectorize(calibrated_function)  # Vectorized for efficient DataFrame operations


# Get the calibration function
calibration_func = get_calibration_function(true_od, measured_od)

# Example DataFrame
df = pd.DataFrame({
    "Sample_ID": ["A", "B", "C", "D", "E", "F"],
    "Measured_OD": [0.7, 0.35, 0.18, 0.1, 0.05, 0.02]  # New measured OD values
})

# Apply the calibration function to the DataFrame
df["Calibrated_OD"] = df["Measured_OD"].apply(calibration_func)

# Display the updated DataFrame
import ace_tools as tools
tools.display_dataframe_to_user(name="Calibrated OD Data", dataframe=df)


# Example Usage
#true_od_data = np.array([1, 2, 4, 8, 16, 32])  # Known true OD values
#measured_od_data = np.array([0.8, 0.42, 0.22, 0.12, 0.06, 0.03])  # Spectrophotometer readings

# Get the best model function and transformation function
best_model_func, transform_dataframe = calibrate_od(true_od, measured_od)

# Example DataFrame
df = pd.DataFrame({
    "Sample_ID": ["A", "B", "C", "D", "E", "F"],
    "Measured_OD": [0.7, 0.35, 0.18, 0.1, 0.05, 0.02]  # New measured OD values
})

# Apply the calibration to the DataFrame
df_calibrated = transform_dataframe(df)







# Display the updated DataFrame
import ace_tools as tools
tools.display_dataframe_to_user(name="Calibrated OD Data", dataframe=df_calibrated)


# Example Usage
#true_od_data = np.array([1, 2, 4, 8, 16, 32])  # Known true OD values
#measured_od_data = np.array([0.8, 0.42, 0.22, 0.12, 0.06, 0.03])  # Spectrophotometer readings

# Get the best model function and transformation function
best_model_func, transform_curve = calibrate_od(true_od, measured_od)

# Example: Convert a single OD measurement
measured_example = 0.15
true_example = best_model_func(measured_example)
print(f"Measured OD: {measured_example}, Estimated True OD: {true_example:.4f}")

# Example: Convert an entire dataset
new_measured_curve = np.array([0.7, 0.35, 0.18, 0.1, 0.05, 0.02])
calibrated_curve = transform_curve(new_measured_curve)
print(f"Measured Curve: {new_measured_curve}")
print(f"Calibrated Curve: {calibrated_curve}")


