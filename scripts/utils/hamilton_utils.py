#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 13:30:19 2025

@author: danbru
"""
import re
import pubchempy as pcp
import pandas as pd

def percentage_to_volume(final_concentration_percent, total_volume, ureg):
    """
    Convert percentage concentration to volume of solute.
    """
    return (final_concentration_percent / 100) * total_volume.to('microliter').magnitude  # Return as float



def convert_to_molar(concentration, molecular_weight):
    """
    Converts a concentration to molar units (mol/L) using the molecular weight.
    Handles units like mg/mL, mM, etc.
    """
    if concentration.check('[mass] / [volume]'):
        # Convert mass concentration to molar concentration
        molar_conc = (concentration / molecular_weight).to('mM')
        return molar_conc
    elif concentration.check('[substance] / [volume]'):
        # Already in molar units
        return concentration.to('mM')
    else:
        raise ValueError("Unsupported concentration units for molar conversion.")

def calculate_percentage_volume(percentage, final_volume, percentage_type, stock_percentage, Q_):
    
    
    """
    Calculates the volume of stock solution needed based on percentage.

    Args:
        percentage (float): Desired final percentage (e.g., 10 for 10%)
        final_volume (Quantity): Final volume in microliters
        percentage_type (str): 'v/v' or 'w/v'
        stock_percentage (float): Percentage concentration of the stock solution

    Returns:
        Volume to add in microliters
    """
    if percentage_type == 'v/v':
        # Volume of solute = (desired_percentage / stock_percentage) * final_volume
        V_solute = (percentage / stock_percentage) * final_volume
        return V_solute.to('microliter').magnitude
    elif percentage_type == 'w/v':
        # Convert final_volume to milliliters
        final_volume_ml = final_volume.to('milliliter').magnitude
        # Mass of solute required in grams
        mass_solute = (percentage / 100.0) * final_volume_ml  # grams
        # Stock concentration in g/mL
        stock_concentration = (stock_percentage / 100.0)  # g/mL
        # Volume of stock solution needed in milliliters
        V_stock_ml = mass_solute / stock_concentration  # mL
        # Convert volume to microliters
        V_stock_ul = Q_(V_stock_ml, 'milliliter').to('microliter').magnitude
        return V_stock_ul
    else:
        raise ValueError("percentage_type must be 'v/v' or 'w/v'.")
        

def parse_concentration(concentration_str, Q_):
    """
    Parses a concentration string and returns a tuple:
    (type, Quantity, percentage_type)
    where type is 'percentage' or 'molar/mass'.
    """
    # Ensure the input is a string and remove any annotations
    if concentration_str is None:
        concentration_str = ''
    else:
        concentration_str = re.sub(r'\s*\(.*?\)', '', str(concentration_str)).strip()

    if '°' in concentration_str:
        concentration_str = '0 mM'
    else:
        concentration_str = re.sub(r'\s*\(.*?\)', '', str(concentration_str)).strip()

    # Handle 'None', 'none', 'null', '0', '0.0', or empty strings
    if concentration_str.lower() in ['none', 'null', '', '0', '0.0']:
        return ('molar/mass', Q_(0, 'dimensionless'), None)

    if '%' in concentration_str:
        # Determine if it's v/v or w/v based on annotations or context
        if 'v/v' in concentration_str.lower():
            percentage_type = 'v/v'
        elif 'w/v' in concentration_str.lower():
            percentage_type = 'w/v'
        else:
            # Default to w/v if not specified
            percentage_type = 'w/v'
        # Extract the numeric value
        matches = re.findall(r"[\d.]+", concentration_str)
        if matches:
            percentage_value = float(matches[0])
        else:
            raise ValueError(f"Unable to parse percentage concentration: '{concentration_str}'.")
        return ('percentage', Q_(percentage_value, '%'), percentage_type)
    else:
        try:
            quantity = Q_(concentration_str)
            return ('molar/mass', quantity, None)
        except Exception as e:
            raise ValueError(f"Unable to parse concentration string: '{concentration_str}'. Error: {e}")

def calculate_volume(C_stock, C_final, V_total):
    """
    Calculates the volume of stock solution needed to achieve the final concentration.
    """
    V1 = (C_final * V_total) / C_stock
    return V1.to('microliter').magnitude

def compute_volume_for_compound(compound_name,final_conc, final_conc_type,final_percentage_type,stock_conc,
    stock_conc_type,stock_percentage,total_volume,Q_,ureg,print_name_for_error=''):
    """
    Computes the volume (in uL) needed from a stock solution to achieve
    a target final concentration in a total reaction volume.

    Returns a float (uL) of how much to add. If final_conc is zero
    or if the compound is not found (and can't be computed),
    returns 0.0 by default.
    """

    # If final concentration is zero, no volume needed
    if final_conc.magnitude <= 0:
        return 0.0

    # If final concentration is a percentage
    if final_conc_type == 'percentage':
        if stock_conc_type != 'percentage':
            raise ValueError(
                f"Stock concentration and final concentration units do not match for {print_name_for_error or compound_name}."
            )

        volume_uL = calculate_percentage_volume(
            final_conc.magnitude,
            total_volume,
            final_percentage_type,
            stock_percentage,
            Q_
        )
        return volume_uL

    # Otherwise, handle molar concentration or similar
    # Attempt PubChem lookup if we need a molecular weight
    import pubchempy as pcp
    compounds = pcp.get_compounds(compound_name, 'name')
    if not compounds:
        print(f"Error: Compound '{compound_name}' not found in PubChem.")
        return 0.0

    molecular_weight = float(compounds[0].molecular_weight) * ureg('g/mol')

    # Convert final concentration from e.g. mg/mL or M to just M for volume calculation
    final_molar = convert_to_molar(final_conc, molecular_weight)
    stock_molar = convert_to_molar(stock_conc, molecular_weight)

    volume_uL = calculate_volume(stock_molar, final_molar, total_volume)
    return volume_uL



def extract_highest_doses(json_data, ureg, Q_):
    
    # Function to convert dose to numeric value
    #def parse_dose(dose):
    #    if dose is None:
    #        return 0 * ureg.mM  # Default to 0 mM
    #    
    #    if '%' in dose:
    #        dose = dose.split('%',' ')[0]
    #    
    #    value, unit = dose.split()
    #    return float(value) * ureg(unit)

    #def parse_dose(dose):
    #    if dose is None:
    #        return 0 * ureg('mM')  # or whatever your default is
   # 
   #     # look for percent with w/v
   #     m = re.match(r'([\d\.]+)\s*%\s*\(w/v\)', dose, flags=re.IGNORECASE)
   #     if m:
   #         value = float(m.group(1))
   #         # grams per 100 milliliters
   #         return value * ureg('g') / (ureg('L'))
   # 
   #     # fallback for other units
   #     value, unit = dose.split(None, 1)
   #     return float(value) * ureg(unit)

    def parse_dose(dose):
        if dose is None:
            return 0 * ureg.mM  # Default to 0 mM
    
        if '%' in dose:
            # percent (w/v) → grams per 100 mL
            value = float(dose.split('%')[0])
            return value * ureg.g / (100 * ureg.mL)
    
        # split into exactly two parts: value and the rest as unit
        value_str, unit_str = dose.split(None, 1)
        return float(value_str) * ureg(unit_str)


    # Extracting highest doses
    highest_doses = {
        "supplement": (None, 0 * ureg.mM),
        "treatment": (None, 0 * ureg.mM),
        "negative_control": (None, 0 * ureg.mM)
    }

    max_supplement = 0 * ureg.mM
    max_negative_control = 0 * ureg.mM
    max_supplement_name = None
    max_negative_control_name = None

    for exp in json_data['experiments']:
        if "negative" in exp['type'].lower():
            category = "negative_control"
        elif "experiment" in exp['type'].lower():
            category = "treatment"
        else:
            category = "supplement"
        
        if exp['media_supplementation']:
            dose_value = parse_dose(exp['media_supplementation_doses'])
            
            try:
                dose_value = convert_to_molar(dose_value, Q_(float(pcp.get_compounds(exp['media_supplementation'], 'name')[0].exact_mass), 'g/mol'))
                dose_value = convert_to_molar(dose_value, Q_(float(pcp.get_compounds(exp['media_supplementation'], 'name')[0].exact_mass), 'g/mol'))
            except:
                molweight = float(input(f"  -Please provide molecular weight of {exp['media_supplementation']} in g/mol."))
                dose_value = convert_to_molar(dose_value, Q_(molweight, 'g/mol'))
                
            if category == "negative_control":
                if dose_value > max_negative_control:
                    max_negative_control = dose_value
                    max_negative_control_name = exp['media_supplementation']
            else:
                if dose_value > max_supplement:
                    max_supplement = dose_value
                    max_supplement_name = exp['media_supplementation']
        
        if exp['treatment'] and category == "treatment":
            dose_value = parse_dose(exp['treatment_parameters'])
            
            try:
                dose_value = convert_to_molar(dose_value, Q_(float(pcp.get_compounds(exp['treatment'].replace(' derivative',''), 'name')[0].exact_mass), 'g/mol'))
            except:
                molweight = float(input(f"  -Please provide molecular weight of {exp['treatment'].replace(' derivative','')} in g/mol."))
                dose_value = convert_to_molar(dose_value, Q_(molweight, 'g/mol'))
                
            if dose_value > highest_doses["treatment"][1]:
                highest_doses["treatment"] = (exp['treatment'], dose_value)

    highest_doses["supplement"] = (max_supplement_name, max_supplement)
    highest_doses["negative_control"] = (max_negative_control_name, max_negative_control)
    
    return highest_doses



def find_rows_by_inchikey(file_path, inchikey):
    """
    Find all rows in an Excel file where the given compound's InChIKey is present.

    :param file_path: Path to the Excel file.
    :param compound_name: Compound name to search for.
    :return: DataFrame of matching rows (or None if no matches found).
    """
    try:

        # Load the Excel file
        df = pd.read_excel(file_path, dtype=str)  # Read all data as strings
        matching_rows = df[df.apply(lambda row: row.astype(str).str.contains(inchikey, case=False, na=False).any(), axis=1)]

        if not matching_rows.empty:
            return matching_rows
        else:
            return None

    except Exception as e:
        print(f"Error: {e}")
        return None

def get_closest_concentration(df, X, ureg):
    
    X = X*5 # Hardcoding this for now, just means that there is a limit to the concentration of the stocks, as we need to adapt to the concentration of the YNB-media
    df = df.copy()
    
    # Convert the "Available concentrations" column to Pint quantities
    df["Available concentrations"] = df["Available concentrations"].apply(lambda x: ureg(x))
    df["Available concentrations"] = df["Available concentrations"].apply(lambda x: x.to(X.units))
    df["Concentration Magnitude"] = df["Available concentrations"].apply(lambda x: x.magnitude)
    df_filtered = df[df["Concentration Magnitude"] >= X.magnitude]

    # If no concentrations are greater or equal, return the smallest available
    if df_filtered.empty:
        return df.loc[df["Concentration Magnitude"].idxmin()]

    # Return the row with the minimum concentration that is still greater than or equal to X
    return df_filtered.loc[df_filtered["Concentration Magnitude"].idxmin()]
