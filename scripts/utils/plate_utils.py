import pandas as pd
import json
import re
import subprocess
import itertools
import io
import os

def copy_excel_file(file_path, n_copies):
    # Load the original Excel file
    original_excel = pd.ExcelFile(file_path)
    
    # Extract the file name and extension
    base_name = file_path.rsplit('.', 1)[0]
    extension = file_path.rsplit('.', 1)[1]

    # Loop to create copies
    for i in range(2, (n_copies-1) + 2):
        # New file name
        new_file_name = f"{base_name}_{i}.{extension}"
        
        # Save a copy of the original Excel file
        with pd.ExcelWriter(new_file_name) as writer:
            for sheet_name in original_excel.sheet_names:
                sheet_df = original_excel.parse(sheet_name)
                sheet_df.to_excel(writer, sheet_name=sheet_name, index=False)
        
        print(f"Copied: {new_file_name}")

def save_series_to_file(series, file_path):
    # Convert the series to a single string
    series_string = ''.join(map(str, series))
    
    # Format the string with the desired output
    formatted_string = f'channels\n"{series_string}"'
    
    # Save the formatted string to a text file
    with open(file_path, 'w') as file:
        file.write(formatted_string)
    
    print(f"Series saved to {file_path}")

def create_folder_structure(base_path, main_folder, subfolders_structure):
    
    # Define the main folder path
    main_folder_path = os.path.join(base_path, main_folder)
    
    # Create the main folder if it doesn't exist
    if not os.path.exists(main_folder_path):
        os.makedirs(main_folder_path)
        print(f"Created main folder: {main_folder_path}")
    else:
        print(f"Main folder already exists: {main_folder_path}")
    
    # Create subfolders and their subsubfolders
    for subfolder, subsubfolders in subfolders_structure.items():
        # Path for each subfolder
        subfolder_path = os.path.join(main_folder_path, subfolder)
        if not os.path.exists(subfolder_path):
            os.makedirs(subfolder_path)
            print(f"Created subfolder: {subfolder_path}")
        else:
            print(f"Subfolder already exists: {subfolder_path}")

        # Create each subsubfolder within the current subfolder
        for subsubfolder in subsubfolders:
            subsubfolder_path = os.path.join(subfolder_path, subsubfolder)
            if not os.path.exists(subsubfolder_path):
                os.makedirs(subsubfolder_path)
                print(f"Created subsubfolder: {subsubfolder_path}")
            else:
                print(f"Subsubfolder already exists: {subsubfolder_path}")

def load_json(file_path):
    """
    Load a JSON file and return its content.
    """
    try:
        with open(file_path, 'r') as file:
            data = json.load(file)
            return data
    except Exception as e:
        print(f"Error loading JSON file: {e}")
        return None

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

def save_json(data, file_path):
    try:
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)  # indent=4 is used for pretty-printing
    except Exception as e:
        print(f"Error saving JSON file: {e}")

def run_minizinc_command(plate_file):
    
    # Use absolute paths for minizinc and files
    minizinc_path = "/Applications/MiniZincIDE.app/Contents/Resources/minizinc"
    mzn_file = "../../plaid/plate-design.mzn"
    json_file = plate_file
    
    command = f"{minizinc_path} --solver Gecode '{mzn_file}' {json_file}"
    
    try:
        result = subprocess.run(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
        
        if result.returncode == 0:
            return result.stdout
        else:
            print(f"Error: {result.stderr}")
            return None
    except Exception as e:
        print(f"Exception: {e}")
        return None

def output_to_dataframe(output_str):
    
    # Clean up the output to remove extra lines if needed
    output_str_cleaned = output_str.strip()  # Remove leading/trailing whitespaces

    # Use io.StringIO to simulate reading a CSV from a string
    data = io.StringIO(output_str_cleaned)

    # Read the string into a Pandas DataFrame
    df = pd.read_csv(data).iloc[:-1,:]
    
    return df

def plate_filler(plate_layout):

    # Define the list of all well positions (A01 to H12)
    all_wells = [f'{row}{col}' for row, col in itertools.product(['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'], [f'{i:02d}' for i in range(1, 13)])]
    
    # Find missing wells
    existing_wells = set(plate_layout['well'])
    missing_wells = set(all_wells) - existing_wells
    
    # Create placeholder rows for missing wells
    placeholder_rows = pd.DataFrame({'plateID': 'plate_1',
                                     'well': list(missing_wells),
                                     'cmpdname': 'Media',
                                     'CONCuM': 'Media control',
                                     'cmpdnum': 'Media',
                                     'VOLuL': 'nan'})
    
    # Append placeholder rows to the original dataframe
    df_filled = pd.concat([plate_layout, placeholder_rows], ignore_index=True)
    
    # Sort the DataFrame based on the well order (A01 to H12)
    df_filled['well'] = pd.Categorical(df_filled['well'], categories=all_wells, ordered=True)
    df_filled = df_filled.sort_values('well').reset_index(drop=True)

    return df_filled


