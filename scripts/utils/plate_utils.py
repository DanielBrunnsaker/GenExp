import pandas as pd
import json
import subprocess
import itertools
import io
import os
import config

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
    
def create_folder_structure(base_path, main_folder, subfolders_structure):
    
    # Define the main folder path
    main_folder_path = os.path.join(base_path, main_folder)
    
    # Create the main folder if it doesn't exist
    if not os.path.exists(main_folder_path):
        os.makedirs(main_folder_path)
    
    # Create subfolders and their subsubfolders
    for subfolder, subsubfolders in subfolders_structure.items():
        # Path for each subfolder
        subfolder_path = os.path.join(main_folder_path, subfolder)
        if not os.path.exists(subfolder_path):
            os.makedirs(subfolder_path)

        # Create each subsubfolder within the current subfolder
        for subsubfolder in subsubfolders:
            subsubfolder_path = os.path.join(subfolder_path, subsubfolder)
            if not os.path.exists(subsubfolder_path):
                os.makedirs(subsubfolder_path)

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


def save_json(data, file_path):
    try:
        with open(file_path, 'w') as file:
            json.dump(data, file, indent=4)  # indent=4 is used for pretty-printing
    except Exception as e:
        print(f"Error saving JSON file: {e}")

def run_minizinc_command(plate_file):
    
    # Use absolute paths for minizinc and files
    minizinc_path = config.MINIZINC_PATH
    mzn_file = "../plaid/plate-design.mzn"    
    command = f"{minizinc_path} --solver Gecode '{mzn_file}' {plate_file}"
    
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


