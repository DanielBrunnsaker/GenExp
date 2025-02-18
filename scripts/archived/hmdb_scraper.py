import requests
from bs4 import BeautifulSoup
import pandas as pd
import time

# Crawl-delay specified in robots.txt

# Define common adduct types and their contributions


def fetch_metabolite_page(hmdb_id):
    """
    Constructs the URL for a metabolite using its HMDB ID.
    """
    url = f"https://hmdb.ca/metabolites/{hmdb_id}"
    headers = {"User-Agent": "Mozilla/5.0"}
    response = requests.get(url, headers=headers)
    time.sleep(CRAWL_DELAY)  # Respect crawl-delay after fetching the page

    if response.status_code != 200:
        print(f"Failed to fetch page for HMDB ID {hmdb_id}. Status code: {response.status_code}")
        return None

    return BeautifulSoup(response.text, "html.parser")  # Return BeautifulSoup object

def fetch_common_name(soup):
    """
    Extracts the Common Name of the compound from the metabolite page.
    """
    header = soup.find("th", string="Common Name")
    if header:
        value_cell = header.find_next("td")
        if value_cell:
            return value_cell.text.strip()
    print("Common Name not found.")
    return None

def fetch_monoisotopic_mass(soup):
    """
    Extracts the Monoisotopic Molecular Weight from the metabolite page.
    """
    header = soup.find("th", string="Monoisotopic Molecular Weight")
    if header:
        value_cell = header.find_next("td")
        if value_cell:
            return float(value_cell.text.strip())
    print("Monoisotopic Molecular Weight not found.")
    return None

def fetch_monoisotopic_mass(soup):
    """
    Extracts the Monoisotopic Molecular Weight from the metabolite page.
    If the value is 'Not Available' or cannot be converted to float, return None.
    """
    header = soup.find("th", string="Monoisotopic Molecular Weight")
    if header:
        value_cell = header.find_next("td")
        if value_cell:
            value = value_cell.text.strip()
            try:
                return float(value)
            except ValueError:
                print(f"Monoisotopic Molecular Weight is not a valid number: '{value}'")
                return None  # Handle cases where value is "Not Available" or another invalid format
    print("Monoisotopic Molecular Weight not found.")
    return None

def parse_ccs_table(table, common_name, ccs_type, monoisotopic_mass):
    """
    Parses a CCS table and extracts relevant data, including m/z.
    """
    rows = []
    tbody = table.find("tbody")
    if not tbody:
        print(f"No table body found in {ccs_type} section for {common_name}.")
        return []

    for tr in tbody.find_all("tr"):  # Iterate through all rows
        cells = tr.find_all("td")
        
        # Handle Experimental and Predicted table structures
        if ccs_type == "Experimental" and len(cells) >= 4:
            adduct_type = cells[0].text.strip()
            source_or_predictor = cells[1].text.strip()
            ccs_value = cells[2].text.strip()
            reference = cells[3].text.strip()
        elif ccs_type == "Predicted" and len(cells) >= 4:
            source_or_predictor = cells[0].text.strip()  # Predictor
            adduct_type = cells[1].text.strip()  # Adduct Type
            ccs_value = cells[2].text.strip()
            reference = cells[3].text.strip()
        else:
            continue  # Skip rows with insufficient data

        # Only calculate m/z if CCS Value is present
        mz = None
        if ccs_value and adduct_type in adducts and monoisotopic_mass is not None:
            adduct = adducts.get(adduct_type)
            mz = (monoisotopic_mass + adduct["mass"]) / abs(adduct["charge"])  # Fix: Use absolute charge

        rows.append({
            "Name": common_name,
            "Monoisotopic Mass": monoisotopic_mass,
            "Adduct Type": adduct_type,
            "CCS Type": ccs_type,
            "Data Source / Predictor": source_or_predictor,
            "m/z": mz,
            "CCS Value (Å²)": ccs_value,
            "Reference": reference,
        })
    return rows

def fetch_ccs_data(hmdb_id):
    """
    Fetches CCS data (experimental and predicted) and additional details from a metabolite page using HMDB ID.
    """
    soup = fetch_metabolite_page(hmdb_id)
    if not soup:
        return []

    # Extract Common Name
    common_name = fetch_common_name(soup)
    
    # Extract Monoisotopic Molecular Weight
    monoisotopic_mass = fetch_monoisotopic_mass(soup)

    # Extract InChI Key
    inchi_key = fetch_inchi_key(soup)

    # Extract Chemical Formula
    chemical_formula = fetch_chemical_formula(soup)
    
    # Extract CCS Data
    data = []

    # Parse Experimental Collision Cross Sections Table
    #experimental_table = soup.find("table", {"class": "table table-bordered ccs"})
    #if experimental_table:
    #    print(f"Parsing Experimental table for HMDB ID {hmdb_id} ({common_name})...")
    #    data += parse_ccs_table(experimental_table, common_name, "Experimental", monoisotopic_mass)

    # Parse Predicted Collision Cross Sections Table
    predicted_section = soup.find("th", string="Predicted Chromatographic Properties")
    if predicted_section:
        predicted_table = predicted_section.find_next("table", {"class": "table table-bordered ccs"})
        if predicted_table:
            print(f"Parsing Predicted table for HMDB ID {hmdb_id} ({common_name})...")
            data += parse_ccs_table(predicted_table, common_name, "Predicted", monoisotopic_mass)

    # Add InChI Key and Chemical Formula to all rows
    for row in data:
        row["InChI Key"] = inchi_key
        row["Chemical Formula"] = chemical_formula

    return data

def fetch_chemical_formula(soup):
    """
    Extracts the Chemical Formula of the compound from the metabolite page.
    """
    header = soup.find("th", string="Chemical Formula")
    if header:
        value_cell = header.find_next("td")
        if value_cell:
            # Extract the full text content, including subscript formatting
            formula = ''.join(value_cell.stripped_strings)
            return formula.strip()
    print("Chemical Formula not found.")
    return None

def fetch_inchi_key(soup):
    """
    Extracts the InChI Key of the compound from the metabolite page.
    """
    header = soup.find("th", string="InChI Key")
    if header:
        value_cell = header.find_next("td")
        if value_cell:
            return value_cell.text.strip()
    print("InChI Key not found.")
    return None

def scrape_metabolites(hmdb_ids):
    """
    Scrapes a list of HMDB IDs and returns their CCS data as a pandas DataFrame.
    """
    df = pd.DataFrame(columns=[
        "Name", "Monoisotopic Mass", "InChI Key", "Chemical Formula", 
        "CCS Type", "Adduct Type", "Data Source / Predictor", "m/z", 
        "CCS Value (Å²)", "Reference"
    ])
    for hmdb_id in hmdb_ids:
        print(f"Scraping HMDB ID {hmdb_id}...")
        data = fetch_ccs_data(hmdb_id)
        if data:
            df = pd.concat([df, pd.DataFrame(data)], ignore_index=True)
        else:
            print(f"No data found for HMDB ID {hmdb_id}.")
        time.sleep(CRAWL_DELAY)  # Respect crawl-delay between queries
    return df

def scrape_metabolites(hmdb_ids, output_csv="metabolites_data.csv"):
    """
    Scrapes a list of HMDB IDs and saves their CCS data to a CSV file incrementally.
    """
    columns = [
        "Name", "Monoisotopic Mass", "InChI Key", "Chemical Formula", 
        "CCS Type", "Adduct Type", "Data Source / Predictor", "m/z", 
        "CCS Value (Å²)", "Reference"
    ]
    
    # Create and initialize the CSV file
    df = pd.DataFrame(columns=columns)
    df.to_csv(output_csv, index=False, mode='w')
    
    for hmdb_id in hmdb_ids:
        print(f"Scraping HMDB ID {hmdb_id}...")
        data = fetch_ccs_data(hmdb_id)
        
        if data:
            new_df = pd.DataFrame(data)
            new_df.to_csv(output_csv, index=False, mode='a', header=False)
        else:
            print(f"No data found for HMDB ID {hmdb_id}.")
        
        time.sleep(5)  # Respect crawl-delay between queries
        
    return df

CRAWL_DELAY = 5  # seconds
# Example usage
adducts = {
    "[M+H]+": {"mass": 1.007276, "charge": 1},
    "[M+Na]+": {"mass": 22.989218, "charge": 1},
    "[M+K]+": {"mass": 38.963158, "charge": 1},
    "[M+NH4]+": {"mass": 18.033823, "charge": 1},
    "[M-H]-": {"mass": -1.007276, "charge": -1},
    "[M+Cl]-": {"mass": 34.969402, "charge": -1},
    "[M+H-H2O]+": {"mass": -18.010565 + 1.007276, "charge": 1},
    "[M+CH3OH+H]+": {"mass": 32.026215 + 1.007276, "charge": 1},
    "[M+HCOO]-": {"mass": 44.998201, "charge": -1},
    "[M+ACN+H]+": {"mass": 41.026547 + 1.007276, "charge": 1},
    "[M+2H]+": {"mass": 2.014552, "charge": 2},
    "[M-H2O+H]+": {"mass": -18.010565 + 1.007276, "charge": 1},
    "[M+Na-2H]-": {"mass": 22.989218 - 2.014552, "charge": -1},
    "[M+H+K]+": {"mass": 39.970434, "charge": 2},
    "[M+H+Na]+": {"mass": 23.996494, "charge": 2},
    "[M+2Na-H]+": {"mass": 45.978436 - 1.007276, "charge": 1},
    "[M+2Na]+": {"mass": 45.978436, "charge": 2},
    "[M-2H]-": {"mass": -1.007276, "charge": -2}
}




# WE read in the ymdb-compounds and get their hmdb_id?

ymdb_compounds = pd.read_csv('/Users/danbru/Downloads/yeast-compounds-YMDB-2025-02-02')

import requests
from bs4 import BeautifulSoup
import time

def get_hmdb_id(ymdb_id):
    # Construct the URL
    url = f"https://www.ymdb.ca/compounds/{ymdb_id}"
    try:
        # Send a GET request
        response = requests.get(url)
        response.raise_for_status()  # Raise an HTTPError for bad responses (4xx and 5xx)

        # Parse the JSON embedded in the response
        soup = response.text.strip()
        
        # Convert the soup to JSON
        data = json.loads(soup)

        # Extract the HMDB_ID
        hmdb_id = data.get("hmdb_id", None)
        return hmdb_id
    except requests.RequestException as e:
        print(f"Error fetching data for YMDB ID {ymdb_id}: {e}")
        return None
    except json.JSONDecodeError as e:
        print(f"Error decoding JSON for YMDB ID {ymdb_id}: {e}")
        return None


# Example usage with a list of YMDB IDs
#ymdb_ids = ['YMDB00001', 'YMDB00002', 'YMDB00003']  # Replace with your actual list of YMDB IDs

ymdb_ids = list(ymdb_compounds['MET_ID'])
#ymdb_ids = ['YMDB00012']

results = {}
for ymdb_id in ymdb_ids:
    hmdb_id = get_hmdb_id(ymdb_id)
    results[ymdb_id] = hmdb_id
    print(f"YMDB ID: {ymdb_id}, HMDB ID: {hmdb_id}")
    
    # Add a 1-second delay
    time.sleep(1)

#hmdb_ids = [get_hmdb_id(ymdb_id) for ymdb_id in ymdb_ids]
hmdb_ids = [value for value in results.values() if value is not None]


# Matched via pubchem and HMDB
'''
hmdb_ids = ["HMDB0000161", # Alanine
            "HMDB0000517", # Arginine
            "HMDB0000168", # Asparagine
            "HMDB0000191", # Aspartic acid
            "HMDB0000148", # Glutamic acid
            "HMDB0000641", # Glutamine
            "HMDB0000123", # Glycine
            "HMDB0000177", # Histidine
            "HMDB0000172", # Isoleucine
            "HMDB0000687", # Leucine
            "HMDB0000182", # Lysine
            "HMDB0000696", # Methionine
            "HMDB0000159", # Phenylalanine
            "HMDB0000162", # Prolinne
            "HMDB0000187", # Serine
            "HMDB0000167", # Threonine
            "HMDB0000929", # Tryptophan
            "HMDB0000158", # Tyrosine
            "HMDB0000883", # Valine
            "HMDB0006471", # Methylisocitric acid
            "HMDB0000125", # Glutathione
            "HMDB0006899", # Ala-Gly
            "HMDB0000193", # Isocitrate
            "HMDB0000227", # Mevalonic acid
            "HMDB0001375", # 3-Hydroxy-3-methylglutaryl-CoA
            ]  

'''
#hmdb_ids = results.items()


hmdb_ids = ["HMDB00159"]

df = scrape_metabolites(hmdb_ids, output_csv="/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/imCCS/ymdb_lib2.csv")
df = df.dropna()


df = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/imCCS/ymdb_lib.csv')
df.columns = ['Name','Monoisotopic Mass','Adduct Type','CCS Type','Data Source / Predictor','m/z','CCS Value (Å²)',
              'Reference','InChI Key','Chemical Formula']

df2 = pd.read_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/imCCS/ymdb_lib2.csv')
df2.columns = ['Name','Monoisotopic Mass','Adduct Type','CCS Type','Data Source / Predictor','m/z','CCS Value (Å²)',
              'Reference','InChI Key','Chemical Formula']


df = pd.concat([df, df2], axis = 0)
df = df.drop_duplicates()
df['Reference'] = df['Reference'].astype(str)
df['Adduct Type'] = df['Adduct Type'].replace('[M-2H]-', '[M-2H]2-')




# Add data from the unified CCS compendium?
ccs_compendium = pd.read_csv('/Users/danbru/Downloads/UnifiedCCSCompendium_FullDataSet_2025-02-03.csv')[['InChiKey','mz','Ion.Species','Ion.Species.Agilent','CCS','Sources']]
ccs_compendium['Ion.Species'] = ccs_compendium['Ion.Species'] + ccs_compendium['Ion.Species.Agilent'].str[-1]
ccs_compendium = ccs_compendium.drop('Ion.Species.Agilent', axis = 1)
ccs_compendium.columns = ['InChI Key','m/z','Adduct Type','CCS Value (Å²)','Reference']
ccs_compendium['Adduct Type'] = ccs_compendium['Adduct Type'].replace('[M-2H+Na]-', '[M+Na-2H]-')
ccs_compendium['Adduct Type'] = ccs_compendium['Adduct Type'].replace('[M+2H]2', '[M+2H]2+')
ccs_compendium['Adduct Type'] = ccs_compendium['Adduct Type'].replace('[M-2H]-', '[M-2H]2-')
ccs_compendium['Adduct Type'] = ccs_compendium['Adduct Type'].replace('[M+H-H2O]]', '[M-H2O+H]+')

ccs_compendium['CCS Type'] = 'Experimental'
ccs_compendium['Data Source / Predictor'] = 'CCS Compendium'
ccs_compendium = ccs_compendium[ccs_compendium['InChI Key'].isin(df['InChI Key'])]
ccs_compendium = ccs_compendium.merge(df[['Name','InChI Key','Monoisotopic Mass', 'Chemical Formula']], left_on  = 'InChI Key', right_on = 'InChI Key', how = 'left').drop_duplicates()

# Combine the predicted and experimental tables
df_total = pd.concat([df, ccs_compendium])


#df_total = df_total.sort_values(by=["Name", "Adduct Type", "CCS Type"], key=lambda col: col != "Experimental")

# Drop duplicates keeping only the first (Experimental) row within each group of Name and Adduct
#filtered_df = df_total.drop_duplicates(subset=["Name", "Adduct Type"], keep="first")


sorting_order_ccs_type = ["Experimental", "Predicted", "Theoretical"]
sorting_order_data_source = ["DarkChem", "DeepCCS", "AllCCS"]

sorted_df = df_total.sort_values(
    by=["Name", "Adduct Type", "CCS Type", "Reference", "Data Source / Predictor"],
    key=lambda col: (
        col.map(lambda x: sorting_order_ccs_type.index(x) if x in sorting_order_ccs_type else len(sorting_order_ccs_type))
        if col.name == "CCS Type" else (
            col.map(lambda x: sorting_order_data_source.index(x) if x in sorting_order_data_source else len(sorting_order_data_source))
            if col.name == "Data Source / Predictor" else (col != 6 if col.name == "Reference" else col)
        )
    )
)
filtered_df = sorted_df.drop_duplicates(subset=["Name", "Adduct Type"], keep="first")

df_total = filtered_df


# Format it so that it works with the transition list in Skyline?
transition_list = pd.DataFrame()
transition_list['Molecule List Name'] = df_total['Data Source / Predictor']+'-'+df_total['Reference'].str.replace(", ", "")
transition_list['Precursor Name'] = df_total['Name']
transition_list['Precursor Formula'] = df_total['Chemical Formula']
transition_list['Precursor Adduct'] = df_total['Adduct Type']
transition_list['Precursor Charge'] = transition_list['Precursor Adduct'].map(lambda x: adducts[x]['charge'] if x in adducts else None)
transition_list['Product m/z'] = df_total['m/z'].round(4)
transition_list['Product Charge'] = transition_list['Precursor Adduct'].map(lambda x: adducts[x]['charge'] if x in adducts else None)
transition_list['CCS'] = df_total['CCS Value (Å²)'].astype(float).round(1)
transition_list = transition_list.drop_duplicates()


transition_list.to_csv('/Users/danbru/Library/CloudStorage/OneDrive-Chalmers/Desktop/GenExp/data/imCCS/transitionList_new.csv')

