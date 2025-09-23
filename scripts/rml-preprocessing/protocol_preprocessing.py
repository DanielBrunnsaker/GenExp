#!/usr/bin/env python3
import json
import sys
import os
from tqdm.auto import tqdm
sys.path.insert(0, os.path.abspath(
    os.path.join(os.path.dirname(__file__), '../..')))

from scripts.hypothesis_formalisation import (
    CHEBI,
    term_from_label
)

# Load the JSON file given as a command line argument
args = sys.argv[1:]
# print(args, file=sys.stderr)
if len(args) != 1:
    print("Usage: python protocol_preprocessing.py <path_to_json_file>", file=sys.stderr)
    sys.exit(1)
json_file_path = args[0]
print(json_file_path, file=sys.stderr)
try:
    with open(json_file_path, 'r') as file:
        data = json.load(file)
        print(f"Loaded JSON data from {json_file_path}", file=sys.stderr)
except FileNotFoundError:
    print(f"Error: The file {json_file_path} does not exist.", file=sys.stderr)
    sys.exit(1)
except json.JSONDecodeError:
    print(f"Error: The file {json_file_path} is not a valid JSON file.", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"An unexpected error occurred: {e}", file=sys.stderr)
    sys.exit(1)

# Create a new data structure to hold the modified data

supp_stats_keys = {
    1: {"treatment": None, "supplement": None, "text": "Treatment=None, Supplement=None"},
    2: {"treatment": None, "supplement": "PosLow", "text": "Treatment=None, Supplement=PosLow"},
    3: {"treatment": None, "supplement": "PosHigh", "text": "Treatment=None, Supplement=PosHigh"},
    4: {"treatment": "Yes", "supplement": None, "text": "Treatment=Yes, Supplement=None"},
    5: {"treatment": "Yes", "supplement": "PosLow", "text": "Treatment=Yes, Supplement=PosLow"},
    6: {"treatment": "Yes", "supplement": "PosHigh", "text": "Treatment=Yes, Supplement=PosHigh"},
    7: {"treatment": "No", "supplement": "NegHigh", "text": "Treatment=No, Supplement=NegHigh"},
    8: {"treatment": "Yes", "supplement": "NegHigh", "text": "Treatment=Yes, Supplement=NegHigh"},
}


for (i, exp) in enumerate(data['experiments']):
    # Add identifier that can link to the growth data tests

    # Add identifier that can link to the metabolomics comparisons
    exp["metabolomics_comparison_id"] = supp_stats_keys.get(exp["experiment"])

    # If the `media_supplementation`field is not None, convert it to a Chebi term
    if exp["media_supplementation"] is not None:
        # Convert the media_supplementation to a Chebi term
        chebi_term = term_from_label(
            f"L-{exp['media_supplementation']}", CHEBI)
        if chebi_term is None:  # Could be a case issue
            print(
                f"Chebi term 'L-{exp['media_supplementation']} not found, trying 'L-{exp['media_supplementation'].lower()}'", file=sys.stderr)
            chebi_term = term_from_label(
            f"L-{exp['media_supplementation'].lower()}", CHEBI)
        if chebi_term is None:  # Could be not an amino acid
            print(
                f"Chebi term 'L-{exp['media_supplementation']} not found, trying '{exp['media_supplementation']}'", file=sys.stderr)
            chebi_term = term_from_label(exp['media_supplementation'], CHEBI)
        if chebi_term is None:  # Could be not an amino acid
            print(
                f"Chebi term '{exp['media_supplementation']} not found, trying '{exp['media_supplementation'].lower()}'", file=sys.stderr)
            chebi_term = term_from_label(exp['media_supplementation'].lower(), CHEBI)
        if chebi_term is None and exp['media_supplementation'].lower() == "aminoadipate":  # Aminoadipate
            chebi_term = term_from_label("L-2-aminoadipate(2-)", CHEBI)
        exp["media_supplementation_chebi"] = chebi_term

    # Separate the `media_supplementation_doses` field into the value and the unit
    if exp["media_supplementation_doses"] is not None:
        if '%' in exp["media_supplementation_doses"]:
            i = exp["media_supplementation_doses"].find("%")
            exp["media_supplementation_doses_value"] = exp["media_supplementation_doses"][:i]
            exp["media_supplementation_doses_unit"] = exp["media_supplementation_doses"][i:]
        else:
            exp["media_supplementation_doses_value"], exp["media_supplementation_doses_unit"] = (
                exp["media_supplementation_doses"].split(" ")
            )

    # Same for `treatment_parameters`
    if exp["treatment_parameters"] is not None:
        if '%' in exp["treatment_parameters"]:
            i = exp["treatment_parameters"].find("%")
            exp["treatment_parameters_value"] = exp["treatment_parameters"][:i]
            exp["treatment_parameters_unit"] = exp["treatment_parameters"][i:]
        else:
            exp["treatment_parameters_value"], exp["treatment_parameters_unit"] = (
                exp["treatment_parameters"].split(" ")
            )

    # If the `treatment` field is not None, convert it to a Chebi term
    if exp["treatment"] is not None:
        # Convert the treatment to a Chebi term
        chebi_term = term_from_label(exp["treatment"], CHEBI)
        exp["treatment_chebi"] = chebi_term

    exp["media_description"] = exp["baseline_media"] + (
        "" if exp["media_supplementation"] is None else f"; supplementation: {exp['media_supplementation']} {exp['media_supplementation_doses']}"
    )

# Write out the modified JSON to stdout
print(json.dumps(data, indent=4))
