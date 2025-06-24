import sys
import pandas

args = sys.argv[1:]
print(args, file=sys.stderr)
if len(args) != 1:
    print("Usage: python growth_preprocessing.py <path_to_csv_file>")
    sys.exit(1)
csv_file_path = args[0]
print(csv_file_path, file=sys.stderr)

# Load in the growth data
gd = pandas.read_csv(csv_file_path, index_col=0)
gd.reset_index(inplace=True)
gd.index.name = "i"
gd.rename(columns={"index": "label"}, inplace=True)
gd["assay_number"] = [
    None,
    "4",
    None,
    "3",
    None,
    "6",
    "7",
    "8"
]
gd["link"] = [
    None,
    'Treatment=Yes, Supplement=None',
    None,
    'Treatment=None, Supplement=PosHigh',
    None,
    'Treatment=Yes, Supplement=PosHigh',
    'Treatment=No, Supplement=NegHigh',
    'Treatment=Yes, Supplement=NegHigh'
]
gd["Treatment"] = [
    None,
    'Yes',
    None,
    'None',
    None,
    'Yes',
    'No',
    'Yes'
]
gd["Supplement"] = [
    None,
    'None',
    None,
    'PosHigh',
    None,
    'PosHigh',
    'NegHigh',
    'NegHigh'
]

# label_dict = {
#     "Intercept": "",
#     f"C(Treatment)[T.{treatment}]": "",
#     f"Q(""Supplementation (per mM)"")": "",
#     f"Supplement at {treatment_parameters}": "",
#     f"C(Treatment)[T.{treatment}]:Q(""Supplementation (per mM)"")": "",
#     f"{treatment}×(Supplement at {treatment_parameters})": "",
#     f"Q(""Negative control"")": "",
#     f"C(Treatment)[T.{treatment}]:Q(""Negative control"")": ""
# }

gd.to_csv(sys.stdout)
