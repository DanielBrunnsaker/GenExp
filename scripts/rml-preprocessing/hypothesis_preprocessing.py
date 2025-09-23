#!/usr/bin/env python3
import json
import sys
import os


# Load the JSON file given as a command line argument
args = sys.argv[1:]
# print(args, file=sys.stderr)
if len(args) != 1:
    print("[!ERR] Usage: python hypothesis_preprocessing.py <path_to_selected_hypothesis_dir>", file=sys.stderr)
    sys.exit(1)
selected_hypothesis_directory = args[0]
json_file_path = os.path.join(
    selected_hypothesis_directory, "hypothesis_details.json")
hypothesis_text_path = os.path.join(selected_hypothesis_directory, "hypothesis.txt")
# print(json_file_path, file=sys.stderr)
try:
    with open(json_file_path, 'r') as file:
        data = json.load(file)
        print(f"[INFO] Loaded JSON data from {json_file_path}", file=sys.stderr)
except FileNotFoundError:
    print(f"[WARN] The file {json_file_path} does not exist.", file=sys.stderr)
    sys.exit()
except json.JSONDecodeError:
    print(f"[!ERR] The file {json_file_path} is not a valid JSON file.", file=sys.stderr)
    sys.exit(1)
except Exception as e:
    print(f"[!ERR] An unexpected error occurred: {e}", file=sys.stderr)
    sys.exit(1)

try:
    with open(hypothesis_text_path, "r") as fi:
        for i, ln in enumerate(fi):
            if i == 4:
                hypothesis_text = ln.rstrip()
                data["hypothesis_text"] = hypothesis_text
                break
except FileNotFoundError:
    data["hypothesis_text"] = ""

# Write out the modified JSON to stdout
print(json.dumps(data, indent=4))
