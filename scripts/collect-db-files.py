#!/usr/bin/env python3
import os
import sys
from rdflib import Dataset

def collect_trig_files(root_dir):
    trig_files = []
    for dirpath, _, filenames in os.walk(root_dir):
        for fname in filenames:
            if fname.endswith(".trig"):
                trig_files.append(os.path.join(dirpath, fname))
    return trig_files

def merge_trig_files(root_dir, output_file):
    ds = Dataset()
    trig_files = collect_trig_files(root_dir)
    for trig_file in trig_files:
        print(f"Parsing {trig_file} ...", file=sys.stderr)
        ds.parse(trig_file, format="trig")
    ds.serialize(destination=output_file, format="trig")
    print(f"Merged dataset written to {output_file} ({len(ds)} quads)", file=sys.stderr)

if __name__ == "__main__":
    top_dir = os.getcwd()
    output_trig = os.path.join(top_dir, "experiments" ,"merged_dataset.trig")
    merge_trig_files(os.path.join(top_dir,"experiments"), output_trig)