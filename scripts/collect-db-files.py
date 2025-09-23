#!/usr/bin/env python3
import os
import sys
from rdflib import Dataset

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

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

    ontology_file = os.path.join(BASE_DIR, "ontology-files", "hypo.ttl")
    print(f"Parsing ontology {ontology_file} ...", file=sys.stderr)
    ontology_graph = ds.graph(identifier="http://hypo.project-genesis.io/ontology")
    ontology_graph.parse(ontology_file, format="turtle")
    
    ds.serialize(destination=output_file, format="trig")
    print(f"Merged dataset written to {output_file} ({len(ds)} quads)", file=sys.stderr)

if __name__ == "__main__":
    top_dir = os.getcwd()
    output_trig = os.path.join(top_dir, "experiments" ,"merged_dataset.trig")
    merge_trig_files(os.path.join(top_dir,"experiments"), output_trig)