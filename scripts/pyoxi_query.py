import pyoxigraph as pox
import sys
import json
import pandas as pd
from pprint import pprint

import builtins
import traceback

def rprint(*objs, **kwargs):
    my_prefix = len(traceback.format_stack())*" "
    builtins.print(my_prefix, *objs, **kwargs)

store = pox.Store()
store.bulk_load(path="experiments/merged_dataset.trig", format=pox.RdfFormat.TRIG)

def transpose_tsv(tsv_string):
    return '\n'.join(map(lambda t : '\t'.join(t), zip(*map(lambda s : s.split('\t'),tsv_string.split('\n')))))


query = sys.stdin.read()
results = store.query(query)
if isinstance(results, pox.QueryBoolean):
    rprint(bool(results))
elif isinstance(results, pox.QuerySolutions):
    try:
        # pprint(json.loads(results.serialize(format=pox.QueryResultsFormat.JSON).decode()))
        tsv_string = results.serialize(format=pox.QueryResultsFormat.TSV).decode().strip()
        for ln in transpose_tsv(tsv_string).split('\n'):
            rprint(ln)

    except AttributeError:
        rprint("No results")