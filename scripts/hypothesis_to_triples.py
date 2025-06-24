# imports
import pickle
from pprint import pprint

# tests
with open("results/patterns/feature_dict.pickle", "rb") as fi:
    fd = pickle.load(fi)
pprint(fd)

# input : (directionality,target,clause)


# output : triples in .ttl format