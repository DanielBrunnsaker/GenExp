from rdflib import Graph, Literal, URIRef, BNode
import rdflib.graph
from rdflib.plugins.stores import sparqlstore
import rdflib
from rdflib.namespace import RDFS, RDF, OWL
import os
import requests
import json


BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

OBO = rdflib.Namespace('http://purl.obolibrary.org/obo/')
HYPO = rdflib.Namespace('http://hypo.project-genesis.io#')
OBOINOWL = rdflib.Namespace('http://www.geneontology.org/formats/oboInOwl#')

term = rdflib.URIRef('http://purl.obolibrary.org/obo/')


def role_between_classes(a, b, r):
    bn = rdflib.BNode()
    trips = []
    trips.append((bn, RDF.type, OWL.Restriction))
    trips.append((bn, OWL.onProperty, r))
    trips.append((bn, OWL.someValuesFrom, b))
    trips.append((a, RDFS.subClassOf, bn))
    return trips


def term_from_label(label, g):
    term = g.value(predicate=RDFS.label,
                            object=rdflib.Literal(label, datatype=rdflib.URIRef('http://www.w3.org/2001/XMLSchema#string')))
    if term == None:
        term = g.value(predicate=RDFS.label,
                                object=rdflib.Literal(label, lang='en'))
    if term == None:
        term = g.value(predicate=RDFS.label,
                                object=rdflib.Literal(label))
    if term == None:
        term = hypo_graph.value(predicate=OBOINOWL.hasExactSynonym,
                                object=rdflib.Literal(label))
    if term == None:
        term = hypo_graph.value(predicate=OBOINOWL.hasRelatedSynonym,
                                object=rdflib.Literal(label))
    return term

hypo_graph = rdflib.Graph()
hypo_graph.parse(os.path.join(BASE_DIR, 'ontology-files/hypo.ttl'), format='turtle')

# Set up Fuseki endpoint to create new SPARQL store
host = "localhost"
chebi_fuseki_port = 3030
try:
    requests.get(f"http://{host}:{chebi_fuseki_port}/hypotheses/query")
except requests.exceptions.ConnectionError:
    host = "host.docker.internal"
    requests.get(f"http://{host}:{chebi_fuseki_port}/hypotheses/query")

chebi_query_endpoint = f'http://{host}:{chebi_fuseki_port}/hypotheses/query'
chebi_update_endpoint = f'http://{host}:{chebi_fuseki_port}/hypotheses/update'
chebi_store = sparqlchebi_store.SPARQLStore()
chebi_store.open((chebi_query_endpoint, chebi_update_endpoint))

chebi = Graph()

"""
Add A subClassOf B

hypo_graph.add((HYPO.A, RDFS.subClassOf, HYPO.B))

Add A <exists>rel.B
for t in role_between_classes(OBO.A, HYPO.B, HYPO.rel):
    hypo_graph.add((t))

"""

hypotheses_dir = os.path.join(BASE_DIR, "experiments/plans")

hypotheses = []

for folder in os.listdir(hypotheses_dir):
    with open(os.path.join(hypotheses_dir, folder, "hypothesis/hypothesis_details.json"), "r") as fi:
        hypotheses.append(json.load(fi))

CHEM_COMP_ACC = term_from_label("chemical compound accumulation", hypo_graph)

def hypothesis_to_triples(h):
    triples = []
    # Match the amino acid
    amino_acid = term_from_label()

    # Generate `observable` concepts for the reference
    triples.append((amino_acid, RDF.type, CHEM_COMP_ACC))
    triples.extend(role_between_classes())
    # Split the formula