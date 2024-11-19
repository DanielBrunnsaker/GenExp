import rdflib
from rdflib.namespace import RDFS, RDF, OWL
import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

OBO = rdflib.Namespace('http://purl.obolibrary.org/obo/')
HYPO = rdflib.Namespace('http://hypo.project-genesis.io#')


term = rdflib.URIRef('http://purl.obolibrary.org/obo/')


def role_between_classes(a, b, r):
    bn = rdflib.BNode()
    trips = []
    trips.append((bn, RDF.type, OWL.Restriction))
    trips.append((bn, OWL.onProperty, r))
    trips.append((bn, OWL.someValuesFrom, b))
    trips.append((a, RDFS.subClassOf, bn))
    return trips



hypo_graph = rdflib.Graph()
hypo_graph.parse(os.path.join(BASE_DIR, 'ontology-files/hypo.ttl'))

def term_from_label(label):
    term = hypo_graph.value(predicate=RDFS.label,
                        object=rdflib.Literal(label, datatype=rdflib.URIRef('http://www.w3.org/2001/XMLSchema#string')))
    if term == None:
        term = hypo_graph.value(predicate=RDFS.label,
                        object=rdflib.Literal(label, lang='en'))
    if term == None:
        term = hypo_graph.value(predicate=RDFS.label,
                        object=rdflib.Literal(label))
    return term


"""
Add A subClassOf B

hypo_graph.add((HYPO.A, RDFS.subClassOf, HYPO.B))

Add A <exists>rel.B
for t in role_between_classes(OBO.A, HYPO.B, HYPO.rel):
    hypo_graph.add((t))

"""



"""
for sparql store:
this graph can be treated as a normal rdflib graph - (with caveat that I haven't written to graph like this, I think it should work as long as backend server allows for it)
"""

from rdflib.plugins.stores import sparqlstore

kg_endpoint='http://localhost:3030/kg'
sp_store = sparqlstore.SPARQLStore(kg_endpoint)
kg = rdflib.Graph(store=sp_store)