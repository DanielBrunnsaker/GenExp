from rdflib import Graph, Literal, URIRef, BNode
import rdflib.graph
from rdflib.plugins.stores import sparqlstore
import rdflib
from rdflib.namespace import RDFS, RDF, OWL
import os
import requests

"""
for sparql store:
this graph can be treated as a normal rdflib graph - (with caveat that I haven't written to graph like this, I think it should work as long as backend server allows for it)
"""


def my_bnode_ext(node):
    if isinstance(node, BNode):
        return f"<bnode:b{node}>"
    return sparqlstore._node_to_sparql(node)


# Set up Fuseki endpoint to create new SPARQL store
host = "localhost"
fuseki_port = 3030
try:
    requests.get(f"http://{host}:{fuseki_port}/hypotheses/query")
except requests.exceptions.ConnectionError:
    host = "host.docker.internal"
    requests.get(f"http://{host}:{fuseki_port}/hypotheses/query")

query_endpoint = f'http://{host}:{fuseki_port}/hypotheses/query'
update_endpoint = f'http://{host}:{fuseki_port}/hypotheses/update'
store = sparqlstore.SPARQLUpdateStore(node_to_sparql=my_bnode_ext)
store.open((query_endpoint, update_endpoint))