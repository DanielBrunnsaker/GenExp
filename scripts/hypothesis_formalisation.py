import rdflib
from rdflib.namespace import RDFS, RDF, OWL

OBO = rdflib.Namespace('http://purl.obolibrary.org/obo/')


def role_between_classes(a, b, r):
    bn = rdflib.BNode()
    trips = []
    trips.append((bn, RDF.type, OWL.Restriction))
    trips.append((bn, OWL.onProperty, r))
    trips.append((bn, OWL.someValuesFrom, b))
    trips.append((a, RDFS.subClassOf, bn))
    return trips
