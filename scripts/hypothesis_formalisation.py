from rdflib import Graph, Literal, URIRef, BNode
import rdflib.graph
from rdflib.plugins.stores import sparqlstore
import rdflib
from rdflib.namespace import RDFS, RDF, OWL
import os
import requests
import json
import re
import uuid6
from typing import Optional
from pprint import pprint

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

OBO = rdflib.Namespace('http://purl.obolibrary.org/obo/')
HYPO = rdflib.Namespace('http://hypo.project-genesis.io#')
OBOINOWL = rdflib.Namespace('http://www.geneontology.org/formats/oboInOwl#')

term = rdflib.URIRef('http://purl.obolibrary.org/obo/')


def create_id(prefix: Optional[str] = None):
    if prefix is not None:
        return prefix + "-{}".format(uuid6.uuid7())
    else:
        return "{}".format(uuid6.uuid7())


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
        term = g.value(predicate=OBOINOWL.hasExactSynonym,
                       object=rdflib.Literal(label))
    if term == None:
        term = g.value(
            predicate=OBOINOWL.hasRelatedSynonym,
            object=rdflib.Literal(label)
        )
    return term


hypo_graph = rdflib.Graph()
hypo_graph.parse(os.path.join(
    BASE_DIR, 'ontology-files/hypo.ttl'), format='turtle')

# Set up Fuseki endpoint to create new SPARQL store
host = "localhost"
chebi_fuseki_port = 3030
try:
    requests.get(f"http://{host}:{chebi_fuseki_port}/genesis/query")
except requests.exceptions.ConnectionError:
    # Perhaps we are in a docker container (ensure that the correct port is forwarded)
    host = "host.docker.internal"
    requests.get(f"http://{host}:{chebi_fuseki_port}/genesis/query")

chebi_query_endpoint = f'http://{host}:{chebi_fuseki_port}/genesis/query'
chebi_store = sparqlstore.SPARQLStore(chebi_query_endpoint)

chebi = Graph(chebi_store, rdflib.graph.DATASET_DEFAULT_GRAPH_ID)

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
CHEM_ACC_OF = term_from_label("accumulationOfChemical", hypo_graph)
STATE_HAS_OBSERVABLE = term_from_label("stateHasObservable", hypo_graph)
REFERENCE_STATE = HYPO[create_id(prefix="S-REF")]
ORGANISM_STATE = term_from_label("organismState", hypo_graph)
hypo_graph.add((REFERENCE_STATE, RDFS.subClassOf, ORGANISM_STATE))

test_lp = "Cell(A):-exhibits_phenotype(A,'decreased metal resistance',B,C),compound_name(B,'zinc dichloride')."
ATOM_PARSER = re.compile(r"^(\w+)\((.+)\)$")


def find_or_create_node(
    graph,
    triples,
    id_prefix: Optional[str] = None,
    namespace: Optional[rdflib.Namespace] = None
):
    '''Given a graph, and a set of triples to uniquely define a node in
    that graph, either return the node if it exists, or create the node.

    Specify
        - `id_prefix` (a string) and 
        - `namespace` (an `rdflib.Namespace`)

    The `triples` should be a set of triples with a NoneType value in 
    the position to be replaced.
    '''

    QUERY_TEMPLATE = """SELECT {} WHERE {{
        {}
    }}"""

    query_triples = [tuple(map(lambda v: "?queried_node" if v is None else v, t))
                     for t in triples]
    query = QUERY_TEMPLATE.format(
        "?queried_node",
        "\n\t".join(map(lambda t: " ".join(map(lambda a: "?{}".format(a) if isinstance(
            a, rdflib.BNode) else a if a.startswith("?") else "<{}>".format(a), t)) + " .", query_triples))
    )
    qres = graph.query(query)

    if len(qres) == 0:
        # print("Need to add this node.")
        node_uri = namespace[create_id(prefix=id_prefix)]
        node_triples = [
            tuple(map(lambda v: node_uri if v is None else v, t)) for t in triples]
        for t in node_triples:
            graph.add(t)
        return node_uri
    else:
        assert len(qres) == 1
        return qres.__iter__().__next__().get("queried_node")


def logic_program_to_state(logic_program):
    triples = []
    atoms = []
    head, body = test_lp.split(":-")
    while len(body) > 0:
        end = body.find(")")
        predicate, arguments = ATOM_PARSER.match(body[:end+1]).groups()
        atoms.append((predicate, [a.replace("'", "")
                     for a in arguments.split(",")]))
        if end < len(body):
            body = body[end+2:]

    stateid = create_id(prefix="S")

    exhibits_phenotype = filter(lambda t: t[0] == "exhibits_phenotype", atoms)
    phenotype_uri_list = []
    for _, args in exhibits_phenotype:
        phid = create_id(prefix="P")
        qualifier, phtype = args[1].split(" ", maxsplit=1)

        phtype = term_from_label(phtype, hypo_graph)

        # Fetch the reference phenotype from the graph if it exists,
        # create it if it doesn't exist.
        ref_query_trips = role_between_classes(
            REFERENCE_STATE, None, STATE_HAS_OBSERVABLE) + [(None, RDFS.subClassOf, phtype)]
        ref_phtype = find_or_create_node(
            hypo_graph, ref_query_trips, id_prefix="P-REF", namespace=HYPO)
        # print(ref_phtype)

        # Fetch the compared phenotype from the graph if it exists,
        # create it if it doesn't exist.
        qualifying_relation = term_from_label(
            qualifier + "ComparedTo", hypo_graph)
        comp_query_trips = role_between_classes(
            None, ref_phtype, qualifying_relation) + [(None, RDFS.subClassOf, phtype)]
        comp_phtype = find_or_create_node(
            hypo_graph, comp_query_trips, id_prefix="P", namespace=HYPO)
        # print(comp_phtype)
        phenotype_uri_list.append(comp_phtype)

        # triples.extend(role_between_classes(
        #     HYPO[stateid], HYPO[phid], STATE_HAS_OBSERVABLE))
        # triples.append((HYPO[phid], RDFS.subClassOf, phtype))
        # triples.append((HYPO[phid], qualifying_relation,
        #                ref_phtype))

    # Fetch the complete state from the graph if it exists,
    # create it if it doesn't exist.
    state_triples = [t for uri in phenotype_uri_list for t in role_between_classes(
        None, uri, STATE_HAS_OBSERVABLE)] + [(None, RDFS.subClassOf, ORGANISM_STATE)]
    state = find_or_create_node(
        hypo_graph, state_triples, id_prefix="S", namespace=HYPO
    )
    return state


def hypothesis_to_triples(h):
    triples = []

    # Match the amino acid
    amino_acid = term_from_label(h["observable"], chebi)  # Case-sensitive?
    print(amino_acid)

    # Generate `observable` concepts for the reference
    triples.append((HYPO.P1, RDF.type, CHEM_COMP_ACC))
    triples.append((HYPO.P1, CHEM_ACC_OF, amino_acid))
    # triples.extend(role_between_classes())
    # Split the formula

    return triples


