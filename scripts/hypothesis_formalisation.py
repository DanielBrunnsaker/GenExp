from copy import copy, deepcopy
from itertools import groupby
import os
import json
import re
from typing import Optional
from pprint import pprint

from rdflib import Graph, Literal, URIRef, BNode
import rdflib.graph
from rdflib.plugins.stores import sparqlstore
import rdflib
from rdflib.namespace import RDFS, RDF, OWL
import requests
import uuid6

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

OBO = rdflib.Namespace("http://purl.obolibrary.org/obo/")
HYPO = rdflib.Namespace("http://hypo.project-genesis.io#")
OBOINOWL = rdflib.Namespace("http://www.geneontology.org/formats/oboInOwl#")

# term = rdflib.URIRef("http://purl.obolibrary.org/obo/")


def create_id(prefix: Optional[str] = None):
    if prefix is not None:
        return prefix + f"-{uuid6.uuid7()}"
    else:
        return f"{uuid6.uuid7()}"


def role_between_classes(a, b, r):
    bn = rdflib.BNode()
    trips = []
    trips.append((bn, RDF.type, OWL.Restriction))
    trips.append((bn, OWL.onProperty, r))
    trips.append((bn, OWL.someValuesFrom, b))
    trips.append((a, RDFS.subClassOf, bn))
    return trips


def term_from_label(label, g):
    term = g.value(
        predicate=RDFS.label,
        object=rdflib.Literal(
            label, datatype=rdflib.URIRef("http://www.w3.org/2001/XMLSchema#string")
        ),
    )
    if term is None:
        term = g.value(predicate=RDFS.label, object=rdflib.Literal(label, lang="en"))
    if term is None:
        term = g.value(predicate=RDFS.label, object=rdflib.Literal(label))
    if term is None:
        term = g.value(predicate=OBOINOWL.hasExactSynonym, object=rdflib.Literal(label))
    if term is None:
        term = g.value(
            predicate=OBOINOWL.hasRelatedSynonym, object=rdflib.Literal(label)
        )
    return term


# Set up Fuseki endpoint to create new SPARQL store
HOST = "localhost"
CHEBI_FUSEKI_PORT = 3030
try:
    requests.get(f"http://{HOST}:{CHEBI_FUSEKI_PORT}/genesis/query", timeout=10)
except requests.exceptions.ConnectionError:
    # Perhaps we are in a docker container (ensure that the correct port is forwarded)
    HOST = "host.docker.internal"
    requests.get(f"http://{HOST}:{CHEBI_FUSEKI_PORT}/genesis/query", timeout=10)

chebi_query_endpoint = f"http://{HOST}:{CHEBI_FUSEKI_PORT}/genesis/query"
chebi_store = sparqlstore.SPARQLStore(chebi_query_endpoint)

chebi = Graph(chebi_store, rdflib.graph.DATASET_DEFAULT_GRAPH_ID)

# Create a local dataset for hypotheses
hypo_ds = rdflib.Dataset()
hypo_ds.bind("hypo", HYPO)

# Create local graph for hypothesis ontology
hypo_graph = hypo_ds.graph(rdflib.URIRef("http://hypo.project-genesis.io"))
hypo_graph.parse(os.path.join(BASE_DIR, "ontology-files/hypo.ttl"), format="turtle")

"""
Add A subClassOf B

hypo_graph.add((HYPO.A, RDFS.subClassOf, HYPO.B))

Add A <exists>rel.B
for t in role_between_classes(OBO.A, HYPO.B, HYPO.rel):
    hypo_graph.add((t))

"""

hypotheses_dir = os.path.join(BASE_DIR, "experiments")

hypotheses = []

for folder in os.listdir(hypotheses_dir):
    try:
        with open(
            os.path.join(hypotheses_dir, folder, "hypothesis/hypothesis_details.json"),
            "r",
            encoding="utf-8",
        ) as fi:
            hypotheses.append(json.load(fi))
    except (FileNotFoundError, NotADirectoryError):
        continue

CHEM_COMP_ACC = term_from_label("chemical compound accumulation", hypo_graph)
CHEM_ACC_OF = term_from_label("accumulationOfChemical", hypo_graph)
RESISTANCE_TO_CHEM = term_from_label("resistanceToChemical", hypo_graph)
METAL_RESISTANCE = term_from_label("metalResistance", hypo_graph)
STATE_HAS_OBSERVABLE = term_from_label("stateHasObservable", hypo_graph)
DEFAULT_REFERENCE_STATE = HYPO[create_id(prefix="S-REF")]
ORGANISM_STATE = term_from_label("organismState", hypo_graph)
hypo_graph.add((DEFAULT_REFERENCE_STATE, RDFS.subClassOf, ORGANISM_STATE))

TEST_LOGIC_PROGRAMS = [
    "Cell(A):-exhibits_phenotype(A,'decreased metal resistance',B,C),compound_name(B,'zinc dichloride').",
    "Cell(A):-exhibits_phenotype(A,'increased resistance to chemicals',B,C),compound_name(B,fenpropimorph).",
]
ATOM_PARSER = re.compile(r"^(\w+)\((.+)\)$")


def find_or_create_node(
    graph,
    triples,
    id_prefix: Optional[str] = None,
    namespace: Optional[rdflib.Namespace] = None,
):
    """Given a graph, and a set of triples to uniquely define a node in
    that graph, either return the node if it exists, or create the node.

    Specify
        - `id_prefix` (a string) and
        - `namespace` (an `rdflib.Namespace`)

    The `triples` should be a set of triples with a NoneType value in
    the position to be replaced.
    """

    query_template = """SELECT {} WHERE {{
        {}
    }}"""

    query_triples = [
        tuple(map(lambda v: "?queried_node" if v is None else v, t)) for t in triples
    ]
    query = query_template.format(
        "?queried_node",
        "\n\t".join(
            map(
                lambda t: " ".join(
                    map(
                        lambda a: (
                            f"?{a}"
                            if isinstance(a, rdflib.BNode)
                            else a if a.startswith("?") else f"<{a}>"
                        ),
                        t,
                    )
                )
                + " .",
                query_triples,
            )
        ),
    )
    qres = graph.query(query)

    if len(qres) == 0:
        # print("Need to add this node.")
        node_uri = namespace[create_id(prefix=id_prefix)]
        node_triples = [
            tuple(map(lambda v: node_uri if v is None else v, t)) for t in triples
        ]
        for t in node_triples:
            graph.add(t)
        return node_uri
    else:
        assert len(qres) == 1
        return next(iter(qres)).get("queried_node")


def logic_program_to_state(logic_program, reference_state=DEFAULT_REFERENCE_STATE):
    triples = []
    atoms = []
    head, body = logic_program.split(":-")
    while len(body) > 0:
        end = body.find(")")
        predicate, arguments = ATOM_PARSER.match(body[: end + 1]).groups()
        atoms.append((predicate, [a.replace("'", "") for a in arguments.split(",")]))
        if end < len(body):
            body = body[end + 2 :]

    def keyfunc(t):
        return t[0]

    atoms = {k: list(g) for k, g in groupby(sorted(atoms, key=keyfunc), keyfunc)}
    phenotype_uri_list = [
        write_exhibits_phenotype_to_graph(args, atoms, hypo_graph)
        for _, args in atoms.get("exhibits_phenotype", [])
    ]

    # Fetch the complete state from the graph if it exists,
    # create it if it doesn't exist.
    state_triples = [
        t
        for uri_pair in phenotype_uri_list
        for t in role_between_classes(None, uri_pair[1], STATE_HAS_OBSERVABLE)
    ] + [(None, RDFS.subClassOf, ORGANISM_STATE)]
    state = find_or_create_node(
        hypo_graph, state_triples, id_prefix="S", namespace=HYPO
    )
    # TODO: Add Label or other property to make querying easier
    return state


def write_exhibits_phenotype_to_graph(
    args,
    atoms,
    graph,
    reference_state=DEFAULT_REFERENCE_STATE,
    phenotype_namespace=HYPO,
):
    # Get qualifier and phenotype from arguments
    qualifier, phtype_label = args[1].split(" ", maxsplit=1)
    phtype = term_from_label(phtype_label, graph)
    qualifying_relation = term_from_label(qualifier + "ComparedTo", graph)

    match phtype_label:
        case "resistance to chemicals":
            compound_name = term_from_label(
                next(filter(lambda t: t[1][0] == args[2], atoms.get("compound_name")))[
                    1
                ][1],
                chebi,
            )
            additional_ref_triples = role_between_classes(
                None, compound_name, RESISTANCE_TO_CHEM
            )
            additional_comp_triples = role_between_classes(
                None, compound_name, RESISTANCE_TO_CHEM
            )
        case "metal resistance":
            compound_name = term_from_label(
                next(filter(lambda t: t[1][0] == args[2], atoms.get("compound_name")))[
                    1
                ][1],
                chebi,
            )
            additional_ref_triples = role_between_classes(
                None, compound_name, METAL_RESISTANCE
            )
            additional_comp_triples = role_between_classes(
                None, compound_name, METAL_RESISTANCE
            )
        case _:
            additional_ref_triples = []
            additional_comp_triples = []

    # Fetch the reference phenotype from the graph if it exists,
    # create it if it doesn't exist. (Link to reference state defined
    # above in variable `REFERENCE_STATE`)
    # TODO: Add Label or other property to make querying easier
    ref_query_trips = (
        role_between_classes(reference_state, None, STATE_HAS_OBSERVABLE)
        + [(None, RDFS.subClassOf, phtype)]
        + additional_ref_triples
    )
    ref_phtype = find_or_create_node(
        hypo_graph, ref_query_trips, id_prefix="P-REF", namespace=phenotype_namespace
    )

    # Fetch the compared phenotype from the graph if it exists,
    # create it if it doesn't exist.
    # TODO: Add Label or other property to make querying easier
    comp_query_trips = (
        role_between_classes(None, ref_phtype, qualifying_relation)
        + [(None, RDFS.subClassOf, phtype)]
        + additional_comp_triples
    )
    comp_phtype = find_or_create_node(
        hypo_graph, comp_query_trips, id_prefix="P", namespace=phenotype_namespace
    )
    return ref_phtype, comp_phtype


def hypothesis_to_triples(h):
    triples = []

    # Represent the logic program
    state_1 = logic_program_to_state(h["logic_program"])

    # Match the amino acid
    amino_acid = term_from_label(h["observable"], chebi)  # Case-sensitive?
    print(amino_acid)
    # state_2 =

    # Generate `observable` concepts for the reference
    # triples.append((HYPO.P1, RDF.type, CHEM_COMP_ACC))
    # triples.append((HYPO.P1, CHEM_ACC_OF, amino_acid))
    # triples.extend(role_between_classes())
    # Split the formula

    # State implies State

    return triples
