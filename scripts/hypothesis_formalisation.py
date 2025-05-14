from copy import copy, deepcopy
import datetime
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
DCT = rdflib.Namespace("http://purl.org/dc/terms/")

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
            label, datatype=rdflib.URIRef(
                "http://www.w3.org/2001/XMLSchema#string")
        ),
    )
    if term is None:
        term = g.value(predicate=RDFS.label,
                       object=rdflib.Literal(label, lang="en"), any=False)
    if term is None:
        term = g.value(predicate=RDFS.label,
                       object=rdflib.Literal(label), any=False)
    if term is None:
        term = g.value(predicate=OBOINOWL.hasExactSynonym,
                       object=rdflib.Literal(label), any=False)
    if term is None:
        term = g.value(
            predicate=OBOINOWL.hasRelatedSynonym, object=rdflib.Literal(
                label), any=False
        )
    return term


# Set up Fuseki endpoint to create new SPARQL store
HOST = "localhost"
CHEBI_FUSEKI_PORT = 3037
try:
    requests.get(
        f"http://{HOST}:{CHEBI_FUSEKI_PORT}/genesis/query", timeout=10)
except requests.exceptions.ConnectionError:
    # Perhaps we are in a docker container (ensure that the correct port is forwarded)
    HOST = "host.docker.internal"
    requests.get(
        f"http://{HOST}:{CHEBI_FUSEKI_PORT}/genesis/query", timeout=10)

chebi_query_endpoint = f"http://{HOST}:{CHEBI_FUSEKI_PORT}/genesis/query"
chebi_store = sparqlstore.SPARQLStore(chebi_query_endpoint)

chebi = Graph(chebi_store, rdflib.graph.DATASET_DEFAULT_GRAPH_ID)

# Create a local dataset for hypotheses
hypo_ds = rdflib.Dataset()
hypo_ds.bind("hypo", HYPO)

# Create local graph for hypothesis ontology
ontology = hypo_ds.graph(
    rdflib.URIRef("http://hypo.project-genesis.io"))
ontology.parse(os.path.join(
    BASE_DIR, "ontology-files/hypo.ttl"), format="turtle")

phenotypes = hypo_ds.graph(
    rdflib.URIRef("http://hypo.project-genesis.io/phenotypes"))
states = hypo_ds.graph(
    rdflib.URIRef("http://hypo.project-genesis.io/states"))
hmeta = hypo_ds.graph(
    rdflib.URIRef("http://hypo.project-genesis.io/hypothesis-metadata"))

# Bind namespaces
for g in hypo_ds.graphs():
    g.bind("obo", OBO)
    g.bind("owl", OWL)
    g.bind("rdf", RDF)
    g.bind("rdfs", RDFS)
    g.bind("hypo", HYPO)
    g.bind("oboInOwl", OBOINOWL)

"""
Add A subClassOf B

ontology.add((HYPO.A, RDFS.subClassOf, HYPO.B))

Add A <exists>rel.B
for t in role_between_classes(OBO.A, HYPO.B, HYPO.rel):
    ontology.add((t))

"""

hypotheses_dir = os.path.join(BASE_DIR, "experiments")

hypotheses = []

for folder in os.listdir(hypotheses_dir):
    try:
        with open(
            os.path.join(hypotheses_dir, folder,
                         "hypothesis/hypothesis_details.json"),
            "r",
            encoding="utf-8",
        ) as fi:
            hypotheses.append(json.load(fi))
    except (FileNotFoundError, NotADirectoryError):
        continue

CHEM_COMP_ACC = term_from_label(
    "chemical compound accumulation", ontology)
CHEM_ACC_OF = term_from_label("accumulationOfChemical", ontology)
RESISTANCE_TO_CHEM = term_from_label("resistanceToChemical", ontology)
# For now, we are using the relation that is defined from the
# `resistance to chemicals` phenotype in APO, of which metal
# resistance is a subclass.
METAL_RESISTANCE = term_from_label("resistanceToChemical", ontology)
STATE_HAS_OBSERVABLE = term_from_label("stateHasObservable", ontology)
DEFAULT_REFERENCE_STATE = HYPO[create_id(prefix="S-REF")]
ORGANISM_STATE = term_from_label("organismState", ontology)
states.add((DEFAULT_REFERENCE_STATE, RDFS.subClassOf, ORGANISM_STATE))

TEST_LOGIC_PROGRAMS = [
    "Cell(A):-exhibits_phenotype(A,'decreased metal resistance',B,C),compound_name(B,'zinc dichloride').",
    "Cell(A):-exhibits_phenotype(A,'increased resistance to chemicals',B,C),compound_name(B,fenpropimorph).",
]
ATOM_FINDER = re.compile(r"([A-z0-9]+\(.+?\))(?:,|.|$)")
ATOM_PARSER = re.compile(r"^(\w+)\((.+)\)$")


def replace_blank_nodes_with_queried_ones(triples, query_result):
    """Replace blank nodes in the triples with the queried ones.
    This function is necessary to avoid duplication of blank nodes
    when using the `find_or_create_node` function."""
    new_trips = []
    for t in triples:
        new_t = []
        for v in t:
            if isinstance(v, rdflib.BNode):
                vq = rdflib.term.Variable(v)
                # Check if the blank node is in the query result
                if vq not in query_result.bindings[0]:
                    raise ValueError(
                        f"Blank node {v} not found in query result.")
                # Replace the blank node with the queried one
                new_t.append(query_result.bindings[0].get(vq))
            else:
                new_t.append(v)
        new_trips.append(tuple(new_t))
    return new_trips


def find_or_create_node(
    dataset: rdflib.Dataset,
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

    Returns the node URI and the triples with the node URI
    replacing the NoneType value.
    """

    query_template = """SELECT {} WHERE {{
        GRAPH ?g {{
            {}
        }}
    }}"""

    query_triples = [
        tuple(map(lambda v: "?queried_node" if v is None else v, t)) for t in triples
    ]
    query = query_template.format(
        "?queried_node ?g {}".format(" ".join(
            set([f"?{v}" for t in query_triples for v in t if isinstance(v, rdflib.BNode)]))),
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
    qres = dataset.query(query)

    if len(qres) == 0:
        # print("Need to add this node.")
        node_uri = namespace[create_id(prefix=id_prefix)]
        node_triples = [
            tuple(map(lambda v: node_uri if v is None else v, t)) for t in triples
        ]
    else:
        assert len(qres) == 1
        node_uri = next(iter(qres)).get("queried_node")
        node_triples = replace_blank_nodes_with_queried_ones([
            tuple(map(lambda v: node_uri if v is None else v, t)) for t in triples
        ], qres)

    return node_uri, node_triples


def logic_program_to_state(logic_program, h_graph, reference_state=DEFAULT_REFERENCE_STATE):
    triples = []
    atoms = []
    head, body = logic_program.split(":-", 1)
    for atom in ATOM_FINDER.findall(body):
        predicate, arguments = ATOM_PARSER.match(atom).groups()
        atoms.append((predicate, [a.replace("'", "")
                     for a in arguments.split(",")]))

    def keyfunc(t):
        return t[0]

    atoms = {k: list(g) for k, g in groupby(
        sorted(atoms, key=keyfunc), keyfunc)}
    # print(atoms)
    phenotype_uri_tuple, phen_quads_list = zip(*[
        write_exhibits_phenotype_to_graph(
            args, atoms, h_graph, reference_state=reference_state)
        for _, args in atoms.get("exhibits_phenotype", [])
    ])

    # Fetch the complete state from the graph if it exists,
    # create it if it doesn't exist.
    state_query_triples = [
        t
        for uri_pair in phenotype_uri_tuple
        for t in role_between_classes(None, uri_pair[1], STATE_HAS_OBSERVABLE)
    ] + [(None, RDFS.subClassOf, ORGANISM_STATE)]
    state, state_triples = find_or_create_node(
        hypo_ds, state_query_triples, id_prefix="S", namespace=HYPO
    )
    for t in state_triples:
        states.add(t)

    quads = [q for l in phen_quads_list for q in l] + \
        [t + (states.identifier,) for t in state_triples]

    # TODO: Add Label or other property to make querying easier
    return state, quads


def write_exhibits_phenotype_to_graph(
    args,
    atoms,
    h_graph,
    ontology=ontology,
    reference_state=DEFAULT_REFERENCE_STATE,
    phenotype_namespace=HYPO,
):
    # Get qualifier and phenotype from arguments
    qualifier, phtype_label = args[1].split(" ", maxsplit=1)
    phtype = term_from_label(phtype_label, ontology)
    qualifying_relation = term_from_label(
        qualifier + "ComparedTo", ontology)

    match phtype_label:
        case "resistance to chemicals":
            try:
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
            except TypeError:
                "Could not find compound name in atoms."

        case "metal resistance":
            try:
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
            except TypeError:
                "Could not find compound name in atoms."
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
    ref_phtype, ref_trips = find_or_create_node(
        hypo_ds, ref_query_trips, id_prefix="P-REF", namespace=phenotype_namespace
    )
    # Now we are adding the reference state to the states graph
    for t in ref_trips[:4]:
        states.add(t)

    for t in ref_trips[4:]:
        phenotypes.add(t)

    # Fetch the compared phenotype from the graph if it exists,
    # create it if it doesn't exist.
    # TODO: Add Label or other property to make querying easier
    comp_query_trips = (
        [(None, RDFS.subClassOf, phtype)]
        + additional_comp_triples
        + role_between_classes(None, ref_phtype, qualifying_relation)
    )
    comp_phtype, comp_trips = find_or_create_node(
        hypo_ds, comp_query_trips, id_prefix="P", namespace=phenotype_namespace
    )
    for t in comp_trips:
        phenotypes.add(t)

    return (ref_phtype, comp_phtype), [
        t + (states.identifier,) for t in ref_trips[:4]] + [
        t + (phenotypes.identifier,) for t in ref_trips[4:] + comp_trips
    ]


def amino_acid_and_qualifier_to_state(dataset, amino_acid, qualifier, reference_state=DEFAULT_REFERENCE_STATE, phenotype_namespace=HYPO):
    states = dataset.graph("http://hypo.project-genesis.io/states")
    phenotypes = dataset.graph("http://hypo.project-genesis.io/phenotypes")

    qualifier = {"lower": "decreased",
                 "higher": "increased"}.get(qualifier, qualifier)

    qualifying_relation = term_from_label(
        qualifier + "ComparedTo", ontology)

    # Fetch the reference phenotype from the graph if it exists,
    # create it if it doesn't exist. (Link to reference state defined
    # above in variable `REFERENCE_STATE`)
    # TODO: Add Label or other property to make querying easier
    ref_query_trips = (
        role_between_classes(reference_state, None, STATE_HAS_OBSERVABLE)
        + [(None, RDFS.subClassOf, CHEM_COMP_ACC)]
        + role_between_classes(None, amino_acid, CHEM_ACC_OF)
    )
    ref_phtype, ref_trips = find_or_create_node(
        hypo_ds, ref_query_trips, id_prefix="P-REF", namespace=phenotype_namespace
    )
    # Now we are adding the reference state to the states graph
    for t in ref_trips[:4]:
        states.add(t)
    for t in ref_trips[4:]:
        phenotypes.add(t)

    # Fetch the compared phenotype from the graph if it exists,
    # create it if it doesn't exist.
    # TODO: Add Label or other property to make querying easier
    comp_query_trips = (
        [(None, RDFS.subClassOf, CHEM_COMP_ACC)]
        + role_between_classes(None, amino_acid, CHEM_ACC_OF)
        + role_between_classes(None, ref_phtype, qualifying_relation)
    )
    comp_phtype, comp_trips = find_or_create_node(
        hypo_ds, comp_query_trips, id_prefix="P", namespace=phenotype_namespace
    )
    for t in comp_trips:
        phenotypes.add(t)

    # Fetch the complete state from the graph if it exists,
    # create it if it doesn't exist.
    state_query_triples = (
        role_between_classes(None, comp_phtype, STATE_HAS_OBSERVABLE)
        + [(None, RDFS.subClassOf, ORGANISM_STATE)]
    )
    state, state_triples = find_or_create_node(
        hypo_ds, state_query_triples, id_prefix="S", namespace=HYPO
    )
    for t in state_triples:
        states.add(t)

    quads = [
        t + (states.identifier,) for t in ref_trips[:4] + state_triples] + [
        t + (phenotypes.identifier,) for t in ref_trips[4:] + comp_trips
    ]

    # TODO: Add Label or other property to make querying easier
    return state, quads


def hypothesis_to_quads(hypo_ds, h, h_graph):
    # Match the amino acid
    amino_acid = term_from_label(h["observable"], chebi)  # Case-sensitive?
    # print(amino_acid)
    state_1, state_1_quads = amino_acid_and_qualifier_to_state(
        hypo_ds, amino_acid, h["qualifier"], reference_state=DEFAULT_REFERENCE_STATE, phenotype_namespace=HYPO)

    # Represent the logic program
    state_2, state_2_quads = logic_program_to_state(
        h["logic_program"], h_graph)

    # Generate `observable` concepts for the reference
    # triples.append((HYPO.P1, RDF.type, CHEM_COMP_ACC))
    # triples.append((HYPO.P1, CHEM_ACC_OF, amino_acid))
    # triples.extend(role_between_classes())
    # Split the formula

    # State implies State
    quads = state_1_quads + state_2_quads
    quads.extend([t + (h_graph.identifier,)
                 for t in role_between_classes(state_1, state_2, HYPO.implies)])

    return quads


def add_hypothesis_as_new_graph_in_hypo_ds(hypo_ds, h, **kwargs):
    """Add a hypothesis as a new graph in the hypothesis dataset. 
    The graph is created with a unique ID and the hypothesis is added to the metadata graph.

    Args:
        hypo_ds (rdflib.Dataset): The hypothesis dataset.
        h (dict): The hypothesis to add.
        **kwargs: Additional arguments for the hypothesis.


    Returns:
        list: A list of quads added to the dataset.
    """
    quads = []

    hid = create_id(prefix="H")
    h_graphid = rdflib.URIRef(f"http://hypo.project-genesis.io/{hid}")
    h_graph = hypo_ds.graph(h_graphid)
    h_graph.bind("hypo", HYPO)
    h_graph.bind("obo", OBO)
    h_graph.bind("owl", OWL)
    h_graph.bind("rdf", RDF)
    h_graph.bind("rdfs", RDFS)

    quads = hypothesis_to_quads(hypo_ds, h, h_graph)
    for q in quads:
        if q[3] == h_graph.identifier:
            h_graph.add(q[:3])

    # Hack for now to remove empty hypothesis graphs
    # This is not good for future, we want to include
    # all triples related to the hypothesis
    if len(h_graph) == 0:
        print("No triples were added to the graph.")
        hypo_ds.remove_graph(h_graphid)

    # Add the hypothesis to the metadata graph
    meta_trips = []
    meta_trips.append((h_graphid, RDF.type, HYPO.hypothesis))

    # If creation date is provided, add it to the metadata graph
    if "creation_date" in kwargs:
        meta_trips.append(
            (h_graphid, DCT.created, Literal(kwargs["creation_date"])))
    # Else, add the current date
    else:
        meta_trips.append((h_graphid, DCT.created, Literal(
            rdflib.Literal(datetime.datetime.now().isoformat()))))

    # If a creator is provided, add it to the metadata graph
    if "creator" in kwargs:
        meta_trips.append((h_graphid, DCT.creator, Literal(kwargs["creator"])))
    # Else, add 'Genesis' as the creator
    else:
        meta_trips.append((h_graphid, DCT.creator, Literal("Genesis")))

    for t in meta_trips:
        hmeta.add(t)
        quads.append(t + (hmeta.identifier,))

    return quads


def load_hypotheses(hypo_ds, root_dir="experiments"):
    """Searches through the `root_dir` for hypothesis folders and loads them into the dataset."""
    for folder in os.listdir(root_dir):
        try:
            load_hypothesis_from_top_folder(hypo_ds, root_dir, folder)
        except FileNotFoundError:
            print(f"Could not load hypothesis from {folder}")
            # Check next level down in the directory tree
            load_hypotheses(hypo_ds, os.path.join(root_dir, folder))
            continue
        except NotADirectoryError:
            print(f"{folder} is not a directory, ending search.")
            continue


def load_hypothesis_from_top_folder(hypo_ds, root_dir, folder):
    with open(
        os.path.join(root_dir, folder,
                     "hypothesis/selected_hypothesis/hypothesis_details.json"),
        "r",
        encoding="utf-8",
    ) as fi:
        h = json.load(fi)

        # Title?

        # Extract creation date from folder name and load into datetime object
        # (Example folder names: "arginine_202503131539" or "arginine_202503131539 copy")
        DATE_EXTRACTOR = re.compile(r"_(\d{8}\d{4})")
        creation_date_str = DATE_EXTRACTOR.search(folder).group(1)
        creation_date = datetime.datetime.strptime(
            creation_date_str, "%Y%m%d%H%M")

        print(
            f"Loading hypothesis from {os.path.join(root_dir, folder)} (created on {creation_date.strftime('%Y-%m-%d %H:%M:%S')})")
        return add_hypothesis_as_new_graph_in_hypo_ds(
            hypo_ds, h, creation_date=creation_date, creator="Genesis")


if __name__ == "__main__":
    with open("experiments/test-patterns.txt", "r") as fi:
        test_patterns = [l.rstrip().split(':', 1) for l in fi]

    test_hypotheses = [{'logic_program': p[1], 'observable': p[0], 'qualifier': 'higher',
                        'number': i, 'reason': '<blank>'} for (i, p) in enumerate(test_patterns)]

    load_hypotheses(hypo_ds, os.path.join(BASE_DIR, "experiments"))
