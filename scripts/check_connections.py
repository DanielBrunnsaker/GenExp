import requests
import sys

try:
    from config import CHEBI_QUERY_ENDPOINT
except ModuleNotFoundError:
    from scripts.config import CHEBI_QUERY_ENDPOINT


# Set up Fuseki endpoint to create new SPARQL store
try:
    print(f"Attempting to connect to ChEBI endpoint at: {CHEBI_QUERY_ENDPOINT}...", end="", file=sys.stderr)
    response = requests.get(CHEBI_QUERY_ENDPOINT, timeout=10)
    response.raise_for_status()
except requests.exceptions.ConnectionError:
    print(f"could not connect.", file=sys.stderr)
    HOST = "localhost"
    CHEBI_FUSEKI_PORT = 3037
    CHEBI_QUERY_ENDPOINT = f"http://{HOST}:{CHEBI_FUSEKI_PORT}/genesis/query"
    print(f"Attempting to connect to ChEBI endpoint at: {CHEBI_QUERY_ENDPOINT}...", end="", file=sys.stderr)
    try:
        response = requests.get(CHEBI_QUERY_ENDPOINT, timeout=10)
        response.raise_for_status()
    except requests.exceptions.RequestException:
        print(f"could not connect.", file=sys.stderr)
        # Perhaps we are in a docker container (ensure that the correct port is forwarded)
        HOST = "host.docker.internal"
        CHEBI_QUERY_ENDPOINT = f"http://{HOST}:{CHEBI_FUSEKI_PORT}/genesis/query"
        print(f"Attempting to connect to ChEBI endpoint at: {CHEBI_QUERY_ENDPOINT}...", end="", file=sys.stderr)
        try:
            response = requests.get(CHEBI_QUERY_ENDPOINT, timeout=10)
            response.raise_for_status()
        except requests.exceptions.RequestException as e:
            # This base class covers ALL requests exceptions
            print(f"""could not connect, out of options.
Please check the `CHEBI_QUERY_ENDPOINT' constant in the `config.py' file
    points to a valid Jena Fuseki endpoint.\n""", file=sys.stderr)
            # print(f"Error: {e}", file=sys.stderr)
            sys.exit(1)  # non-zero exit code signals failure
        else:
            print("Connection successful.", file=sys.stderr)
            sys.exit(0)



# CHEBI_STORE = sparqlstore.SPARQLStore(CHEBI_QUERY_ENDPOINT)
# CHEBI = Graph(CHEBI_STORE, rdflib.graph.DATASET_DEFAULT_GRAPH_ID)