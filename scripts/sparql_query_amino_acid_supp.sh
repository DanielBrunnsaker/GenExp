# Ask user what amino acid they want to search for
# and take input as a ChEBI ID
echo "Enter the ChEBI ID of the supplemented amino acid you want to search for (e.g., 17203 for Proline):"
read AMINO_ACID_CHEBI

# Replace colon with underscore for SPARQL query
export AA_OBO_URI="obo:CHEBI_"$AMINO_ACID_CHEBI
echo "Amino acid supplement OBO URI: $AA_OBO_URI"

# Do the same but for treatment chemical
echo "Enter the ChEBI ID of the treatment chemical you want to search for (e.g., 30751 for Formic Acid):"
read SUPP_CHEBI

export SUPP_OBO_URI="obo:CHEBI_"$SUPP_CHEBI
echo "Treatment chemical OBO URI: $SUPP_OBO_URI"

# Define the SPARQL query with user input
SPARQL_QUERY_ASK="$(envsubst <scripts/sparql-templates/ask_aa_supp_template.rq)"
SPARQL_QUERY_SELECT="$(envsubst <scripts/sparql-templates/select_aa_supp_template.rq)"

# Ask first
echo "Running ASK query to check for experiments with $AMINO_ACID_CHEBI and $SUPP_CHEBI..."
ASK_RESULT=$(echo $SPARQL_QUERY_ASK|python scripts/pyoxi_query.py)
echo -e "ASK query result:\\n$ASK_RESULT"


# If true, run the SELECT query (note need to trim whitespace frmo the $ASK_RESULT variable)
if [ "$(echo $ASK_RESULT | xargs)" = "True" ]; then
    echo "Experiments found. Running SELECT query..."
    SELECT_RESULT=$(echo $SPARQL_QUERY_SELECT|python scripts/pyoxi_query.py)
    echo -e "SELECT query result:\\n$SELECT_RESULT"
else
    echo "No experiments found with the supplemented amino acid and treatment chemical."
fi