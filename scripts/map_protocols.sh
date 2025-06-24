#!/bin/bash
export RMLMAPPER_JAR=/opt/tools/rmlmapper-7.3.3-r374-all.jar
# place first argument in SOURCE_JSON
export ROOT_DIR="$1"
export LAST_DIR=$(basename "$ROOT_DIR")

export HYPOTHESIS_JSON="tmp/hypothesis_details.json"
export PROTOCOL_JSON="tmp/protocol.json"
export GROWTH_DATA_AUC_FOREST="tmp/data.csv"

# Preprocessing the json files
python scripts/rml-preprocessing/protocol_preprocessing.py \
    $ROOT_DIR"/protocol/protocol.json" >$PROTOCOL_JSON".tmp"

python scripts/rml-preprocessing/hypothesis_preprocessing.py \
    $ROOT_DIR"/hypothesis/selected_hypothesis" >$HYPOTHESIS_JSON".tmp"

# Preprocessing the growth csv file
python scripts/rml-preprocessing/growth_preprocessing.py \
    $ROOT_DIR"/results/growth/tests/AUC_forest.csv" >$GROWTH_DATA_AUC_FOREST

jq 'with_entries(select(.value != null))' \
    $HYPOTHESIS_JSON".tmp" \
    >$HYPOTHESIS_JSON
jq '.experiments |= map(with_entries(select(.value != null)))' \
    $PROTOCOL_JSON".tmp" \
    >$PROTOCOL_JSON
# jq 'with_entries(select(.value != null))' \
#     $ROOT_DIR"/results/growth/.json" \
#     >$GROWTH_DATA_JSON
# echo $LAST_DIR

# Get trig file name based on the input file name SOURCE_JSON (keep directory structure)
OUTPUT_FILE=$ROOT_DIR"/protocol/study.trig"
mkdir -p "$(dirname "$OUTPUT_FILE")"
# Fill in the template
envsubst <ontology-files/rml-maps/protocol-mapper.ttl >tmp/protocol-mapper.ttl

# # Show what was written
# echo "Generated tmp/protocol-mapper.ttl:"
# cat tmp/protocol-mapper.ttl

# Run the RML mapper
# replacing the input file from SOURCE_JSON and outputting to .trig format
java -jar $RMLMAPPER_JAR -m tmp/protocol-mapper.ttl \
    -o $OUTPUT_FILE \
    -s trig

# rm tmp/*
