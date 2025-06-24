EXPERIMENT_DIRECTORY="experiments"

find $EXPERIMENT_DIRECTORY -name "study.trig" -delete
find $EXPERIMENT_DIRECTORY -name "hypotheses.trig" -delete
find $EXPERIMENT_DIRECTORY -name "merged_dataset.trig" -delete

echo "Creating hypothesis database from \`./$EXPERIMENT_DIRECTORY\` directory..." >&2
HYPO_DS_TRIG=$(python scripts/hypothesis_formalisation.py)
ls $HYPO_DS_TRIG >/dev/null &&
    echo "Hypotheses stored in: "$HYPO_DS_TRIG" ("$(ls -lh $HYPO_DS_TRIG | awk '{print $5}')"B)" >&2
python scripts/collect-db-files.py
