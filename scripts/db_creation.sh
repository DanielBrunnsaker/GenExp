EXPERIMENT_DIRECTORY="experiments"
# Check connections before proceeding
if ! python scripts/check_connections.py; then
    echo "check_connections.py failed. Exiting."
    exit 1
fi

echo "check_connections.py succeeded. Continuing..."

# Look for the files first (without deleting)
files_found=$(find "$EXPERIMENT_DIRECTORY" -type f \( -name "study.trig" -o -name "hypotheses.trig" -o -name "merged_dataset.trig" \))

if [ -n "$files_found" ]; then
    echo "The following files were found:"
    while IFS= read -r file; do
        echo -e "\t- $file"
    done <<<"$files_found"
    read -p "Do you want to delete these files? (y/N): " confirm
    if [[ "$confirm" =~ ^[Yy]$ ]]; then
        # Proceed with deletion
        find "$EXPERIMENT_DIRECTORY" -type f \( -name "study.trig" -o -name "hypotheses.trig" -o -name "merged_dataset.trig" \) -delete
        echo "Files deleted."
    else
        echo "Database creation aborted. No files were deleted."
        exit 0
    fi
else
    echo "No matching files found. Nothing to delete."
fi

echo "Creating hypothesis database from \`./$EXPERIMENT_DIRECTORY\` directory..." >&2
HYPO_DS_TRIG=$(python scripts/hypothesis_formalisation.py)
ls $HYPO_DS_TRIG >/dev/null &&
    echo "Hypotheses stored in: "$HYPO_DS_TRIG" ("$(ls -lh $HYPO_DS_TRIG | awk '{print $5}')"B)" >&2
python scripts/collect-db-files.py
