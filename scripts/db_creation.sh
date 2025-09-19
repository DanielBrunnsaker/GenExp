EXPERIMENT_DIRECTORY="experiments"
TMP_DELETED_DIR=$(mktemp -d)
mkdir -p tmp

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
        # Move files to temp directory instead of deleting
        while IFS= read -r file; do
            mv "$file" "$TMP_DELETED_DIR/"
        done <<<"$files_found"
        echo "Files moved to temporary directory: $TMP_DELETED_DIR"
    else
        echo "Database creation aborted. No files were deleted."
        exit 0
    fi
else
    echo "No matching files found. Nothing to delete."
fi

echo "Creating hypothesis database from \`./$EXPERIMENT_DIRECTORY\` directory..." >&2
python scripts/hypothesis_formalisation.py
status=$?
HYPO_DS_TRIG=$(cat tmp/ds_filename.txt)
if [ $status -ne 0 ]; then
    echo "Error: Database creation failed. Restoring deleted files..." >&2
    # Move files back to their original locations
    for file in "$TMP_DELETED_DIR"/*; do
        orig_path="$EXPERIMENT_DIRECTORY/$(basename "$file")"
        mv "$file" "$orig_path"
    done
    rm -rf "$TMP_DELETED_DIR"
    exit 1
else
    ls $HYPO_DS_TRIG >/dev/null &&
        echo "Hypotheses stored in: "$HYPO_DS_TRIG" ("$(ls -lh $HYPO_DS_TRIG | awk '{print $5}')"B)" >&2
    python scripts/collect-db-files.py
    # If successful, remove the temp directory and its contents
    rm -rf "$TMP_DELETED_DIR"
fi
