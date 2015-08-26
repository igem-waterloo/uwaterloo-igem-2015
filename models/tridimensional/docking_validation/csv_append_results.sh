#!/bin/bash
# Accept a single directory as an argument, find all CSVs in that directory (& sub-directories) and append them
# to a single file in the directory with an additional column of the original filename

csv_dir=$1
if [ -z "$1" ]; then
    echo "ERROR: requires csv_dir argument"
    echo "Example: csv_append_results.sh /results/csvs_July30"
    exit 1
fi

outcsv=$csv_dir"/"$(basename "$csv_dir").csv
if [ -e $outcsv ]; then
    echo "ERROR: Output file $outcsv already exists!"
    exit 1
else
    echo "Appending all CSVs found in $csv_dir"
    echo "to $outcsv..."
fi

# Find constants.py, assuming it is stored in the same directory as this script
docking_validation_dir=$( cd "$(dirname "${BASH_SOURCE[0]}" ) " && pwd)
constants="$(dirname $docking_validation_dir)/constants.py"

# Combine the line that starts "CSV HEADER" & the subsequent line, keep only what is between the square brackets, then
# remove all the quote and whitespace characters -> gives CSV header as a single CSV line
csv_cols=$(grep -A 1 "CSV_HEADER" $constants | tr '\n' ' '  | sed 's/.*\[\(.*\)\].*/\1/' | sed $'s/\'//g' | sed -e 's/\s//g')

# Initialize output CSV
echo $outcsv
echo "$csv_cols,run_folder" > "$outcsv"

# Find all CSV in supplied directory
find "$csv_dir" -name '*.csv' -print0 | while read -d '' -r csv_file; do
    # Assume the dir containing the CSVs has a decsriptive name
    label=$(basename "$(dirname "$csv_file")")
    echo "Appending from $label..."

    # Pass contents of CSV file except header to awk, which adds column to the end containing label variable
    cat "$csv_file" | tail -n+2 | awk -v label="$label" -F, '{$(NF+1)=label;}1' OFS=, >> "$outcsv"
done

echo "File $outcsv finished!"
