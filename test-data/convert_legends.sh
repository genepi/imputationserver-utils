#!/bin/bash

# Check if the correct number of arguments is provided
if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <input_file> <output_file> <chr>"
    exit 1
fi

# Read input and output file names from arguments
INPUT_FILE="$1"
OUTPUT_FILE="$2"
TEMP_FILE="${OUTPUT_FILE%.gz}.new"
CHR="$3"

# Check if the input file exists
if [ ! -f "$INPUT_FILE" ]; then
    echo "Input file $INPUT_FILE not found!"
    exit 1
fi

# Step 1: Use csvtk to modify the file
echo "Running csvtk mutate2..."
csvtk mutate2 -d " " --after "id" -s -n chr -e "'$CHR'" -T "$INPUT_FILE" > "$TEMP_FILE"

# Step 2: Compress the new file with bgzip
echo "Running bgzip..."
bgzip "$TEMP_FILE"

# Rename the compressed file to the desired output name
mv "$TEMP_FILE.gz" "$OUTPUT_FILE"

# Step 3: Index the compressed file with tabix
echo "Running tabix..."
tabix -s 2 -b 3 -e 3 -S 1 "$OUTPUT_FILE" -f

# Clean up temporary file
rm "$TEMP_FILE"

echo "Processing complete. Output file: $OUTPUT_FILE"
