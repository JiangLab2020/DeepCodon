
SLIDE=$(grep -E "^SLIDE\s*=\s*True" src/attached_code/config.py)


if [ -n "$SLIDE" ]; then
    OUTPUT_FILE="output_sequences_slide.fasta"
    TSV_NAME="matches_slide.tsv"
else
    OUTPUT_FILE="output_sequences.fasta"
    TSV_NAME="matches.tsv"
fi

echo "Using output file: $OUTPUT_FILE"
echo "Using tsv file: $TSV_NAME"

diamond blastp -d reference -q $OUTPUT_FILE -o $TSV_NAME --top 20
