#!/bin/bash

input_file="Enterobacteriaceae_SeqSpe.txt"
output_dir="./../rawdata"
command_file="Fetch_Cds_BaseSpeHttp.pl"
mkdir -p "$output_dir"


split -d -l $(( $(wc -l < "$input_file") / 10 + 1 )) "$input_file" "$output_dir/split_"


split_files=$(ls "$output_dir"/split_*)


count=0
for file in $split_files; do
    mkdir -p "$output_dir/$count"
    mv "$file" "$output_dir/$count/$(basename $file)"
    cp "$command_file" "$output_dir/$count/$command_file"
    ((count++))
done

