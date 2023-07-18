#!/bin/bash

output_dir=$1
shift

genome_path=$1
shift

mkdir -p "$output_dir"
# Set path to motif database
motifsDB='/path/to/gimme.vertebrate.v5.0.meme'

for bed_file in "$@"
do
    fa_file="${bed_file%.bed}.fa"
    bedtools getfasta -fi ${genome_path} -bed ${bed_file} -fo ${fa_file}
    filename=$(basename "$bed_file")
    extension="${filename##*.}"
    filename="${filename%.*}"

    fimo --thresh 0.0001 --o "$output_dir/$filename" "$motifsDB" "$fa_file"
done

