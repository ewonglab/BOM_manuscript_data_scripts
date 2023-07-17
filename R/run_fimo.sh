#!/bin/bash

output_dir=$1
shift

genome_path=$1
shift

mkdir -p "$output_dir"
# Set path to motif database
motifsDB='/g/data/zk16/useful/gimmemotifs/gimme.vertebrate.v5.0.meme'
genome_path=/g/data/zk16/cc3704/mouse_data/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa

for bed_file in "$@"
do
    fa_file="${bed_file%.bed}.fa"
    bedtools getfasta -fi ${genome_path} -bed ${bed_file} -fo ${fa_file}
    filename=$(basename "$bed_file")
    extension="${filename##*.}"
    filename="${filename%.*}"

    /g/data/zk16/cc3704/meme/meme-5.4.0/src/fimo --thresh 0.0001 \
    --o "$output_dir/$filename" "$motifsDB" "$fa_file"
done

