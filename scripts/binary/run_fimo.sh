#!/bin/bash

O=/g/data/zk16/xzhang/BOM/Tutorial/motifs
M=/g/data/zk16/useful/gimmemotifs/gimme.vertebrate.v5.0.meme
G=/g/data/zk16/cc3704/mouse_data/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa
B=/g/data/zk16/xzhang/BOM/Tutorial/bed_files

while getopts o:g:b:m: op
do 
    case $op in
        o)
            echo "Output file folder is: $OPTARG"
            O=$OPTARG;;
        g)
            echo "Path to genome reference fasta file is: $OPTARG"
            G=$OPTARG;;
        b)
            echo "Path to bed files: $OPTARG"
            B=$OPTARG;;
        m)
            echo "Path to motif database: $OPTARG"
            M=$OPTARG;;
        \?)
            echo "Usage: args [-o] [-g] [-b] [-m] "
            echo "-m means path to motif database"
            echo "-g means path to genome reference fasta file"
            echo "-o set a output path (Default: ./Tutorial/motifs) "
            echo "-b means path to bed file folder (Default: ./Tutorial/bed_files)"
            exit 1;;
    esac
done
  
mkdir -p ${O}
  
for bed_file in ${B}/*
do
    fa_file="${bed_file%.bed}.fa"
    filename=$(basename "$bed_file")
    extension="${filename##*.}"
    filename="${filename%.*}"
    
    # get DNA sequence from bed files
    bedtools getfasta -fi ${G} -bed ${bed_file} -fo ${fa_file}
    
    # use fimo to search motifs
    fimo --thresh 0.0001 --o ${O}/${filename} ${M} ${fa_file}
done
