#!/bin/bash
read -r -p "Genomes: " genome_dir
read -r -p "Outputs dir: " out_dir
read -r -p "Keep temps? (y/n): " temp_choice

mkdir "$out_dir"
mkdir tmp_dir

for genome in "$genome_dir"/*.fna;do
      
      genome_name=$(basename "$genome" .fna)
      prodigal -i "$genome" -a "tmp_dir/$genome_name".faa
    done
find tmp_dir -type f -name "*.faa" -exec cat {} + >> "$out_dir/Prodigal_annotation_prot.fasta"

if [ "$temp_choice" == 'n' ];then
    rm -r tmp_dir
fi