#!/bin/bash
set -euo pipefail

zless SARS-CoV2.1k.fa.gz | head > SARS-CoV2.5.fna
bgzip -@ 8 SARS-CoV2.5.fna
samtools faidx SARS-CoV2.5.fna.gz

# PGGB fix
# nvim /home/balaz/.local/share/micromamba/envs/maria/bin/pggb:516
#
# smoothxg_xpoa_cmd="-S"                -> smoothxg_xpoa_cmd=""
# if [[ "$run_abpoa" == true ]]; then
# smoothxg_xpoa_cmd=""                  -> smoothxg_xpoa_cmd="-A"
# fi
#
# https://github.com/pangenome/smoothxg/releases
# Run SPOA by default; the -S flag is replaced by -A to run abPOA

pggb -i SARS-CoV2.5.fna.gz -n 5
mv SARS-CoV2.5.fna.gz.*.final.gfa SARS-CoV2.5.gfa
rm SARS-CoV2.5.fna.gz.fd8c760.*

gunzip SARS-CoV2.5.fna.gz
iss generate --genomes SARS-CoV2.5.fna --n_reads 10 --model HiSeq --output reads
../../tools/moni-0.2.0-Linux/bin/moni build -f -r SARS-CoV2.5.fna -o index
../../tools/moni-0.2.0-Linux/bin/moni mems -i index -p reads_R1.fastq -o mems
cargo run SARS-CoV2.5.gfa > tags.txt

