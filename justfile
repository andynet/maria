fasta_base := "data/2023-11-13_example_run/SARS-CoV2.5"
reads_base := "data/2023-11-13_example_run/reads_R1"
# mamba activate maria

revcomp:
    seqkit seq -i {{fasta_base}}.fna > tmp.fna
    seqtk seq -r tmp.fna | sed "/^>/ s/$/_rev/" > tmp.rev.fna
    cat tmp.fna tmp.rev.fna | seqkit sort -w 0 > {{fasta_base}}.fna

build_graph_from_fasta:
    bgzip -@ 16 -c {{fasta_base}}.fna > {{fasta_base}}.fna.gz
    samtools faidx {{fasta_base}}.fna.gz
    pggb -i {{fasta_base}}.fna.gz -n 5 -d tmp -o {{fasta_base}}
    sort -k2 -n {{fasta_base}}/*.smooth.final.gfa > {{fasta_base}}.gfa

# build_fasta_from_graph:
#     gfatk path --all {{fasta_base}}.gfa > {{fasta_base}}.fna

# add_revcomp:
#     seqtk seq -r {{fasta_base}}.fna | sed "/^>/ s/$/_rev/" > {{fasta_base}}.rev.fna
#     cat {{fasta_base}}.fna {{fasta_base}}.rev.fna > tmp.fna
#     mv tmp.fna {{fasta_base}}.fna

run_moni:
    tools/moni-0.2.0-Linux/bin/moni build -f -r {{fasta_base}}.fna -o {{fasta_base}}
    tools/moni-0.2.0-Linux/bin/moni ms -i {{fasta_base}} -p {{reads_base}}.fastq -o {{reads_base}}
    tools/moni-0.2.0-Linux/bin/moni mems -i {{fasta_base}} -p {{reads_base}}.fastq -o {{reads_base}}

make_grammar:
    less {{fasta_base}}.fna | grep -v "^>" | tr -d "\n" > {{fasta_base}}.fnajoin
    ./tools/bigrepair/bigrepair {{fasta_base}}.fnajoin
    ./scripts/print_plain_slp {{fasta_base}}.fnajoin

run_maria:
    cargo run --bin main3 -- \
        -t data/pftag/triggers.txt \
        -f {{fasta_base}}.gfa \
        -g {{fasta_base}}.fnajoin.plainslp \
        -m {{reads_base}}.mems \
        -p {{reads_base}}.pointers

run_alternative:
    bwa index {{fasta_base}}.fna
    bwa fastmap {{fasta_base}}.fna {{reads_base}}.fastq -l 2 > {{reads_base}}.fastmap
    python ./scripts/fastmap2sam.py {{fasta_base}}.fna {{reads_base}}.fastmap > {{reads_base}}.sam
    samtools view {{reads_base}}.sam -o {{reads_base}}.bam -b
    ./tools/gfainject --gfa {{fasta_base}}.gfa --bam {{reads_base}}.bam  > {{reads_base}}_alternative.gaf
#
# cargo run --bin main -- -g data/real/SARS-CoV2.5.gfa -m data/real/reads_R1.mems -p data/real/reads_R1.pointers

