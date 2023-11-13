fasta_base := "data/real/SARS-CoV2.5"
reads_base := "data/real/reads_R1"

build_graph_from_fasta:
    bgzip -@ 16 -c {{fasta_base}}.fna > {{fasta_base}}.fna.gz
    samtools faidx {{fasta_base}}.fna.gz
    pggb -i {{fasta_base}}.fna.gz -n 5 -d tmp -o {{fasta_base}}
    sort {{fasta_base}}/*.smooth.final.gfa > {{fasta_base}}.gfa

build_fasta_from_graph:
    echo "skdgj"

run_moni:
    tools/moni-0.2.0-Linux/bin/moni build -f -r {{fasta_base}}.fna -o {{fasta_base}}
    tools/moni-0.2.0-Linux/bin/moni ms -i {{fasta_base}} -p {{reads_base}}.fastq -o {{reads_base}}
    tools/moni-0.2.0-Linux/bin/moni mems -i {{fasta_base}} -p {{reads_base}}.fastq -o {{reads_base}}

make_grammar:
    # mamba activate maria
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
    gfatk path --all $gfa > tmp-genomes.fa
    bwa index tmp-genomes.fa
    bwa fastmap tmp-genomes $reads -l 2 > fastmap.out
    python fastmap2sam.py
    samtools view -b tmp-sam -o tmp-bam
    gfainject --gfa tmp-sorted-graph.gfa --bam tmp-bam  > results-mem.gaf
#
# cargo run --bin main -- -g data/real/SARS-CoV2.5.gfa -m data/real/reads_R1.mems -p data/real/reads_R1.pointers

