from Bio import SeqIO
import sys


def print_header(fasta: str):
    for record in SeqIO.parse(fasta, "fasta"):
        print(f"@SQ\tSN:{record.id}\tLN:{len(record.seq)}")


if __name__ == "__main__":
    fasta = sys.argv[1]
    mapping = sys.argv[2]

    print_header(fasta)

    # QNAME FLAG RNAME POS MAPQ CIGAR RNEXT PNEXT TLEN SEQ QUAL TAGS
    # QNAME:    ENA|MW565758|MW565758.1_0_1/1
    # FLAG:     16  // 0 for standard strand, 16 for reverse strand
    # RNAME:    ENA|MW565758|MW565758.1
    # POS:      21551   // 1-based
    # MAPQ:     60
    # CIGAR:    126M
    # RNEXT:    *
    # PNEXT:    0
    # TLEN:     0
    # SEQ:      TCTTATTTATAGATGGTATCTACCTGAAAGATGTCATTATTTTTCTGAAATCATACTGTTGTTTTCCAACAGATATCTACCCAAGGCCAAAGCCACGCCCTCAACCCCAGCCTGGCAATTCCGGCAACAGTGGAGGTAATGAGTATTTATT
    # QUAL:     FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF,FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF
    # TAGS:     NM:i:0	MD:Z:151	AS:i:151	XS:i:0
    with open(mapping) as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("SQ"):
            (_, qname, qlen) = line.split()
        if line.startswith("EM"):
            (_, qstart, qend, n, rest) = line.split(maxsplit=4)
            for rrecord in rest.split():
                (rname, rest) = rrecord.split(':')
                flag = 0 if rest.startswith("+") else 16
                pos = rest[1:]

                mapq = 60
                tlen = int(qend) - int(qstart)
                cigar = f"{tlen}M"
                rnext = "*"
                pnext = "0"
                tags = ""
                seq = "*"
                qual = "*"
                print(
                    f"{qname}\t{flag}\t{rname}\t{pos}\t{mapq}\t{cigar}"
                    f"\t{rnext}\t{pnext}\t{tlen}\t{seq}\t{qual}\t{tags}"
                )
