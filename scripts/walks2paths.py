import sys
import re


def to_pathnode(walknode):
    number = walknode[1:]
    if walknode[0] == ">":
        return number + "+"
    return number + "-"


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: walks2paths.py <gfa_with_walks> > <output_gfa>")
        exit()

    gfa_filename = sys.argv[1]
    with open(gfa_filename) as f:
        lines = f.readlines()

    for line in lines:
        if line.startswith("H"):
            print(line, end="")
        if line.startswith("S"):
            print(line, end="")
        if line.startswith("L"):
            print(line, end="")
        if line.startswith("W"):
            # RecordType SampleId HapIndex SeqId SeqStart SeqEnd Walk
            (_, sample, _, _, _, _, walk) = line.split()
            walk = re.findall("[><][0123456789]+", walk)
            walk = list(map(to_pathnode, walk))
            walk = ",".join(walk)
            print(f"P\t{sample}\t{walk}\t*")
