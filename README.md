https://github.com/plasma-umass/scalene

# MONI installation
```
cd tools
wget https://github.com/maxrossi91/moni/releases/download/v0.2.0/moni-0.2.0-Linux.tar.gz
tar -xvzf moni-0.2.0-Linux.tar.gz
./moni-0.2.0-Linux/bin/moni -h
```

# Grammar for LCE queries

## How to get it running
```
git clone https://gitlab.com/manzai/bigrepair.git
cd bigrepair
make
```

We also need the helper script `scripts/print_plain_slp.c` compiled, we can do it with
```
cd scripts
cc -Wall -g -O9    print_plain_slp.c   -o print_plain_slp
```

## Usage
```
./tools/bigrepair/bigrepair data/grammar/tiny
./scripts/print_plain_slp data/grammar/tiny
```
There are 256 terminals in the grammar representing the char values 0-255. The terminal for the char `x` is associated with the identifier of `x`.
Since not all of the terminals are printable and may cause formatting problems (e.g. `\n`), they are ommited from the plaintext output.
The rest of the grammar, i.e. non-terminals and associated rules, is contained in the file `tiny.plainslp`.
The single production rule for every non-terminal is defined in a separate line numbered `L` (starting from 1), the id of the non-terminal is `id = L+255`. The rule is in the form `X Y`, where `X < id` and `Y < id`.
Hence the first line of `tiny.plainslp` represents the rule `256` and so on.
