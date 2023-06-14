use std::ops::Index;
use std::fs::File;
use std::io::{prelude::*, BufReader};

const NTERM: usize = 256;

pub struct Grammar {
     root: usize,
     left: Vec<usize>,
    right: Vec<usize>,
    sizes: Vec<usize>,
    terminals: Vec<u8>,
}

impl Grammar {
    pub fn from_file(filename: &str) -> Self { 
        let file = File::open(filename).expect("Cannot read grammar file.");
        let reader = BufReader::new(file);

        let mut left = Vec::new();
        let mut right = Vec::new();
        let mut sizes = Vec::new();
        for line in reader.lines() {
            let line = line.expect("Cannot read line in grammar file.");
            let symbols: Vec<usize> = line
                .trim()
                .split_whitespace()
                .map(|x| x.parse().expect("Cannot parse symbol."))
                .collect();

            assert_eq!(symbols.len(), 2, "Incorrect formatting of grammar file.");
            left.push(symbols[0]);
            right.push(symbols[1]);

            let left_size;
            if symbols[0] < NTERM { left_size = 1; }
            else { left_size = sizes[symbols[0] - NTERM] }

            let right_size;
            if symbols[1] < NTERM { right_size = 1; }
            else { right_size = sizes[symbols[1] - NTERM]; }

            sizes.push(left_size + right_size);
        }

        let max_term = (NTERM - 1) as u8;
        let terminals = (0..=max_term).collect();
        Grammar { root: left.len() - 1, left, right, sizes, terminals}
    }

    pub fn len(&self) -> usize { self.sizes[self.root] }

}

impl Index<usize> for Grammar {
    type Output = u8;

    fn index(&self, index: usize) -> &Self::Output {
        let mut symbol = self.root + NTERM;
        let mut skipped = 0;
        while symbol >= NTERM {
            let left_symbol = self.left[symbol - NTERM];
            let left_size;
            if left_symbol < NTERM { left_size = 1; }
            else { left_size = self.sizes[left_symbol - NTERM]; }
            if skipped + left_size > index {
                symbol = left_symbol;
            } else {
                symbol = self.right[symbol - NTERM];
                skipped += left_size;
            }
        }
        return &self.terminals[symbol];
    }
}

