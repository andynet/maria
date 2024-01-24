use std::ops::Index;
use std::fs::File;
use std::io::{prelude::*, BufReader};
use std::path::Path;

#[cfg(test)]
use std::collections::HashMap;

const NTERM: usize = 256;

pub struct Grammar {
     root: usize,
     left: Vec<usize>,
    right: Vec<usize>,
    sizes: Vec<usize>,
    terminals: Vec<u8>,
}

impl Grammar {
    pub fn from_file<P: AsRef<Path> + ?Sized>(filename: &P) -> Self { 
        let file = File::open(filename).expect("Cannot read grammar file.");
        let reader = BufReader::new(file);

        let mut left = Vec::new();
        let mut right = Vec::new();
        let mut sizes = Vec::new();
        for line in reader.lines() {
            let line = line.expect("Cannot read line in grammar file.");
            let symbols: Vec<usize> = line
                .split_whitespace()
                .map(|x| x.parse().expect("Cannot parse symbol."))
                .collect();

            assert_eq!(symbols.len(), 2, "Incorrect formatting of grammar file.");
            left.push(symbols[0]);
            right.push(symbols[1]);

            let left_size  = if symbols[0] < NTERM { 1 } else { sizes[symbols[0] - NTERM] };
            let right_size = if symbols[1] < NTERM { 1 } else { sizes[symbols[1] - NTERM] };
            sizes.push(left_size + right_size);
        }

        let max_term = (NTERM - 1) as u8;
        let terminals = (0..=max_term).collect();
        Grammar { root: left.len() - 1, left, right, sizes, terminals}
    }

    #[cfg(test)]
    pub fn from_bytes(s: &[u8]) -> Self {
        let s: Vec<_> = s.iter().map(|&x| x as usize).collect();
        let mut rule_id: HashMap<(usize, usize), usize> = HashMap::new();
        let mut id = 256;

        let mut v1 = s;
        let mut v2 = Vec::new();
        while v1.len() > 1 {
            let n = (v1.len() / 2) * 2; // continuing with even lengths only
            for i in (0..n).step_by(2) {
                let rule = rule_id.get(&(v1[i], v1[i+1]));
                match rule {
                    None => { 
                        rule_id.insert((v1[i], v1[i+1]), id);
                        v2.push(id);
                        id += 1;
                    }
                    Some(&id) => {
                        v2.push(id);
                    }
                }
            }
            if v1.len() > n { v2.push(*v1.last().unwrap()) }
            v1 = v2;
            v2 = Vec::new();
        }
        // println!("{:?}", rule_id);
        let mut g: Vec<_> = rule_id.iter().map(|(k, v)| (*v, *k)).collect();
        g.sort_unstable();
        println!("{:?}", g);

        let mut left = Vec::new();
        let mut right = Vec::new();
        let mut sizes = Vec::new();
        for (_, (l, r)) in g {
            left.push(l);
            right.push(r);

            let left_size = if l < NTERM { 1 } else { sizes[l - NTERM] };
            let right_size = if r < NTERM { 1 } else { sizes[r - NTERM] };
            sizes.push(left_size + right_size);
        }

        let max_term = (NTERM - 1) as u8;
        let terminals = (0..=max_term).collect();
        Grammar { root: left.len() - 1, left, right, sizes, terminals }
    }

    #[cfg(test)]
    pub fn print<T: Write>(&self, mut output: T) {
        for i in 0..self.left.len() {
            writeln!(output, "{} {}", self.left[i], self.right[i])
                .expect("Error writting to file.");
        }
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
            let left_size = if left_symbol < NTERM { 1 } else { self.sizes[left_symbol - NTERM] };
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

