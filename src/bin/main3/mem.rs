use std::fs::File;
use std::io::{BufReader, BufRead, Lines};

pub struct MEMReader {
    mem_lines: Lines<BufReader<File>>,
    ptr_lines: Lines<BufReader<File>>,
}

impl MEMReader {
    pub fn new(mems_filename: &str, ptrs_filename: &str) -> Self {
        let mem_lines = BufReader::new(File::open(mems_filename).expect("Cannot open MEM file.")).lines();
        let ptr_lines = BufReader::new(File::open(ptrs_filename).expect("Cannot open PTR file.")).lines();

        Self{mem_lines, ptr_lines}
    }
}

impl Iterator for MEMReader {
    type Item = (String, Vec<(usize, usize, usize)>);

    fn next(&mut self) -> Option<Self::Item> {
        let id1 = self.mem_lines.next();
        let id2 = self.ptr_lines.next();
        if id1.is_none() || id2.is_none() { return None; }

        let id1 = id1.unwrap().expect("Error reading line in MEM.");
        let id2 = id2.unwrap().expect("Error reading line in PTR.");
        if id1 != id2 { panic!("IDs are not the same!\n{id1}\n{id2}") }
        let id = id1.strip_prefix('>').expect("Id does not start with >").to_owned();

        let mem = self.mem_lines.next()
            .expect("Expected mem line.")
            .expect("Error reading line in MEM.");
        let ptr = self.ptr_lines.next()
            .expect("Expected ptr line.")
            .expect("Error reading line in PTR.");

        let mems: Vec<(usize, usize)> = mem.split_whitespace()
            .map(parse_MEM).collect();
        let ptrs: Vec<usize> = ptr.split_whitespace()
            .map(|x| x.parse().expect("Cannot parse PTR.")).collect();

        let mut resulting_mems = Vec::new();
        for mem in mems {
            // length, read position, reference position
            resulting_mems.push((mem.1, mem.0, ptrs[mem.0]));
        }
        return Some((id, resulting_mems));
    }
}

#[allow(non_snake_case)]
fn parse_MEM(mem: &str) -> (usize, usize) {
    let mem = &mem[1..mem.len()-1];
    let v: Vec<usize> = mem.split(',')
        .map(|x| x.parse().expect("Cannot parse MEM."))
        .collect();
    return (v[0], v[1]);
}

