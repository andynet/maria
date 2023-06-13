// #![allow(non_snake_case)]

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
        let id = id1.strip_prefix(">").expect("Id does not start with >").to_owned();

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
    let v: Vec<usize> = mem.split(",")
        .map(|x| x.parse().expect("Cannot parse MEM."))
        .collect();
    return (v[0], v[1]);
}

/*
let mem_file = File::open(&args.mems_filename).expect("Cannot open MEM file.");
let ptr_file = File::open(&args.ptr_filename).expect("Cannot open PTR file.");

let mut mem_reader = BufReader::new(mem_file);
let mut ptr_reader = BufReader::new(ptr_file);

let mut mem_line = String::new();
let mut ptr_line = String::new();

let mut b1 = mem_reader.read_line(&mut mem_line).expect("Cannot read mem line.");
let mut b2 = ptr_reader.read_line(&mut ptr_line).expect("Cannot read ptr line.");

let mut read_id = String::new();
while b1 != 0 && b2 != 0 {
    if !mem_line.starts_with(">") {
        let mems: Vec<MEM> = mem_line.split_whitespace()
            .map(|x| x.parse().expect("Cannot parse MEM")).collect();
        let ptrs: Vec<usize> = ptr_line.split_whitespace()
            .map(|x| x.parse().expect("Cannot parse PTR")).collect();

        for mem in &mems {
            // ads +1 for every $ inserted to seq
            let mut adj = 0;
            while lengths[adj] < ptrs[mem.0] { adj += 1; }
            let ref_mem = (ptrs[mem.0]+adj, mem.1);
            let ref_pos;
            if adj != 0 { ref_pos = ptrs[mem.0]+adj-lengths[adj-1]; }
            else { ref_pos = ptrs[mem.0]; }


            

        }
    } else {
        read_id = mem_line.clone();
        read_id = read_id.strip_prefix(">").unwrap().trim().to_owned();
    }

    mem_line.clear();
    ptr_line.clear();

    b1 = mem_reader.read_line(&mut mem_line).expect("Cannot read mem line.");
    b2 = ptr_reader.read_line(&mut ptr_line).expect("Cannot read ptr line.");
}
*/
