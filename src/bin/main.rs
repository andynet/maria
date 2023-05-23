use gfa::parser::GFAParser;
use gfa::gfa::{Path, SegmentId};
use gfa::optfields::*;
use maria::{create_map,create_sequence,inverse_suffix_array,tag_array};
use bio::data_structures::suffix_array::suffix_array;
use std::str;
use std::str::FromStr;
use clap::Parser;
use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::collections::{HashMap, HashSet};
// use std::io::prelude::*;

/// Find MEMs in a graph
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// GFA file
    #[arg(short)]
    gfa_filename: String,
    /// MEMs file
    #[arg(short)]
    mems_filename: String,
    /// pointers
    #[arg(short)]
    ptr_filename: String,
}

#[derive(Debug, PartialEq, Eq, Hash, Clone)]
struct GraphPos {
    node_id: Vec<u8>,
    pos: usize,
}

impl fmt::Display for GraphPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", str::from_utf8(&self.node_id).unwrap(), self.pos)
    }
}

impl From<(Vec<u8>, usize)> for GraphPos {
    fn from(value: (Vec<u8>, usize)) -> Self {
        GraphPos{ node_id: value.0, pos: value.1 }
    }
}

#[derive(Debug, PartialEq, Eq)]
struct ParseMEMError;

#[derive(Debug)]
struct MEM(usize, usize);

impl FromStr for MEM {
    type Err = ParseMEMError;

    fn from_str(value: &str) -> Result<Self, Self::Err> {
        let value = &value[1..value.len()-1];
        let v: Vec<usize> = value.split(",").map(|x| x.parse().expect("Error.")).collect();
        Ok(Self(v[0], v[1]))
    }
}

fn print_tag(tag: &[GraphPos]) {
    print!("{}: ", tag.len());
    for i in 0..tag.len() { print!("{} ", tag[i]); }
    println!();
}

fn get_sampled_arrays<N, T>(
    seq: &[u8], paths: &[Path<N, T>], map: &HashMap<Vec<u8>, Vec<u8>>
) -> (Vec<GraphPos>, Vec<usize>) 
where
    N: SegmentId + Clone,
    T: OptFields + Clone
{
    let sa  = suffix_array(&seq);
    let isa = inverse_suffix_array(&sa);
    let tag = tag_array(&sa, &isa, &paths, &map);

    let mut tag_sample = Vec::new();
    let mut sa_sample = Vec::new();
    tag_sample.push(tag[0].clone().into());
    sa_sample.push(sa[0]);
    for i in 1..tag.len() {
        if tag[i] != tag[i-1] {
            tag_sample.push(tag[i-1].clone().into());
            sa_sample.push(sa[i-1]);
            tag_sample.push(tag[i].clone().into());
            sa_sample.push(sa[i])
        }
    }
    tag_sample.push(tag[tag.len()-1].clone().into());
    sa_sample.push(sa[tag.len()-1]);
    return (tag_sample, sa_sample);
}

fn get_mems(mems: &str, ptrs: &str) -> Vec<(usize, usize)> {
    let mut res = Vec::new();

    let mem_file = File::open(mems).expect("Cannot open MEM file.");
    let ptr_file = File::open(ptrs).expect("Cannot open PTR file.");

    let mut mem_reader = BufReader::new(mem_file);
    let mut ptr_reader = BufReader::new(ptr_file);

    let mut mem_line = String::new();
    let mut ptr_line = String::new();

    let mut b1 = mem_reader.read_line(&mut mem_line).expect("Cannot read mem line.");
    let mut b2 = ptr_reader.read_line(&mut ptr_line).expect("Cannot read ptr line.");

    while b1 != 0 && b2 != 0 {
        if !mem_line.starts_with(">") {
            let mems: Vec<MEM> = mem_line.split_whitespace()
                .map(|x| x.parse().expect("Cannot parse MEM")).collect();
            let ptrs: Vec<usize> = ptr_line.split_whitespace()
                .map(|x| x.parse().expect("Cannot parse PTR")).collect();

            for mem in &mems {
                // because corona, TODO: fix this
                // ads +1 for every $ inserted to seq
                let adj = ptrs[mem.0] / 29850;
                res.push((ptrs[mem.0]+adj, mem.1));

            }
        }

        mem_line.clear();
        ptr_line.clear();

        b1 = mem_reader.read_line(&mut mem_line).expect("Cannot read mem line.");
        b2 = ptr_reader.read_line(&mut ptr_line).expect("Cannot read ptr line.");
    }
    return res;
}

fn list_unique(tag: &[GraphPos]) -> Vec<GraphPos> {
    let mut set: HashSet<GraphPos> = HashSet::new();
    for i in 0..tag.len() { set.insert(tag[i].clone()); }
    let mut res = Vec::new();
    for item in set { res.push(item); }
    return res;
}

/// returns (l, f) such that: 
/// seq[s1..s1+l] == seq[s2..s2+l]
/// if seq[s1+l] < seq[s2+l]: f = True
fn lce(seq: &[u8], s1: usize, s2: usize) -> (usize, bool) {
    let n = seq.len();
    if s1 == s2 { return (n-s1, false); }    // the same is not smaller

    let mut l = 0;
    while s1+l < n && s2+l < n && seq[s1+l] == seq[s2+l] { l += 1; }
    if s1+l == n { return (l, true); }
    if s2+l == n { return (l, false); }
    if seq[s1+l] < seq[s2+l] { return (l, true); }
    if seq[s1+l] > seq[s2+l] { return (l, false); }
    panic!();
}

fn get_lower(seq: &[u8], mem: &(usize, usize), sa: &[usize]) -> usize {
    let mut l = 0;
    let mut r = sa.len();

    while l < r-1 {
        let m = (l + r) / 2;
        let (e, sa_smaller) = lce(seq, sa[m], mem.0);
        if e < mem.1 && sa_smaller { l = m; }
        else { r = m; }
    }
    return r;
}

fn get_upper(seq: &[u8], mem: &(usize, usize), sa: &[usize]) -> usize {
    let mut l = 0;
    let mut r = sa.len();

    while l < r-1 {
        let m = (l + r) / 2;
        let (e, sa_smaller) = lce(seq, sa[m], mem.0);
        if e < mem.1 && !sa_smaller { r = m; }
        else { l = m; }
    }

    return r;
}

fn get_graph_positions(seq: &[u8], mem: &(usize, usize), tag: &[GraphPos], sa: &[usize]) -> Vec<GraphPos> {
    let lower = get_lower(seq, mem, sa);          // included
    let upper = get_upper(seq, mem, sa);          // excluded

    return list_unique(&tag[lower..upper]);
}

fn main() {
    let args = Args::parse();

    let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();
    let gfa = parser.parse_file(&args.gfa_filename).expect("Error parsing file.");

    let map = create_map(&gfa.segments);
    let seq = create_sequence(&gfa.paths, &map); 
    // ^- This is a problem, because it is linear size.
    // How to get rid of it?

    let (stag, ssa) = get_sampled_arrays(&seq, &gfa.paths, &map);
    print_tag(&stag[..10]);
    println!("{:?}", &ssa[..10]);

    let mems = get_mems(&args.mems_filename, &args.ptr_filename);
    println!("{:?}", mems);

    for mem in &mems {
        let gp: Vec<GraphPos> = get_graph_positions(&seq, &mem, &stag, &ssa);
        println!("{:?}: {:?}", mem, gp);
    }
}


#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_lce() {
        let seq = b"AAAAAAAAAA$";
        assert_eq!(lce(seq, 0, 5), (5, false));
        let seq = b"AAAA$AAAA$";
        assert_eq!(lce(seq, 0, 5), (5, false));
        let seq = b"AGCTGCTGCTTGATGCTGATCG$";
        assert_eq!(lce(seq, 1, 4), (6, true));
        assert_eq!(lce(seq, 0, 0), (23, false));
    }

    #[test]
    fn test_binsearch() {
        let seq = b"AGGTTAGTAC$AGTAACGTTAAC$";
        let sa = suffix_array(seq);

        let mem = (17, 2);
        let expected_result = (14, 18);
        let lower = get_lower(seq, &mem, &sa);
        let upper = get_upper(seq, &mem, &sa);
        assert_eq!((lower, upper), expected_result);

        let mem = (11, 2);
        let expected_result = (7, 10);
        let lower = get_lower(seq, &mem, &sa);
        let upper = get_upper(seq, &mem, &sa);
        assert_eq!((lower, upper), expected_result);
    }
}

