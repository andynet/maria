use gfa::parser::GFAParser;
use gfa::gfa::Segment;
use gfa::gfa::{Path, SegmentId};
use gfa::optfields::*;
use maria::inverse_suffix_array;
use bio::data_structures::suffix_array::suffix_array;
use bio::data_structures::suffix_array::lcp as lcp_array;
use std::str;
use std::str::FromStr;
use clap::Parser;
use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::io::prelude::*;
use std::collections::{HashMap, HashSet};

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

#[derive(PartialEq, Eq, Hash, Clone)]
struct GraphPos {
    node_id: Vec<u8>,
    pos: usize,
}

impl fmt::Display for GraphPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", str::from_utf8(&self.node_id).unwrap(), self.pos)
    }
}

impl fmt::Debug for GraphPos {
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

fn create_sequence(paths: &[Path<usize, ()>], map: &HashMap<usize, Vec<u8>>) 
    -> (Vec<u8>, Vec<usize>)
{
    let mut seq = Vec::new();
    let mut lengths = Vec::new();
    for p in paths {
        todo!();
//        let segments = parse_segments(&p.segment_names);
//        let sequence = get_sequence(&segments, &map);

//         seq.extend_from_slice(&sequence);
        lengths.push(seq.len());
        seq.push(b'$');
    }
    return (seq, lengths);
}

fn get_sampled_arrays(paths: &[Path<usize, ()>], dict: &Dictionary) -> (Vec<GraphPos>, Vec<usize>) {
    let mut tag_sample = Vec::new();
    let mut sa_sample = Vec::new();

    let n = dict.content.len();

    let sa  = suffix_array(&dict.content);
    let lcp = lcp_array(&dict.content, &sa).decompress();
    let isa = inverse_suffix_array(&sa);

    // let tag = tag_array(&dict, &isa);
    // let rs  = remaining_suffix(&isa, &dict);
    //
    // for mut i in 0..n {
    //     if rs[i] <= overlap { continue; }
    //
    //     if lcp[i] < rs[i] && lcp[i+1] < rs[i] { 
    //         tag_sample.push(tags[i]); 
    //         sa_sample.push(sa_begin);
    //         tag_sample.push(tags[i]);
    //         sa_sample.push(sa_end);
    //     }
    //
    //     while same_prefix {
    //         cumulate();
    //     }
    //     resolve_same();
    //
    // }

    for i in 0..sa.len() {
        println!("{i}\t{}\t{}\t{}", sa[i], lcp[i],
            str::from_utf8(&dict.content[sa[i]..]).unwrap());
    }

    return (tag_sample, sa_sample);
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

fn parse_node_id(node_id: Vec<u8>) -> (String, char) {
    let n = node_id.len();
    let id = str::from_utf8(&node_id[0..n-1]).unwrap().to_owned();
    let sign = node_id[n-1] as char;
    return (id, sign);
}

fn create_map(segments: &[Segment<usize, ()>]) -> HashMap<usize, Vec<u8>> {
    let mut map = HashMap::new();
    for s in segments { map.insert(s.name, s.sequence.clone()); }
    return map;
}

struct ContentSlice { start: usize, end: usize }

struct Dictionary {
    content: Vec<u8>,                   // Segment sequences joined on "$"
    map: HashMap<usize, ContentSlice>   // Segment names to sequences
    // ^- Is this possible to do with &[u8] slices? It does not look possible now. [2023-05-31]
}

impl From<Vec<Segment<usize, ()>>> for Dictionary {
    fn from(segments: Vec<Segment<usize, ()>>) -> Self {
        let mut content = Vec::new();
        let mut map = HashMap::new();

        for segment in segments {
            let slice = ContentSlice {
                start: content.len(),
                end:   content.len() + segment.sequence.len()
            };
            map.insert(segment.name, slice);
            content.extend_from_slice(&segment.sequence);
            // content.push(b'$');
        }
        content.push(b'$');
        return Dictionary{ content, map };
    }
}


fn main() {
    let args = Args::parse();

    let parser: GFAParser<usize, ()> = GFAParser::new();
    let gfa = parser
        .parse_file(&args.gfa_filename)
        .expect("Error parsing GFA file.");

    let dictionary = Dictionary::from(gfa.segments);
    // for (k, v) in &dictionary.map {
    //     println!("{}\t{:?}", k, str::from_utf8(&dictionary.content[v.start..v.end]));
    // }
    println!();
    let (stag, ssa) = get_sampled_arrays(&gfa.paths, &dictionary);
}

