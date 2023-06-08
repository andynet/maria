use bio::data_structures::suffix_array::lcp as lcp_array;
use clap::Parser;
use gfa::gfa::GFA;
use gfa::gfa::Path;
use gfa::gfa::Segment;
use gfa::parser::GFAParser;
use maria::arrays::SuffixArray;
use std::collections::HashMap;
use std::fmt::Display;
use std::num::ParseIntError;
use std::ops::Add;
use std::str;
use std::str::FromStr;
use std::usize;

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

#[derive(Clone)]
enum Direction {
    Forward,
    RevComp
}

use std::fmt;
impl Display for Direction {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Direction::Forward => { write!(f, "+") },
            Direction::RevComp => { write!(f, "-") }
        }
    } 
}

#[derive(Clone)]
struct GraphPos {
    node_id: usize,
    direction: Direction,
    node_pos: usize,
}

#[derive(Debug)]
struct ParseGraphPosError;
impl From<ParseIntError> for ParseGraphPosError {
    fn from(_: ParseIntError) -> Self { ParseGraphPosError }
}

impl FromStr for GraphPos {
    type Err = ParseGraphPosError;

    fn from_str(s: &str) -> Result<Self, Self::Err> {
        let node_id: usize = s[..s.len()-1].parse()?;
        let direction;
        match s.bytes().last().ok_or(ParseGraphPosError)? {
            b'+' => { direction = Direction::Forward; },
            b'-' => { direction = Direction::RevComp; },
            _    => { return Err(ParseGraphPosError); }
        }
        let node_pos = 0;
        Ok(GraphPos { node_id, direction, node_pos })
    }
}


fn parse_graph(graph: &GFA<usize, ()>) -> (Vec<usize>, Vec<GraphPos>) {
    let mut len = HashMap::new();
    for seg in &graph.segments { len.insert(seg.name, seg.sequence.len()); }

    let mut result = Vec::new();
    for path in &graph.paths {
        let p = str::from_utf8(&path.segment_names).unwrap();
        let p: Vec<GraphPos> = p.split(',').map(|x| x.parse().unwrap()).collect();
        result.extend_from_slice(&p);
    }

    let mut start = Vec::new();
    let mut s = 0;
    for gp in &result {
        start.push(s);
        let id = gp.node_id;
        s += *len.get(&id).unwrap();
    }

    return (start, result);
}

trait Predecessor {
    fn argpred(&self, item: usize) -> usize;
}

impl Predecessor for Vec<usize> {
    fn argpred(&self, item: usize) -> usize {
        let mut i = self.len() - 1;
        while self[i] > item { i -= 1; }
        return i;
    }
}

fn main() {
    // let args = Args::parse();

    let parser: GFAParser<usize, ()> = GFAParser::new();
    let arbitrary_graph = parser.parse_file("data/pftag/test_arbitrary.gfa")
        .expect("Error parsing GFA file.");

    let (start, graph_pos) = parse_graph(&arbitrary_graph);

    let gfa = parser.parse_file("data/pftag/test.gfa")
        .expect("Error parsing GFA file.");

    let overlap = 1;
    let segments: Vec<Vec<u8>> = parse_segments(&gfa.segments);
    let paths: Vec<Vec<usize>> = parse_paths(&gfa.paths);

    let segment_join = join(&segments);

    let sa  = SuffixArray::create(&*segment_join);
    let lcp = lcp_array(&segment_join, &sa).decompress();
    let isa = permutation_invert(&sa);
    let id  = get_node_ids(&segment_join, &isa);
    let pos = get_node_pos(&segment_join, &isa);

    let path_join = join(&paths);

    let len = get_lengths(&segments);
    let mut rc_rank = get_right_context_rank(&path_join, segments.len());
    let mut seq_pos = get_sequence_position(&path_join, &len, overlap);

    for i in 0..segments.len() {
        let iperm = argsort(&rc_rank[i]);
        rc_rank[i] = permutation_apply(&iperm, &rc_rank[i]);
        seq_pos[i] = permutation_apply(&iperm, &seq_pos[i]);
    }

    let mut i = len.len() + 1;
    while i < id.len() {
        let remaining = len[id[i]] - pos[i];
        if remaining <= overlap { i += 1; continue; }
        let mut j = i+1;
        while lcp[j] >= remaining as isize { j += 1; }

        let block = Block::new(&id[i..j], &pos[i..j], &seq_pos, &rc_rank);
        for (sa, id, pos) in block {
            let idx = start.argpred(sa);
            let node_start = start[idx];
            let node = &graph_pos[idx];
            println!("{}\t{}\t{}\t{}\t{}\t{}", 
                sa, id, pos, node.node_id, node.direction, sa - node_start);
        }
        i = j;
    }

    // let grammar = parse_grammar()
    // for MEM in MEMs:
    //      
}

fn parse_segments(segments: &[Segment<usize, ()>]) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    for segment in segments {
        result.push(segment.sequence.clone());
    }
    return result;
}

fn parse_paths(paths: &[Path<usize, ()>]) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    for path in paths {
        let p = str::from_utf8(&path.segment_names).unwrap();
        let to_usize = |x: &str| (x[0..x.len()-1].parse::<usize>()
            .expect("Cannot parse path."));

        let p: Vec<_> = p.split(',').map(to_usize).collect();
        result.push(p);
    }
    return result;
}

fn join<T>(slice: &[Vec<T>]) -> Vec<T> 
where
    T: From<u8> + Add<T, Output = T> + Copy
{
    let mut result = Vec::new();
    for vector in slice {
        let v: Vec<T> = vector.iter().map(|&x| x + 2.into()).collect();
        result.extend_from_slice(&v);
        result.push(1.into());
    }
    result.push(0.into());
    return result;
}

fn get_lengths(segments: &[Vec<u8>]) -> Vec<usize> { 
    let mut len = Vec::new();
    for segment in segments {
        len.push(segment.len());
    }
    len.push(1);
    return len;
}

fn get_node_ids(segment_join: &[u8], isa: &[usize]) -> Vec<usize> {
    let mut result = vec![0; segment_join.len()];
    let mut id = 0;
    for i in 0..segment_join.len() {
        result[isa[i]] = id;
        if segment_join[i] == 1 { id += 1; }
    }
    return result;
}

fn get_node_pos(segment_join: &[u8], isa: &[usize]) -> Vec<usize> {
    let mut result = vec![0; segment_join.len()];
    let mut pos = 0;
    for i in 0..segment_join.len() {
        result[isa[i]] = pos;
        pos += 1;
        if segment_join[i] == 1 { pos = 0; }
    }
    return result;
}

fn get_sequence_position(path_join: &[usize], len: &[usize], overlap: usize) -> Vec<Vec<usize>> {
    let mut result = vec![Vec::new(); len.len()];
    let mut start = 0;
    for id in path_join {
        if *id >= 2 {
            result[id-2].push(start);
            start += len[id-2] - overlap;
        }
    }
    return result;
}

fn get_right_context_rank(path_join: &[usize], size: usize) -> Vec<Vec<usize>> {
    let mut result = vec![Vec::new(); size];

    let sa = SuffixArray::create(&path_join[..]);
    let isa = permutation_invert(&sa);

    for (i, id) in path_join.iter().enumerate() {
        if *id >= 2 {
            result[id-2].push(isa[i+1]);
        }
    }
    return result;
}

struct Block<'a> {
         id: &'a[usize],
        pos: &'a[usize],
    seq_pos: &'a[Vec<usize>],
    rc_rank: &'a[Vec<usize>],
        idx: Vec<usize>,
        rcr: Vec<usize>
}

impl<'a> Block<'a> {
    fn new(
        id: &'a[usize], pos: &'a[usize], 
        seq_pos: &'a[Vec<usize>], rc_rank: &'a[Vec<usize>]
    ) -> Block<'a> {

        let mut rcr = Vec::new();
        for i in id {
            rcr.push(rc_rank[*i][0]);
        }
        Block { 
            id, pos, seq_pos, rc_rank, 
            idx: vec![0; id.len()],
            rcr
        }
    }
}

impl<'a> Iterator for Block<'a> {
    type Item = (usize, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let k = argmin(&self.rcr);
        if self.rcr[k] == usize::MAX { return None; }
        let result = (
            self.seq_pos[self.id[k]][self.idx[k]] + self.pos[k],
            self.id[k],
            self.pos[k]
        );

        self.idx[k] += 1;
        let val = self.rc_rank[self.id[k]].get(self.idx[k]);
        match val {
            None => { self.rcr[k] = usize::MAX; },
            Some(v) => { self.rcr[k] = *v; }
        }

        return Some(result);
    }
}

fn permutation_invert(perm: &[usize]) -> Vec<usize> {
    let mut inverse = vec![0; perm.len()];
    for i in 0..perm.len() { inverse[perm[i]] = i; }
    return inverse;
}

fn permutation_apply(perm: &[usize], data: &[usize]) -> Vec<usize> {
    let mut result = vec![0; perm.len()];
    for i in 0..perm.len() {
        result[perm[i]] = data[i];
    }
    return result;
}

fn argsort(data: &[usize]) -> Vec<usize> {
    let mut v: Vec<(usize, usize)> = data.iter().enumerate().map(|(i, x)| (*x, i)).collect();
    v.sort_unstable();
    let v: Vec<usize> = v.iter().map(|(_, i)| *i).collect();
    let v = permutation_invert(&v);
    return v;
}

fn argmin(data: &[usize]) -> usize {
    let mut min_pos = 0;
    let mut min_val = usize::MAX;
    for (p, &v) in data.iter().enumerate() {
        if v < min_val { min_pos = p; min_val = v; }
    }
    return min_pos;
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let parser: GFAParser<usize, ()> = GFAParser::new();
        let gfa = parser.parse_file("data/pftag/test.gfa")
            .expect("Error parsing GFA file.");

        let overlap = 1;
        let segments: Vec<Vec<u8>> = parse_segments(&gfa.segments);
        let paths: Vec<Vec<usize>> = parse_paths(&gfa.paths);

        let segment_join = join(&segments);
        let path_join = join(&paths);

        let sa  = SuffixArray::create(&*segment_join);
        let lcp = lcp_array(&segment_join, &sa).decompress();
        let isa = permutation_invert(&sa);

        let id  = get_node_ids(&segment_join, &isa);
        let pos = get_node_pos(&segment_join, &isa);

        let mut segment_join_repr = segment_join.clone();
        for c in segment_join_repr.iter_mut() {
            match *c {
                0   => { *c = b'$'; },
                1   => { *c = b'#'; },
                2.. => { *c -= 2; }
            }
        }
        for i in 0..sa.len() {
            println!("{}\t{}\t{}\t{}\t{}\t{}", i, sa[i], lcp[i], id[i], pos[i],
                str::from_utf8(&segment_join_repr[sa[i]..]).unwrap());
        }
        println!();

        let len     = get_lengths(&segments);
        let mut rc_rank = get_right_context_rank(&path_join, segments.len());
        let mut seq_pos = get_sequence_position(&path_join, &len, overlap);

        for i in 0..segments.len() {
            let iperm = argsort(&rc_rank[i]);
            rc_rank[i] = permutation_apply(&iperm, &rc_rank[i]);
            seq_pos[i] = permutation_apply(&iperm, &seq_pos[i]);
        }

        for id in 0..segments.len() {
            println!("{}\t{}\t{:?}\t{:?}", 
                id, len[id], seq_pos[id], rc_rank[id]
            );
        }
        println!();

        // print_tag_array(&lcp, &id, &pos, overlap, &len, &seq_pos, &rc_rank);
    }

    #[test]
    fn test_arbitrary() {
        let parser: GFAParser<usize, ()> = GFAParser::new();
        let arbitrary_graph = parser.parse_file("data/pftag/test_arbitrary.gfa")
            .expect("Error parsing GFA file.");

        let (seq_pos, graph_pos) = parse_graph(&arbitrary_graph);

        for i in 0..seq_pos.len() {
            println!("{}\t{}", seq_pos[i], graph_pos[i].node_id);
        }
    }
}

