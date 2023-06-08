use bio::data_structures::suffix_array::lcp as lcp_array;
use clap::Parser;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use maria::arrays::SuffixArray;
use std::collections::HashMap;
use std::str;
use std::usize;

mod gp;
mod pred;
mod block;
mod pf;

use gp::GraphPos as GraphPos;
use pred::Predecessor;
use block::Block;

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

fn main() {
    // let args = Args::parse();

    let parser: GFAParser<usize, ()> = GFAParser::new();
    let arbitrary_graph = parser.parse_file("data/pftag/test_arbitrary.gfa")
        .expect("Error parsing GFA file.");

    let (start, graph_pos) = parse_graph(&arbitrary_graph);

    let gfa = parser.parse_file("data/pftag/test.gfa")
        .expect("Error parsing GFA file.");

    let overlap = 1;
    let segments: Vec<Vec<u8>> = pf::parse_segments(&gfa.segments);
    let paths: Vec<Vec<usize>> = pf::parse_paths(&gfa.paths);

    let segment_join = pf::join(&segments);

    let sa  = SuffixArray::create(&*segment_join);
    let lcp = lcp_array(&segment_join, &sa).decompress();
    let isa = pf::permutation_invert(&sa);
    let id  = pf::get_node_ids(&segment_join, &isa);
    let pos = pf::get_node_pos(&segment_join, &isa);

    let path_join = pf::join(&paths);

    let len = pf::get_lengths(&segments);
    let mut rc_rank = pf::get_right_context_rank(&path_join, segments.len());
    let mut seq_pos = pf::get_sequence_position(&path_join, &len, overlap);

    for i in 0..segments.len() {
        let iperm = pf::argsort(&rc_rank[i]);
        rc_rank[i] = pf::permutation_apply(&iperm, &rc_rank[i]);
        seq_pos[i] = pf::permutation_apply(&iperm, &seq_pos[i]);
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

#[cfg(test)]
mod tests;

