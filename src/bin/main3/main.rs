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
mod pf;
mod grammar;

use gp::GraphPos as GraphPos;
use pred::Predecessor;
use grammar::Grammar;

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
    let args = Args::parse();

    let (start, graph_pos) = parse_graph("data/pftag/test_arbitrary.gfa");

    // let (sampled_sa, sampled_tag) = tag_array(
    //     "data/pftag/test_arbitrary.gfa",
    //     "../../../data/pftag/triggers.txt"
    // );

    let grammar = Grammar::from_file("data/pftag/test_join.txt.plainslp");

    // for MEM in MEMs:
    //      
}

fn parse_graph(graph_file: &str) -> (Vec<usize>, Vec<GraphPos>) {
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(graph_file)
        .expect("Error parsing GFA file.");

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

