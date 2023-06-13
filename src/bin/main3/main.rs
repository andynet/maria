use clap::Parser;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use std::collections::HashMap;
use std::str;
use std::usize;
use maria::pf;

mod gp;
mod pred;
mod grammar;
mod mem;

use gp::GraphPos as GraphPos;
use grammar::Grammar;

/// Find MEMs in a graph
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// GFA file
    #[arg(short)]
    gfa_filename: String,
    /// Trigger file used for prefix-free suffix array construction
    #[arg(short)]
    trigger_filename: String,
    /// MEMs file
    #[arg(short)]
    mems_filename: String,
    /// pointers
    #[arg(short)]
    ptr_filename: String,
}

fn main() {
    let args = Args::parse();

    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(args.gfa_filename)    // "data/pftag/test_arbitrary.gfa"
        .expect("Error parsing GFA file.");

    let (start, graph_pos) = parse_graph(&graph);

    let (trigs, trigs_size) = pf::load_trigs(&args.trigger_filename);
    let triggers = pf::get_triggers(&trigs, trigs_size);

    let mut segments = HashMap::new();
    let mut paths = Vec::new();

    for path in &graph.paths {
        let mut seq = pf::reconstruct_path(path, &graph);
        let v = vec![b'.'; trigs_size];
        seq.extend_from_slice(&v);
        pf::split_prefix_free(&seq, &triggers, &mut segments, &mut paths);
    }

    let (segments, paths) = pf::normalize(segments, paths);

    let pfdata = pf::PFData::new(&segments, &paths, trigs_size);

    let mut tag = Vec::new();
    for (sa, _, _) in pfdata.iter() {
        let i = start.argpred(sa);
        let gpos = &graph_pos[i];
        tag.push((sa, gpos.node_id, gpos.direction.clone(), sa - i));
    }

    let grammar = Grammar::from_file("data/pftag/test_join.txt.plainslp");

    // let mems = MEMReader::new(&args.mems_filename, &args.ptr_filename);
    // for mem in mems {
    //
    //     let gp: Vec<GraphPos> = get_graph_positions(&seq, &ref_mem, &stag, &ssa);
    //     for gpos in gp {
    //         let (node_id, sign): (String, char) = parse_node_id(gpos.node_id);
    //         println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
    //             read_id, mem.0, adj, ref_pos, mem.1, node_id, sign, gpos.pos
    //         )
    //     }
    // }
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

#[cfg(test)]
mod tests;

