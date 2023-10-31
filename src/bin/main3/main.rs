use clap::Parser;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use std::collections::HashMap;
use std::str;
use std::usize;
use pfg;

mod gp;
mod pred;
mod grammar;
mod mem;

use gp::GraphPos as GraphPos;
use grammar::Grammar;
use mem::MEMReader;

/// Find MEMs in a graph
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// GFA file
    #[arg(short = 'f')]
    gfa_filename: String,
    /// MEMs file
    #[arg(short)]
    mems_filename: String,
    /// pointers
    #[arg(short)]
    ptr_filename: String,
    /// grammar file
    #[arg(short)]
    grammar_filename: String,
    /// Trigger file used for prefix-free suffix array construction
    #[arg(short)]
    trigger_filename: String,
}

fn main() {
    let args = Args::parse();

    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(&args.gfa_filename)
        .expect("Error parsing GFA file.");

    let (start, graph_pos) = parse_graph(&graph);
    let (path_starts, path_ids) = (
        vec![0, 29849, 59696, 89512, 119326],
        vec!["0", "1", "2", "3", "4"]
    );

    let pfdata = pfg::pf::PFData::from_graph(&args.gfa_filename, &args.trigger_filename);

    let mut sampled_tag = Vec::new();
    let mut sampled_sa = Vec::new();
    for (sa, _, _) in pfdata.iter() {
        let i = start.argpred(sa);
        let graph_position = GraphPos{pos: sa - start[i], ..graph_pos[i]};
        sampled_tag.push(graph_position);
        sampled_sa.push(sa);
    }

    let grammar = Grammar::from_file(&args.grammar_filename);

    let mem_reader = MEMReader::new(&args.mems_filename, &args.ptr_filename);
    for (read_id, mems) in mem_reader {
        for mem in mems {
            let graph_positions = get_graph_positions(
                &grammar, &mem, &sampled_tag, &sampled_sa
            );

            let i = path_starts.argpred(mem.2);
            let ref_id = path_ids[i];
            let pos_in_ref = mem.2 - path_starts[i];
            for gp in graph_positions {
                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    read_id, mem.1, ref_id, pos_in_ref, mem.0,
                    gp.id, gp.sign, gp.pos // id, sign, pos
                )
            }
        }
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
        let id = gp.id;
        s += *len.get(&id).unwrap();
    }

    return (start, result);
}

trait Predecessor {
    fn argpred(&self, item: usize) -> usize;
}

impl Predecessor for Vec<usize> {
    fn argpred(&self, item: usize) -> usize {
        let mut l = 0;
        let mut r = self.len();

        while l < r-1 {
            let m = (l + r) / 2;
            if item > self[m] { l = m; }
            else { r = m; }
        }
        return l;
    }
}

fn get_graph_positions(
    grammar: &Grammar, mem: &(usize, usize, usize), tag: &[GraphPos], sa: &[usize]
) -> Vec<GraphPos> {
    let lower = get_lower(grammar, mem, sa);          // included
    let upper = get_upper(grammar, mem, sa);          // excluded

    return list_unique(&tag[lower..upper]);
}

fn get_lower(grammar: &Grammar, mem: &(usize, usize, usize), sa: &[usize]) -> usize {
    let mut l = 0;
    let mut r = sa.len();

    while l < r-1 {
        let m = (l + r) / 2;
        let (e, sa_smaller) = lce(grammar, sa[m], mem.2);
        if e < mem.0 && sa_smaller { l = m; }
        else { r = m; }
    }
    return r;
}

fn get_upper(grammar: &Grammar, mem: &(usize, usize, usize), sa: &[usize]) -> usize {
    let mut l = 0;
    let mut r = sa.len();

    while l < r-1 {
        let m = (l + r) / 2;
        let (e, sa_smaller) = lce(grammar, sa[m], mem.2);
        if e < mem.0 && !sa_smaller { r = m; }
        else { l = m; }
    }

    return r;
}

/// returns (l, f) such that: 
/// seq[s1..s1+l] == seq[s2..s2+l]
/// if seq[s1+l] < seq[s2+l]: f = True
fn lce(grammar: &Grammar, s1: usize, s2: usize) -> (usize, bool) {
    let n = grammar.len();
    if s1 == s2 { return (n-s1, false); }    // the same is not smaller

    let mut l = 0;
    while s1+l < n && s2+l < n && grammar[s1+l] == grammar[s2+l] { l += 1; }
    if s1+l == n { return (l, true); }
    if s2+l == n { return (l, false); }
    if grammar[s1+l] < grammar[s2+l] { return (l, true); }
    if grammar[s1+l] > grammar[s2+l] { return (l, false); }
    panic!("Incorrect grammar is used.");
}

use std::collections::HashSet;
fn list_unique(tag: &[GraphPos]) -> Vec<GraphPos> {
    let mut set: HashSet<GraphPos> = HashSet::new();
    for i in 0..tag.len() { set.insert(tag[i].clone()); }
    let mut res = Vec::new();
    for item in set { res.push(item); }
    return res;
}

#[cfg(test)]
mod tests;

