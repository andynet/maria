use clap::Parser;
use gfa::gfa::GFA;
use gfa::parser::GFAParser;
use std::collections::HashMap;
use std::collections::HashSet;
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
use pred::Predecessor;

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

    let (path_starts, path_names, node_starts, node_names) = process_graph2(
        &args.gfa_filename
    );

    let (ssa, stag) = get_sampled_arrays(
        &args.gfa_filename, &args.trigger_filename, &node_starts, &node_names
    );

    let grammar = Grammar::from_file(&args.grammar_filename);

    let mem_reader = MEMReader::new(&args.mems_filename, &args.ptr_filename);

    for (read_id, mems) in mem_reader {
        for mem in mems {
            let graph_positions = get_graph_positions(
                &grammar, &mem, &stag, &ssa
            );

            let i = path_starts.argpred(mem.2);
            let ref_id = &path_names[i];
            let pos_in_ref = mem.2 - path_starts[i];
            for gp in graph_positions {
                println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                    read_id, mem.1, ref_id, pos_in_ref, mem.0,
                    gp.id, gp.sign, gp.pos // id, sign, pos
                );
                // println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                //     read_id,        // string:  Query sequence name
                //     qlen,           // int:     Query sequence length
                //     mem.1,          // int:     Query start (0-based; closed)
                //     mem.1 + mem.0,  // int:     Query end (0-based; open)
                //     '+',            // char:    Strand relative to the path: "+" or "-"
                //     path,           // string:  Path matching
                //     plen,           // int:     Path length
                //     pstart,         // int:     Start position on the path (0-based; closed)
                //     pend,           // int:     End position on the path (0-based; open)
                //     nres,           // int:     Number of residue matches
                //     alen,           // int:     Alignment block length
                //     qual            // int:     Mapping quality (0-255; 255 for missing)
                // );
            }
        }
    }

}

fn process_graph(filename: &str) -> (Vec<usize>, Vec<String>, Vec<usize>, Vec<GraphPos>) {
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(filename)
        .expect("Error parsing GFA file.");

    let (start, graph_pos) = parse_graph(&graph);
    let (path_starts, path_ids) = (
        vec![0, 29850, 59697, 89513, 119327],
        vec![
            "ENA|MW565758|MW565758.1".to_string(),
            "ENA|MW565759|MW565759.1".to_string(),
            "ENA|MW565760|MW565760.1".to_string(),
            "ENA|MW565761|MW565761.1".to_string(),
            "ENA|LR883856|LR883856.1".to_string()
        ]
    );
    return (path_starts, path_ids, start, graph_pos);
}

fn process_graph2(filename: &str) -> (
    Vec<usize>, Vec<String>, Vec<usize>, Vec<GraphPos>
) {
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(filename).expect("Error parsing GFA file.");

    let mut seg_len = HashMap::new();
    for seg in &graph.segments { seg_len.insert(seg.name, seg.sequence.len()); }

    let mut path_starts = Vec::new();
    let mut path_names = Vec::new();
    let mut node_starts = Vec::new();
    let mut node_names = Vec::new();

    let mut start = 0;
    for path in &graph.paths {
        path_starts.push(start);
        path_names.push(str::from_utf8(&path.path_name).unwrap().to_string());

        let segments: Vec<GraphPos> =
            str::from_utf8(&path.segment_names).unwrap()
            .split(',').map(|x| x.parse().unwrap()).collect();

        for node in segments {
            node_starts.push(start);
            start += *seg_len.get(&node.id).unwrap();
            node_names.push(node);
        }
    }
    return (path_starts, path_names, node_starts, node_names);
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

fn get_sampled_arrays(
    gfa: &str, triggers: &str, node_starts: &Vec<usize>, node_names: &[GraphPos]
) -> (Vec<usize>, Vec<GraphPos>) {
    let pfdata = pfg::pf::PFData::from_graph(gfa, triggers);

    let mut sampled_tag = Vec::new();
    let mut sampled_sa = Vec::new();
    for (sa, _, _) in pfdata.iter() {
        // TODO: sample tag array
        let i = node_starts.argpred(sa);
        let graph_position = GraphPos{pos: sa - node_starts[i], ..node_names[i]};
        sampled_tag.push(graph_position);
        sampled_sa.push(sa);
    }

    return (sampled_sa, sampled_tag);
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

fn list_unique(tag: &[GraphPos]) -> Vec<GraphPos> {
    let mut set: HashSet<GraphPos> = HashSet::new();
    for i in 0..tag.len() { set.insert(tag[i].clone()); }
    let mut res = Vec::new();
    for item in set { res.push(item); }
    return res;
}

#[cfg(test)]
mod tests;

