use clap::Parser;
use gfa::parser::GFAParser;
use std::collections::HashMap;
use std::str;
use std::usize;
use std::iter::zip;

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
    let (_, _, node_starts, node_names) = process_graph2(&args.gfa_filename);
    let (ssa, stag) = get_sampled_arrays(&args.gfa_filename, &args.trigger_filename, &node_starts, &node_names);
    let grammar = Grammar::from_file(&args.grammar_filename);
    let mem_reader = MEMReader::new(&args.mems_filename, &args.ptr_filename);

    for (read_id, mems) in mem_reader {
        for mem in mems {
            let (sa_values, positions) = get_graph_positions(&grammar, &mem, &stag, &ssa);
            for (sa, _) in zip(sa_values, positions) {
                let (path, plen, pstart, pend) = extract_path(sa, mem.0, &node_starts, &node_names);
                #[allow(clippy::print_literal)] {
                    println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
                        read_id,        // string:  Query sequence name
                        150,            // int:     Query sequence length
                        mem.1,          // int:     Query start (0-based; closed)
                        mem.1 + mem.0,  // int:     Query end (0-based; open)
                        '+',            // char:    Strand relative to the path: "+" or "-"
                        path,           // string:  Path matching
                        plen,           // int:     Path length
                        pstart,         // int:     Start position on the path (0-based; closed)
                        pend,           // int:     End position on the path (0-based; open)
                        mem.0,          // int:     Number of residue matches
                        mem.0,          // int:     Alignment block length
                        60              // int:     Mapping quality (0-255; 255 for missing)
                    );
                }
            }
        }
    }
}

fn extract_path(
    sa_value: usize, seq_len: usize,
    node_starts: &Vec<usize>, node_names: &[GraphPos]
) -> (String, usize, usize, usize) {
    let mut i = node_starts.argpred(sa_value);
    let start = node_starts[i];

    let pstart = sa_value - start;
    let pend = pstart + seq_len;

    let mut path = String::new();
    while node_starts[i] < start + pend {
        path.push_str(&node_names[i].to_path());
        i += 1;
    }
    let plen = node_starts[i] - start;

    return (path, plen, pstart, pend);
}

#[test]
fn extract_path_correctly_handles_the_end() {
    let sa_value = 5;
    let seq_len = 3;
    let node_starts = vec![0, 10];
    let node_names = vec![GraphPos::default()];

    extract_path(sa_value, seq_len, &node_starts, &node_names);
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
    node_starts.push(start); // sentinel
    return (path_starts, path_names, node_starts, node_names);
}

/// Returns sampled suffix array and sampled tag array.
/// Both array are sampled at the starts and ends of run boundaries of the tag array.
fn get_sampled_arrays(
    gfa: &str, triggers: &str, node_starts: &Vec<usize>, node_names: &[GraphPos]
) -> (Vec<usize>, Vec<GraphPos>) {
    let pfdata = pfg::pf::PFData::from_graph(gfa, triggers);
    let mut iterator = pfdata.iter();

    let mut sampled_tag = Vec::new();
    let mut sampled_suf = Vec::new();

    let (sa, _, _) = iterator.next().unwrap();
    let i = node_starts.argpred(sa);
    let gp = GraphPos{pos: sa - node_starts[i], ..node_names[i]};
    sampled_tag.push(gp);
    sampled_suf.push(sa);

    let mut old_gp = gp;
    let mut old_sa = sa;

    for (sa, _, _) in iterator {
        let i = node_starts.argpred(sa);
        let gp = GraphPos{pos: sa - node_starts[i], ..node_names[i]};
        if old_gp != gp {
            sampled_tag.push(old_gp);
            sampled_suf.push(old_sa);
            sampled_tag.push(gp);
            sampled_suf.push(sa);
            old_gp = gp;
        }
        old_sa = sa;
    }

    sampled_tag.push(old_gp);
    sampled_suf.push(old_sa);

    return (sampled_suf, sampled_tag);
}

fn get_graph_positions(
    grammar: &Grammar, mem: &(usize, usize, usize), tag: &[GraphPos], sa: &[usize]
) -> (Vec<usize>, Vec<GraphPos>) {
    let lower = get_lower(grammar, mem, sa);          // included
    let upper = get_upper(grammar, mem, sa);          // excluded

    return list_unique(&sa[lower..upper], &tag[lower..upper]);
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

fn list_unique(sa: &[usize], tag: &[GraphPos]) -> (Vec<usize>, Vec<GraphPos>) {
    let mut map: HashMap<GraphPos, usize> = HashMap::new();
    for i in 0..tag.len() { map.insert(tag[i], sa[i]); }
    let mut suf_uniq = Vec::new();
    let mut tag_uniq = Vec::new();
    for (&k, &v) in map.iter() {
        tag_uniq.push(k);
        suf_uniq.push(v);
    }
    return (suf_uniq, tag_uniq);
}

#[cfg(test)]
mod tests;

