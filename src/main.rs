use clap::Parser;
use gfa::parser::GFAParser;
use core::panic;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Write;
use std::io::stdout;
use std::path::Path;
use std::path::PathBuf;
use std::str;
use std::usize;
use std::iter::zip;

mod cli;
mod gp;
mod pred;
mod grammar;
mod mem;

use gp::GraphPos as GraphPos;
use grammar::Grammar;
use mem::MEMReader;
use pred::Predecessor;
use cli::Args;

fn main() {
    let args = Args::parse();

    match &args.command {
        cli::Commands::Index { gfa, triggers } => {
            let gfa = PathBuf::from(gfa);
            if !gfa.exists() { panic!("File {} does not exist.", gfa.display()); }
            let triggers = PathBuf::from(triggers);
            if !triggers.exists() { panic!("File {} does not exist.", triggers.display()); }

            let tag = gfa.with_extension("tag");
            // println!("f: {gfa:?} {triggers:?} -> {tag:?}");
            create_tag(&gfa, &triggers, &tag)
        },
        cli::Commands::Align { gfa, reads, output } => {
            let gfa = PathBuf::from(gfa);
            if !gfa.exists() { panic!("File {} does not exist.", gfa.display()); }
            let tag = gfa.with_extension("tag");
            if !tag.exists() { panic!(
                "File {} does not exist. Create it with:\n\n\tmaria index {} -t <triggers.txt>.\n",
                tag.display(), gfa.display()
            )}
            let slp = gfa.with_extension("slp");
            if !slp.exists() { panic!("File {} does not exist.", slp.display()) }

            let reads = PathBuf::from(reads);
            let mems = reads.with_extension("mems");
            let ptrs = reads.with_extension("pointers");

            if !mems.exists() { panic!("File {} does not exist.", mems.display()) }
            if !ptrs.exists() { panic!("File {} does not exist.", ptrs.display()) }

            // println!("f: {gfa:?} {tag:?} {slp:?} {mems:?} {ptrs:?} -> {output:?}");
            if let Some(filename) = output {
                let out = BufWriter::new(
                    File::create(filename).expect("Cannot create output file.")
                );
                align(&gfa, &tag, &slp, &mems, &ptrs, out);
            } else {
                let out = stdout().lock();
                align(&gfa, &tag, &slp, &mems, &ptrs, out);
            }
        }
    }
}

/// f: gfa triggers -> tag
fn create_tag(gfa: &Path, triggers: &Path, tag: &Path) {
    println!("Creating tag array {}", tag.display());
    let (_, _, node_starts, node_names) = process_graph(gfa);
    let (ssa, stag) = get_sampled_arrays(&gfa, &triggers, &node_starts, &node_names);

    let mut writer: BufWriter<File> = BufWriter::new(File::create(tag)
        .unwrap_or_else(|_| panic!("Cannot open file {}", tag.display())));
    for i in 0..ssa.len() {
        writeln!(writer, "{}\t{}{}:{}", ssa[i],
            stag[i].id, stag[i].sign, stag[i].pos // implement Display for tag?
        ).expect("Error while writing tags.");
    }
    writer.flush().expect("Error writing.");
    println!("Tag array successfully created.");
}

/// f: gfa tag slp mems ptrs -> output
fn align<T>(
    gfa: &Path, tag: &Path, grammar: &Path, mems: &Path, ptrs: &Path, mut output: T
) where
    T: Write
{
    let (_, _, node_starts, node_names) = process_graph(gfa);
    let (ssa, stag) = read_tag_array(tag);
    let grammar = Grammar::from_file(grammar);
    let mem_reader = MEMReader::new(mems, ptrs);

    for (read_id, mems) in mem_reader {
        for mem in mems {
            let (sa_values, positions) = get_graph_positions(&grammar, &mem, &stag, &ssa);
            for (sa, _) in zip(sa_values, positions) {
                let (path, plen, pstart, pend) = extract_path(sa, mem.0, &node_starts, &node_names);
                #[allow(clippy::write_literal)] {
                    writeln!(output, "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
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
                    ).expect("Error writing output");
                }
            }
        }
    }
}

// fn get_output(output: Option<String>) -> impl Write {
//     if let Some(x) = output {
//         return File::create(x).unwrap();
//     }
//     // let output = match output {
//     //     None => stdout(),
//     //     Some(name) => File::create(name),
//     // };
//     return stdout();
// }

fn read_tag_array(tag: &Path) -> (Vec<usize>, Vec<GraphPos>) {
    use std::fs::read_to_string;
    let res: (Vec<_>, Vec<_>) = read_to_string(tag).unwrap().lines().map(|line| {
        let z: Vec<_> = line.split('\t').collect();
        if z.len() != 2 { panic!("Tag array line {line} is incorrect.") }
        let sa_value: usize = z[0].parse().unwrap();
        let tag_value: GraphPos = z[1].parse().unwrap();
        (sa_value, tag_value)
    }).unzip();
    return res;
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

fn process_graph<P: AsRef<Path>>(filename: P) -> (
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
fn get_sampled_arrays<P: AsRef<Path>>(
    gfa: P, triggers: P, node_starts: &Vec<usize>, node_names: &[GraphPos]
) -> (Vec<usize>, Vec<GraphPos>) {
    let pfdata = pfg::pf::PFData::from_graph(gfa.as_ref().to_str().unwrap(), triggers.as_ref().to_str().unwrap());
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

