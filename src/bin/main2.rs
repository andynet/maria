use gfa::parser::GFAParser;
use bio::data_structures::suffix_array::lcp as lcp_array;
use maria::inverse_suffix_array;
use clap::Parser;
use std::str;
use maria::arrays::SuffixArray;
use gfa::gfa::Segment;
use gfa::gfa::Path;
use std::ops::Add;

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
    let gfa = parser.parse_file("data/pftag/test.gfa")
        .expect("Error parsing GFA file.");

    let overlap = 1;
    let segments: Vec<Vec<u8>> = parse_segments(&gfa.segments);
    let paths: Vec<Vec<usize>> = parse_paths(&gfa.paths);

    let segment_join = join(&segments);

    let sa  = SuffixArray::create(&*segment_join);
    let lcp = lcp_array(&segment_join, &sa).decompress();
    let isa = inverse_suffix_array(&sa);
    let id  = get_node_ids(&segment_join, &isa);
    let pos = get_node_pos(&segment_join, &isa);

    let path_join = join(&paths);

    let len     = get_lengths(&segments);
    let freq    = get_frequency(&path_join, segments.len());
    let seq_pos = get_sequence_position(&path_join, &len, overlap);
    let rc_rank = get_right_context_rank(&path_join, segments.len());

    print_tag_array(&lcp, &id, &pos, overlap, &len, &freq, &seq_pos, &rc_rank);
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

fn get_frequency(path_join: &[usize], size: usize) -> Vec<usize> {
    let mut result = vec![0; size];
    for id in path_join {
        if *id >= 2 { result[id-2] += 1; }
    }
    return result;
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
    let isa = inverse_suffix_array(&sa);

    for (i, id) in path_join.iter().enumerate() {
        if *id >= 2 {
            result[id-2].push(isa[i+1]);
        }
    }
    return result;
}

fn print_tag_array(
        lcp: &[isize],
         id: &[usize],
        pos: &[usize],
    overlap: usize,
        len: &[usize],
       freq: &[usize],
    seq_pos: &[Vec<usize>],
    rc_rank: &[Vec<usize>],
) {
    let mut i = len.len() + 1;
    while i < id.len() {
        let remaining = len[id[i]] - pos[i];
        if remaining <= overlap { i += 1; continue; }
        let mut j = i+1;
        while lcp[j] >= remaining as isize { j += 1; }

        let mut identical_suffixes = Vec::new();
        for k in i..j {
            for l in 0..freq[id[k]] {
                identical_suffixes.push((
                    rc_rank[id[k]][l],
                    id[k],
                    pos[k],
                    seq_pos[id[k]][l] + pos[k]
                ))
            }
        }
        identical_suffixes.sort_unstable();
        for (_, id, pos, sa) in &identical_suffixes {
            println!("{}\t{}\t{}", sa, id, pos);
        }
        i = j;
    }
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
        let isa = inverse_suffix_array(&sa);

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
        let freq    = get_frequency(&path_join, segments.len());
        let seq_pos = get_sequence_position(&path_join, &len, overlap);
        let rc_rank = get_right_context_rank(&path_join, segments.len());

        for id in 0..segments.len() {
            println!("{}\t{}\t{}\t{:?}\t{:?}", 
                id, len[id], freq[id], seq_pos[id], rc_rank[id]
            );
        }
        println!();

        print_tag_array(&lcp, &id, &pos, overlap, &len, &freq, &seq_pos, &rc_rank);
    }
}

