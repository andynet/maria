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
    let args = Args::parse();

    let parser: GFAParser<usize, ()> = GFAParser::new();
    let gfa = parser
        .parse_file(&args.gfa_filename)
        .expect("Error parsing GFA file.");

    let v: Vec<u8> = vec![45,5,7,6,68,8,0];
    let sa = SuffixArray::create(&v[..]);
    // let dictionary = Dictionary::from(gfa.segments);
    // for (k, v) in &dictionary.map {
    //     println!("{}\t{:?}", k, str::from_utf8(&dictionary.content[v.start..v.end]));
    // }
    // let (stag, ssa) = get_sampled_arrays(&gfa.paths, &dictionary);
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

// these two joins can be merged into one
fn join_segments(segments: &[Vec<u8>]) -> Vec<u8> {
    let mut segment_join = Vec::new();
    for segment in segments {
        segment_join.extend_from_slice(segment);
    }
    segment_join.push(b'$');

    return segment_join; 
}

fn join_paths(paths: &[Vec<usize>]) -> Vec<usize> { 
    let mut path_join = Vec::new();
    for p in paths {
        let p: Vec<usize> = p.iter().map(|x| x + 1).collect();
        path_join.extend_from_slice(&p);
        path_join.push(0);
    }
    return path_join; 
}

fn join_u8(slice: &[Vec<u8>]) -> Vec<u8> {
    let mut result = Vec::new();
    for vector in slice {
        let v: Vec<u8> = vector.iter().map(|x| x + 2).collect();
        result.extend_from_slice(&v);
        result.push(1);
    }
    result.push(0);
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

fn join_usize(slice: &[Vec<usize>]) -> Vec<usize> {
    let mut result = Vec::new();
    for vector in slice {
        let v: Vec<usize> = vector.iter().map(|x| x + 2).collect();
        result.push(1);
    }
    result.push(0);
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

fn get_pfgraph_pos(segment_join: &[u8], isa: &[usize], len: &[usize]) 
    -> (Vec<usize>, Vec<usize>) 
{ 
    let n = segment_join.len();
    let mut nodeid = vec![0; n];
    let mut nodepos = vec![0; n];
    let mut id = 0;
    let mut pos = 0;
    for i in 0..n { 
        nodeid[isa[i]] = id;
        nodepos[isa[i]] = pos;

        pos += 1;
        if pos == len[id] {
            id += 1;
            pos = 0;
        }
    }

    return (nodeid, nodepos); 
}

fn get_frequency(path_join: &[usize], m: usize) -> Vec<usize> {
    let mut freq = vec![0; m];
    for id in path_join {
        let true_id = (id + (m-1)) % m;
        freq[true_id] += 1;
    }
    return freq;
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

        let segment_join = join_segments(&segments);
        let path_join = join_paths(&paths);

        let n = segment_join.len();
        let m = &gfa.segments.len() + 1; // one fake sentinel segment

        let sa  = SuffixArray::create(&*segment_join);
        let lcp = lcp_array(&segment_join, &sa).decompress();
        let isa = inverse_suffix_array(&sa);
        let len = get_lengths(&segments);

        let (id, npos) = get_pfgraph_pos(&segment_join, &isa, &len);

        for i in 0..n {
            println!("{}\t{}\t{}\t{}\t{}\t{}", i, sa[i], lcp[i], id[i], npos[i],
                str::from_utf8(&segment_join[sa[i]..]).unwrap())
        }
        println!();

        let freq = get_frequency(&path_join, m);
        println!("{:?}", path_join);
        println!("{:?}", len);
        // let seq_pos = get_sequence_position(&path_join, &len, m);

        let mut orig_pos = vec![Vec::new(); m];
        let mut pos = 0;
        for p in paths.iter() {
            for node in p {
                orig_pos[*node].push(pos);
                pos += len[*node] - 1;
            }
            orig_pos[m-1].push(pos);
            pos += 1;
        }

        let mut str_pos = vec![Vec::new(); m];
        let mut pos = 0;
        for p in paths.iter() {
            for node in p {
                str_pos[*node].push(pos);
                pos += 1;
            }
            let last = str_pos.last_mut().unwrap();
            last.push(pos);
            pos += 1;
        }

        let path_sa = SuffixArray::create(&path_join[..]);
        let path_isa = inverse_suffix_array(&path_sa);
        let mut rcr = Vec::new();   // right context rank
        for s_pos in &str_pos {
            let i = rcr.len();
            rcr.push(Vec::new());
            for p in s_pos {
                let rc_start = (p+1)%path_isa.len();
                rcr[i].push(path_isa[rc_start]);
            }
        }

        for i in 0..m {
            println!("{}\t{}\t{}\t{:?}\t{:?}\t{:?}",
                i, len[i], freq[i], str_pos[i], orig_pos[i], rcr[i]
            );
        }
        println!();

        let mut i = 0;
        while i < n {
            let remaining = len[id[i]] - npos[i];
            if remaining <= overlap { i+=1; continue; }
            let mut j = i+1;
            while lcp[j] >= remaining as isize {
                let rem = len[id[j]] - npos[j];
                if rem <= overlap { break; }
                j += 1;
            }

            // process table from i to j
            let mut identical_suffixes = Vec::new();
            for k in i..j {
                let nid = id[k];
                for l in 0..freq[nid] {
                    identical_suffixes.push((
                        rcr[nid][l],
                        nid,
                        npos[k],
                        orig_pos[nid][l] + npos[k]
                    ))
                }
            }
            identical_suffixes.sort_unstable();
            for (_, nid, np, sa) in &identical_suffixes {
                println!("{}\t{}\t{}", sa, nid, np);
            }
            identical_suffixes.clear();
            // end processing
            i = j;
        }
    }
}

