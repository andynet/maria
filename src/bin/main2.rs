use gfa::parser::GFAParser;
use bio::data_structures::suffix_array::lcp as lcp_array;
use maria::inverse_suffix_array;
use clap::Parser;
use std::str;
use maria::arrays::SuffixArray;

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

// fn get_sampled_arrays(paths: &[Path<usize, ()>], dict: &Dictionary) -> (Vec<GraphPos>, Vec<usize>) {
//     let mut tag_sample = Vec::new();
//     let mut sa_sample = Vec::new();
//
//     // let tag = tag_array(&dict, &isa);
//     // let rs  = remaining_suffix(&isa, &dict);
//     //
//     // for mut i in 0..n {
//     //     if rs[i] <= overlap { continue; }
//     //
//     //     if lcp[i] < rs[i] && lcp[i+1] < rs[i] { 
//     //         tag_sample.push(tags[i]); 
//     //         sa_sample.push(sa_begin);
//     //         tag_sample.push(tags[i]);
//     //         sa_sample.push(sa_end);
//     //     }
//     //
//     //     while same_prefix {
//     //         cumulate();
//     //     }
//     //     resolve_same();
//     //
//     // }
//
//     return (tag_sample, sa_sample);
// }
//
// struct ContentSlice { start: usize, end: usize }
//
// struct Dictionary {
//     content: Vec<u8>,                   // Segment sequences joined on "$"
//     map: HashMap<usize, ContentSlice>   // Segment names to sequences
//     // ^- Is this possible to do with &[u8] slices? It does not look possible now. [2023-05-31]
// }
//
// impl From<Vec<Segment<usize, ()>>> for Dictionary {
//     fn from(segments: Vec<Segment<usize, ()>>) -> Self {
//         let mut content = Vec::new();
//         let mut map = HashMap::new();
//
//         for segment in segments {
//             let slice = ContentSlice {
//                 start: content.len(),
//                 end:   content.len() + segment.sequence.len()
//             };
//             map.insert(segment.name, slice);
//             content.extend_from_slice(&segment.sequence);
//             // content.push(b'$');
//         }
//         content.push(b'$');
//         return Dictionary{ content, map };
//     }
// }
//

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

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let parser: GFAParser<usize, ()> = GFAParser::new();
        let gfa = parser.parse_file("data/pftag/test.gfa")
            .expect("Error parsing GFA file.");

        // get content
        let mut dict_content = Vec::new();
        let mut dict_lengths = Vec::new();
        for segment in &gfa.segments {
            dict_content.extend_from_slice(&segment.sequence);
            dict_lengths.push(segment.sequence.len());
        }
        dict_content.push(b'$');
        dict_lengths.push(1);

        let n = dict_content.len();
        let sa  = SuffixArray::create(&*dict_content);
        let lcp = lcp_array(&dict_content, &sa).decompress();
        let isa = inverse_suffix_array(&sa);

        let mut nodeid = vec![0; n];
        let mut nodepos = vec![0; n];
        let mut id = 0;
        let mut pos = 0;
        for i in 0..n { 
            nodeid[isa[i]] = id;
            nodepos[isa[i]] = pos;

            pos += 1;
            if pos == dict_lengths[id] {
                id += 1;
                pos = 0;
            }
        }

        for i in 0..n {
            println!("{}\t{}\t{}\t{}\t{}\t{}", i, nodeid[i], nodepos[i], sa[i], lcp[i],
                str::from_utf8(&dict_content[sa[i]..]).unwrap())
        }
        println!();

        let m = &gfa.segments.len() + 1; // one fake sentinel segment
        let node_len = dict_lengths;

        let mut paths: Vec<Vec<u64>> = Vec::new();
        for path in &gfa.paths {
            let p = str::from_utf8(&path.segment_names).unwrap();
            let to_number = |x: &str| (x[0..x.len()-1].parse::<u64>()
                .expect("Cannot parse path."));

            let p: Vec<_> = p.split(',').map(to_number).collect();
            paths.push(p);
        }

        let mut node_freq = vec![0; m];
        for p in paths.iter() {
            for node in p {
                node_freq[*node as usize] += 1;
            }
            let last = node_freq.last_mut().unwrap();
            *last += 1;
        }

        let mut str_pos = vec![Vec::new(); m];
        let mut pos = 0;
        for p in paths.iter() {
            for node in p {
                str_pos[*node as usize].push(pos);
                pos += 1;
            }
            let last = str_pos.last_mut().unwrap();
            last.push(pos);
            pos += 1;
        }

        let mut orig_pos = vec![Vec::new(); m];
        let mut pos = 0;
        for p in paths.iter() {
            for node in p {
                orig_pos[*node as usize].push(pos);
                pos += node_len[*node as usize] - 1;
            }
            orig_pos[m-1].push(pos);
            pos += 1;
        }

        let mut path = Vec::new();
        for p in paths.iter() {
            let p: Vec<u64> = p.iter().map(|x| x + 1).collect();
            path.extend_from_slice(&p);
            path.push(0);
        }

        let path_sa = SuffixArray::create(&path[..]);
        let path_isa = inverse_suffix_array(&path_sa);
        let mut right_content_sa = Vec::new();
        for s_pos in &str_pos {
            let i = right_content_sa.len();
            right_content_sa.push(Vec::new());
            for p in s_pos {
                let rc_start = (p+1)%path_isa.len();
                right_content_sa[i].push(path_isa[rc_start]);
            }
        }

        for i in 0..m {
            println!("{}\t{}\t{}\t{:?}\t{:?}\t{:?}",
                i, node_len[i], node_freq[i], str_pos[i], orig_pos[i], right_content_sa[i]
            );
        }

        // https://github.com/rust-bio/rust-bio/issues/3
        // SAIS is generic, but suffix_array is not... WTF

    }
}

