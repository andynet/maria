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

        // let paths: Vec<Vec<u64>> = Vec::new();
        // for path in &gfa.paths {
        //     let p = str::from_utf8(&path.segment_names).unwrap();
        //     println!("{p}");
        //
        //     let to_number = |x: &str| (x[0..x.len()-1].parse::<u64>()
        //         .expect("Cannot parse path."));
        //
        //     let p: Vec<_> = p.split(',').map(to_number).collect();
        //     paths.push(p);
        // }

        // let v: Vec<u64> = vec![1215, 5468, 45];
        // https://github.com/rust-bio/rust-bio/issues/3
        // SAIS is generic, but suffix_array is not... WTF
        // let sa2 = SuffixArray::create(&v[..]);

    }
}

