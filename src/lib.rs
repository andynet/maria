use gfa::gfa::{Path, SegmentId, Segment};
use gfa::optfields::*;
use std::collections::HashMap;
use std::hash::Hash;

fn parse_segments(s: &[u8]) -> Vec<(Vec<u8>, u8)> {
    let mut start = 0;
    let mut segments = Vec::new();
    for i in 0..s.len() {
        if s[i] == b',' {
            let name = s[start..i-1].to_vec();
            let sign = s[i-1];
            segments.push((name, sign));
            start = i+1;
        }
    }
    let name = s[start..s.len()-1].to_vec();
    let sign = s[s.len()-1];
    segments.push((name, sign));
    return segments;
}

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    let n = s.len();
    let mut result = vec![0; n];
    for (i, c) in s.iter().enumerate() {
        match *c {
            b'A' => { result[n-1-i] = b'T'; },
            b'C' => { result[n-1-i] = b'G'; },
            b'G' => { result[n-1-i] = b'C'; },
            b'T' => { result[n-1-i] = b'A'; },
            b'N' => { result[n-1-i] = b'N'; },
            _ => { panic!("Unexpected letter in s."); }
        }
    }
    return result;
}

fn get_sequence(
    segments: &[(Vec<u8>, u8)], map: &HashMap<Vec<u8>, Vec<u8>>
) -> Vec<u8> {
    let mut sequence = Vec::new();
    for (name, sign) in segments {
        let node_string = map.get(name).unwrap();
        if *sign == b'+' {
            sequence.extend_from_slice(node_string);
        } else {
            let rc = reverse_complement(node_string);
            sequence.extend_from_slice(&rc);
        }
    }
    return sequence;
}

pub fn inverse_suffix_array(sa: &[usize]) -> Vec<usize> {
    let mut isa = vec![0; sa.len()];
    for i in 0..sa.len() {
        isa[sa[i]] = i;
    }
    return isa;
}

pub fn tag_array<N: SegmentId, T: OptFields>(
    sa: &[usize],
    isa: &[usize],
    paths: &[Path<N, T>],
    map: &HashMap<Vec<u8>, Vec<u8>>
) -> Vec<(Vec<u8>, usize)> {

    let mut tag = vec![(Vec::new(), 0); sa.len()];
    let mut j = 0;
    for path in paths {
        let path = parse_segments(&path.segment_names);
        for (node_name, sign) in path {
            let node = map.get(&node_name).unwrap();
            for i in 0..node.len() {
                let mut v = node_name.clone();
                v.push(sign);
                tag[isa[j]] = (v, i);
                j += 1;
            }
        }
        j += 1;
    }
    return tag;
}

pub fn create_sequence<N, T>(paths: &[Path<N, T>], map: &HashMap<Vec<u8>, Vec<u8>>) 
    -> (Vec<u8>, Vec<usize>)
where
    N: SegmentId,
    T: OptFields
{
    let mut seq = Vec::new();
    let mut lengths = Vec::new();
    for p in paths {
        let segments = parse_segments(&p.segment_names);
        let sequence = get_sequence(&segments, &map);

        seq.extend_from_slice(&sequence);
        lengths.push(seq.len());
        seq.push(b'$');
    }
    return (seq, lengths);
}

pub fn create_map<N, T>(segments: &[Segment<N, T>]) -> HashMap<N, Vec<u8>>
where
    N: SegmentId + Clone + Eq + Hash,
    T: OptFields + Clone
{
    let mut map = HashMap::new();
    for s in segments { 
        map.insert(s.name.clone(), s.sequence.clone()); 
    }
    return map;
}

#[cfg(test)]
mod tests {
    use super::*;
    use bio::data_structures::suffix_array::suffix_array;
    use gfa::parser::GFAParser;
    use std::{str, usize};

    fn print_table(seq: &[u8], sa: &[usize], tag: &[(Vec<u8>, usize)]) {
        for i in 0..sa.len() {
            print!("{}\t", sa[i]);
            print!("{}\t", str::from_utf8(&tag[i].0).unwrap());
            print!("{}\t", tag[i].1);
            print!("{}", str::from_utf8(&seq[sa[i]..]).unwrap());
            print!("{}\n", str::from_utf8(&seq[..sa[i]]).unwrap());
        }
    }

    #[test]
    fn test() {
        const FILE: &str = "data/test_gfa.gfa";

        let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();
        let gfa = parser.parse_file(FILE).expect("Error parsing file.");
        let version = gfa.header.version.unwrap();
        println!("{}", str::from_utf8(&version).unwrap());

        let map = create_map(&gfa.segments);
        let (seq, _) = create_sequence(&gfa.paths, &map);
        println!("{}", str::from_utf8(&seq).unwrap());

        let sa  = suffix_array(&seq);
        let isa = inverse_suffix_array(&sa);
        let tag = tag_array(&sa, &isa, &gfa.paths, &map);
        println!("{:?}", tag);
        print_table(&seq, &sa, &tag);
    }

    #[test]
    fn test_small() {
        const FILE: &str = "data/test_small_gfa.gfa";

        let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();
        let gfa = parser.parse_file(FILE).expect("Error parsing file.");

        let map = create_map(&gfa.segments);
        let (seq, _) = create_sequence(&gfa.paths, &map);

        let sa  = suffix_array(&seq);
        let isa = inverse_suffix_array(&sa);
        let tag = tag_array(&sa, &isa, &gfa.paths, &map);
        println!("{:?}", tag);
        let expected_tag: Vec<(Vec<u8>, usize)> = vec![
            (  b"".to_vec(), 0), (  b"".to_vec(), 0), (b"2+".to_vec(), 1), (b"2-".to_vec(), 1),
            (b"4+".to_vec(), 0), (b"4+".to_vec(), 0), (b"1-".to_vec(), 0), (b"0+".to_vec(), 0),
            (b"0+".to_vec(), 0), (b"2+".to_vec(), 1), (b"4+".to_vec(), 1), (b"4+".to_vec(), 1),
            (b"1-".to_vec(), 1), (b"0+".to_vec(), 1), (b"0+".to_vec(), 1), (b"3+".to_vec(), 0),
            (b"3+".to_vec(), 0), (b"1+".to_vec(), 0), (b"2+".to_vec(), 0), (b"2-".to_vec(), 0),
            (b"3+".to_vec(), 1), (b"2+".to_vec(), 0), (b"3+".to_vec(), 1), (b"1+".to_vec(), 1)
        ];

        assert_eq!(tag, expected_tag);
    }
}
