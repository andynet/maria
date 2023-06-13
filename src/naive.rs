use gfa::gfa::{Path, SegmentId, Segment};
use gfa::optfields::*;
use std::collections::HashMap;
use std::hash::Hash;

use crate::reverse_complement;

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
    for i in 0..sa.len() { isa[sa[i]] = i; }
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


