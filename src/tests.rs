use super::naive::*;
use bio::data_structures::suffix_array::suffix_array;
use gfa::parser::GFAParser;
use gfa::optfields::OptField;
use std::{str, usize};

fn print_table(seq: &[u8], sa: &[usize], tag: &[(Vec<u8>, usize)]) {
    for i in 0..sa.len() {
        println!("{}\t{}\t{}\t{}{}",
            sa[i],
            str::from_utf8(&tag[i].0).unwrap(),
            tag[i].1,
            str::from_utf8(&seq[sa[i]..]).unwrap(),
            str::from_utf8(&seq[..sa[i]]).unwrap(),
        );
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

