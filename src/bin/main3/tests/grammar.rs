use crate::grammar::Grammar;
use std::io::{BufReader, BufRead};
use std::{process::Command, fs::File};
use std::str;
use std::io::Write;
use proptest::prelude::*;

#[test]
fn test() {
    let grammar = Grammar::from_file("data/pftag/test_join.txt.plainslp");

    for i in 0..grammar.len() {
        print!("{}", grammar[i] as char);
    }
    println!();
}

#[test]
fn test_my_grammar() {
    let s = b"TGACGGGCAGT".to_owned();
    let g = Grammar::from_bytes(&s);

    assert_eq!(g.len(), s.len());
    for i in 0..g.len() {
        assert_eq!(g[i], s[i]);
    }

    let output = File::create("data/temporary/seq.txt").expect("Cannot create file.");
    g.print(output);
}

#[test]
fn generate_real() {
    let input = File::open("data/real/SARS-CoV2.5.fnajoin")
        .expect("Cannot open input file.");
    let line = BufReader::new(input).lines().next().unwrap().unwrap();
    let g = Grammar::from_bytes(line.as_bytes());
    let output = File::create("data/real/SARS-CoV2.5.fnajoin.plainslp")
        .expect("Cannot open output file");
    g.print(output);
}

fn test_grammar_for_string(s: &[u8]) {
    {
        let mut out = File::create("data/temporary/seq.txt").expect("Cannot create file.");
        write!(out, "{}", str::from_utf8(&s).unwrap()).expect("Cannot write to file.");
    }

    println!("Building grammar...");
    Command::new("./tools/bigrepair/bigrepair")
        .args(["data/temporary/seq.txt"])
        .output()
        .expect("Failer to run bigrepair");

    println!("Converting rules to plaintext...");
    Command::new("./scripts/print_plain_slp")
        .args(["data/temporary/seq.txt"])
        .output()
        .expect("Failed to convert to plain text");

    println!("Loading grammar...");
    let g = Grammar::from_file("data/temporary/seq.txt.plainslp");

    println!("Testing lengths...");
    assert_eq!(g.len(), s.len());
    println!("Testing characters...");
    for i in 0..g.len() {
        assert_eq!(g[i], s[i]);
    }
}

#[test]
fn complex_test() {
    // let s = b"AGCTTGCGTAGCTAGCTGAGCTGATCG".to_owned();
    let s = b"TGACGGGCAGT".to_owned();
    test_grammar_for_string(&s);
}

proptest! {
    #[test]
    fn proptesting(s in "[ACGT]+") {
        let s: Vec<_> = s.bytes().collect();
        test_grammar_for_string(&s);
    }
}
