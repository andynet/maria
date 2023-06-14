use crate::grammar::Grammar;
use std::{process::Command, fs::File};
use std::str;
use std::io::Write;

#[test]
fn test() {
    let grammar = Grammar::from_file("data/pftag/test_join.txt.plainslp");

    for i in 0..grammar.len() {
        print!("{}", grammar[i] as char);
    }
    println!();
}

#[test]
fn complex_test() {
    let s = b"AGCTTGCGTAGCTAGCTGAGCTGATCG".to_owned();
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
