use gfa::parser::GFAParser;
use gfa::optfields::*;
use maria::{create_map,create_sequence,inverse_suffix_array,tag_array};
use bio::data_structures::suffix_array::suffix_array;
use std::str;
use clap::Parser;
use std::fmt;

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

#[derive(Debug)]
struct GraphPos {
    node_id: Vec<u8>,
    pos: usize,
}

impl fmt::Display for GraphPos {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(f, "({}, {})", str::from_utf8(&self.node_id).unwrap(), self.pos)
    }
}

impl From<(Vec<u8>, usize)> for GraphPos {
    fn from(value: (Vec<u8>, usize)) -> Self {
        GraphPos{ node_id: value.0, pos: value.1 }
    }
}

fn print_tag(tag: &[GraphPos]) {
    print!("{}: ", tag.len());
    for i in 0..tag.len() { print!("{} ", tag[i]); }
    println!();
}

fn get_compressed_tag_array(gfa: &str) -> Vec<GraphPos> {
    let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();
    let gfa = parser.parse_file(gfa).expect("Error parsing file.");

    let map = create_map(&gfa.segments);
    let seq = create_sequence(&gfa.paths, &map);

    let sa  = suffix_array(&seq);
    let isa = inverse_suffix_array(&sa);
    let tag = tag_array(&sa, &isa, &gfa.paths, &map);

    let mut result = Vec::new();
    result.push(tag[0].clone().into());
    for i in 1..tag.len() {
        if tag[i] != tag[i-1] {
            result.push(tag[i-1].clone().into());
            result.push(tag[i].clone().into());
        }
    }
    result.push(tag[tag.len()-1].clone().into());
    return result;
}

// fn get_mems(mems: &str, ptrs: &str) -> Vec<(usize, usize)> {
//     let mut file = File::open(mems).expect("Cannot read MEM file.");
//     return Vec::new();
// }

fn main() {
    let args = Args::parse();
    let tag = get_compressed_tag_array(&args.gfa_filename);
    print_tag(&tag[..10]);

    // let mems = get_mems(&args.mems_filename, &args.ptr_filename);
}


