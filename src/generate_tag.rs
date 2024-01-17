use gfa::parser::GFAParser;
use gfa::optfields::*;
use maria::naive::{create_map,create_sequence,inverse_suffix_array,tag_array};
use bio::data_structures::suffix_array::suffix_array;
use std::str;
use std::io::{BufWriter, stdout, Stdout, Write};
use clap::Parser;

/// Generate tag array from gfa file
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// GFA file
    file_name: String,

    /// Print the Wheeler matrix
    #[arg(short, action)]
    print: bool,
}

fn main() {
    let args = Args::parse();

    let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();
    let gfa = parser.parse_file(args.file_name).expect("Error parsing file.");

    let map = create_map(&gfa.segments);
    let (seq, _) = create_sequence(&gfa.paths, &map);

    let sa  = suffix_array(&seq);
    let isa = inverse_suffix_array(&sa);
    let tag = tag_array(&sa, &isa, &gfa.paths, &map);
 
    let mut buffer = BufWriter::new(stdout());
    let s = if args.print { Some(&seq[..]) } else { None };

    write_table(&sa, &tag, s, &mut buffer).expect("Error writting_file");
    buffer.flush().expect("Error writting file.");
}

fn write_table(
    sa: &[usize], tag: &[(Vec<u8>, usize)], seq: Option<&[u8]>, 
    output: &mut BufWriter<Stdout>
) -> std::io::Result<()> {
    for i in 0..sa.len() {
        output.write_all(format!("{}\t", sa[i]).as_bytes())?;
        output.write_all(format!("{}\t", str::from_utf8(&tag[i].0).unwrap()).as_bytes())?;
        output.write_all(format!("{}\t", tag[i].1).as_bytes())?;

        if let Some(s) = seq {
            output.write_all(str::from_utf8(&s[sa[i]..]).unwrap().as_bytes())?;
            output.write_all(str::from_utf8(&s[..sa[i]]).unwrap().as_bytes())?;
        }
        output.write_all(b"\n")?;
    }
    Ok(())
}

