use std::env;
use gfa::parser::GFAParser;
use gfa::optfields::*;
use maria::{create_map,create_sequence,inverse_suffix_array,tag_array};
use bio::data_structures::suffix_array::suffix_array;
use std::str;
use std::io::{BufWriter, stdout, Stdout, Write};

fn main() {
    let args: Vec<String> = env::args().collect();
    if args.len() != 2 {
        println!{"Usage: {} input.gfa > output.txt", args[0]};
        return;
    }
    let file_name = &args[1];

    let parser: GFAParser<Vec<u8>, Vec<OptField>> = GFAParser::new();
    let gfa = parser.parse_file(file_name).expect("Error parsing file.");

    let map = create_map(&gfa.segments);
    let seq = create_sequence(&gfa.paths, &map);

    let sa  = suffix_array(&seq);
    let isa = inverse_suffix_array(&sa);
    let tag = tag_array(&sa, &isa, &gfa.paths, &map);
 
    let buffer = BufWriter::new(stdout());
    write_table(&sa, &tag, None, buffer).expect("Error writing file.");
}

fn write_table(
    sa: &[usize], tag: &[(Vec<u8>, usize)], seq: Option<&[u8]>, 
    mut output: BufWriter<Stdout>
) -> std::io::Result<()> {
    for i in 0..sa.len() {
        output.write(format!("{}\t", sa[i]).as_bytes())?;
        output.write(format!("{}\t", str::from_utf8(&tag[i].0).unwrap()).as_bytes())?;
        output.write(format!("{}\t", tag[i].1).as_bytes())?;

        if let Some(s) = seq {
            output.write(format!("{}", str::from_utf8(&s[sa[i]..]).unwrap()).as_bytes())?;
            output.write(format!("{}", str::from_utf8(&s[..sa[i]]).unwrap()).as_bytes())?;
        }
        output.write(b"\n")?;
    }
    Ok(())
}

