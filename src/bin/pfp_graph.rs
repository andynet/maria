use bio::io::fasta;
use clap::Parser;
use std::collections::HashMap;
use std::fs::{self, File};
use std::io::{self, Write};
use std::str;

/// Build prefix-free graph
#[derive(Parser, Debug)]
#[command(author, version, about)]
struct Args {
    /// fasta file
    fasta_file: String,
    /// trigger file, containing one trigger per line
    trigger_file: String,
    /// output in GFA format v1.1
    output: String,
}

fn main() { 
    let args = Args::parse();

    let (trigs, trigs_size) = load_trigs(&args.trigger_file);
    let triggers = get_triggers(&trigs, trigs_size);

    let mut segments = HashMap::new();
    let mut paths = Vec::new();

    let fasta = File::open(args.fasta_file).expect("Cannot open fasta file.");
    let mut records = fasta::Reader::new(fasta).records();
    while let Some(Ok(record)) = records.next() {
        // record.seq().len();
        let mut seq = record.seq().to_owned();
        seq.push(b'$');
        split_prefix_free(&seq, &triggers, &mut segments, &mut paths);
    }

    let (segments, paths) = normalize(segments, paths);

    let output = File::create(args.output).expect("Cannot create output file.");
    print_gfa(&segments, &paths, trigs_size, output)
        .expect("Error writting GFA");
}

fn get_triggers(trigs: &[u8], size: usize) -> Vec<&[u8]> {
    let mut result = Vec::new();
    for i in (0..trigs.len()).step_by(size+1) {
        result.push(&trigs[i..i+size]);
    }
    return result;
}

fn load_trigs(filename: &str) -> (Vec<u8>, usize) {
    let trigs = fs::read_to_string(filename)
        .expect("Unable to read the triggers file")
        .trim().as_bytes().to_owned();

    let trigs_size = trigs.iter().position(|&x| x == b'\n').unwrap();
    return (trigs, trigs_size);
}

fn normalize(segments: HashMap<Vec<u8>, usize>, mut paths: Vec<Vec<usize>>)
    -> (Vec<(Vec<u8>, usize)>, Vec<Vec<usize>>) 
{
    let mut segments: Vec<_> = segments.into_iter().collect();
    segments.sort_unstable();

    let mut mapping = vec![0; segments.len()];
    for (i, (_, id)) in segments.iter().enumerate() { mapping[*id] = i; }

    for (_, id) in segments.iter_mut() { *id = mapping[*id]; }
    for path in paths.iter_mut() {
        for id in path { *id = mapping[*id]; }
    }
    return (segments, paths);
}

fn add_segment(
         seq: &[u8],
    segments: &mut HashMap<Vec<u8>, usize>,
        path: &mut Vec<usize>
) {
    let segment_id = segments.get(seq);
    match segment_id {
        Some(&id) => { path.push(id); },
        None => {
            path.push(segments.len());
            segments.insert(seq.to_owned(), segments.len());
        }
    }
}

fn split_prefix_free(
         seq: &[u8],    // must end with sentinel
    triggers: &[&[u8]], // must be non-empty
    segments: &mut HashMap<Vec<u8>, usize>,
       paths: &mut Vec<Vec<usize>>
) {
    let n = seq.len();
    let k = triggers.first().expect("No triggers found.").len();

    let mut path = Vec::new();
    let mut i = 0;
    for j in 1..n-k {
        if triggers.contains(&&seq[j..j+k]) {
            let segment_seq = &seq[i..j+k];
            add_segment(segment_seq, segments, &mut path);
            i = j;
        }
    }
    add_segment(&seq[i..n], segments, &mut path);
    paths.push(path);
}

fn print_gfa<T: Write>(
      segments: &Vec<(Vec<u8>, usize)>,
         paths: &Vec<Vec<usize>>,
             k: usize,  // size of the trigger words
    mut output: T
) -> io::Result<()> {

    writeln!(output, "H\tVN:Z:1.1")?;
    for (seq, id) in segments {
        let seq = str::from_utf8(seq).expect("Cannot convert seq to UTF8");
        writeln!(output, "S\t{}\t{}", id, seq)?;
    }

    for (i, path) in paths.iter().enumerate() {
        let path_str = path.iter().map(|x| format!("{}+", x)).collect::<Vec<_>>().join(",");
        writeln!(output, "P\t{}\t{}\t*", i, path_str)?;

        for j in 0..path.len()-1 {
            writeln!(output, "L\t{}\t+\t{}\t+\t{}M", path[j], path[j+1], k)?;
        }
    }
    return Ok(());
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::stdout;

    #[test]
    fn test_only_prefix_free() {
        let triggers: Vec<&[u8]> = vec![b"T"];
        let k = 1;
        let seq1 = b"ATCTGTTAATG$";
        let seq2 = b"AACGTGTACGTACGAAC$";

        let mut segments = HashMap::new();
        let mut paths = Vec::new();

        split_prefix_free(seq1, &triggers, &mut segments, &mut paths);
        split_prefix_free(seq2, &triggers, &mut segments, &mut paths);

        let (segments, paths) = normalize(segments, paths);

        print_gfa(&segments, &paths, k, stdout()).expect("Error writting GFA");
    }
}
