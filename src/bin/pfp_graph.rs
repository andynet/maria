use bio::io::fasta;
use clap::Parser;
use std::collections::HashMap;
use std::fs::File;
use maria::pf;

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

    let (trigs, trigs_size) = pf::load_trigs(&args.trigger_file);
    let triggers = pf::get_triggers(&trigs, trigs_size);

    let mut segments = HashMap::new();
    let mut paths = Vec::new();

    let fasta = File::open(args.fasta_file).expect("Cannot open fasta file.");
    let mut records = fasta::Reader::new(fasta).records();
    while let Some(Ok(record)) = records.next() {
        // record.seq().len();
        let mut seq = record.seq().to_owned();
        seq.push(b'.');
        pf::split_prefix_free(&seq, &triggers, &mut segments, &mut paths);
    }

    let (segments, paths) = pf::normalize(segments, paths);

    let output = File::create(args.output).expect("Cannot create output file.");
    pf::print_gfa(&segments, &paths, trigs_size, output)
        .expect("Error writting GFA");
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

        pf::split_prefix_free(seq1, &triggers, &mut segments, &mut paths);
        pf::split_prefix_free(seq2, &triggers, &mut segments, &mut paths);

        let (segments, paths) = pf::normalize(segments, paths);

        pf::print_gfa(&segments, &paths, k, stdout()).expect("Error writting GFA");
    }
}
