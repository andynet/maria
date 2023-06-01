use std::collections::HashMap;
use std::io::stdout;
use std::io::Write;
use std::io;
use std::str;

fn main() { 
    let triggers: Vec<&[u8]> = vec![b"T"];
    let k = 1;
    let seq = b"ATCTGTTAATG$";

    let mut segments = HashMap::new();
    let mut paths = Vec::new();

    split_prefix_free(seq, &triggers, &mut segments, &mut paths);

    let (segments, paths) = normalize(segments, paths);

    print_gfa(&segments, &paths, k, stdout()).expect("Error writting GFA");
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
