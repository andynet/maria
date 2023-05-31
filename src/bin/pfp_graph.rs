use std::collections::HashMap;
use std::str;

fn main() { todo!(); }

fn split_prefix_free(
         seq: &[u8],
    triggers: &[&[u8]],
    segments: &mut HashMap<Vec<u8>, usize>,
       paths: &mut Vec<Vec<usize>>
) {
    let n = seq.len();
    let k = triggers.first().expect("No triggers found.").len();

    let mut path = Vec::new();
    let mut i = 0;
    for j in 0..n-k+1 {
        if triggers.contains(&&seq[j..j+k]) {
            let segment_seq = &seq[i..j+k];
            let  segment_id = segments.get(segment_seq);

            match segment_id {
                Some(&id) => { path.push(id); },
                None => {
                    path.push(segments.len());
                    segments.insert(segment_seq.to_owned(), segments.len());
                }
            }

            i = j;
        }
    }
    paths.push(path);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let triggers: Vec<&[u8]> = vec![b"$$", b"AC"];
        let seq1 = b"$$AGACACGTACGAAC$$";
        let seq2 = b"$$AACGTGTACGTACGAAC$$";

        let mut segments = HashMap::new();
        let mut paths    = Vec::new();

        split_prefix_free(seq1, &triggers, &mut segments, &mut paths);
        split_prefix_free(seq2, &triggers, &mut segments, &mut paths);

        println!("{:?}", segments);
        println!("{:?}", paths);
        for (k, v) in segments.iter() {
            println!("{}: {}", v, str::from_utf8(k).unwrap());
        }

        // reorder phrases
        let mut map: Vec<(&Vec<u8>, &usize)> = segments.iter().collect();
        map.sort_unstable();
        println!("{:?}", map);

        let mut mapping = vec![0; map.len()];
        for (i, (_, &v)) in map.iter().enumerate() {
            mapping[v] = i;
        }
        println!("{:?}", mapping);

        for (_, v) in segments.iter_mut() { *v = mapping[*v]; }
        for path in &mut paths {
            println!("{:?}", path);
            for i in 0..path.len() {
                path[i] = mapping[path[i]];
            }
        }
        println!("{:?}", segments);
        println!("{:?}", paths);

        // print GFA
        println!("H\tVN:Z:1.1");
        for (k, v) in segments {
            println!("S\t{}\t{}", v, str::from_utf8(&k[2..]).unwrap())
        }
        for (i, path) in paths.iter().enumerate() {
            print!("P\t{}\t", i);
            for val in path { print!("{}+,", val); }
            println!("\t*");

            for i in 1..path.len() {
                println!("L\t{}\t+\t{}\t+\t0M", path[i-1], path[i]);
            }
        }
    }
}
