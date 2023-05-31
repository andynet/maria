use std::collections::HashMap;
use std::str;

fn main() { todo!(); }

fn split_prefix_free(seq: &[u8], triggers: &[&[u8]]) 
    -> (HashMap<Vec<u8>, usize>, Vec<Vec<usize>>) 
{
    let mut segments = HashMap::new();
    let mut paths = Vec::new();

    return (segments, paths);
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test() {
        let triggers: Vec<&[u8]> = vec![b"$$", b"AC"];
        let seq = b"$$AGACACGTACGAAC$$";

        let mut segments: HashMap<Vec<u8>, usize> = HashMap::new();
        let mut paths: Vec<Vec<usize>> = Vec::new();

        let k = triggers[0].len();

        // create phrases
        let mut start = 0;
        let mut segid = 0;
        paths.push(Vec::new());
        let n = seq.len();
        for i in 0..n-k+1 {
            if triggers.contains(&&seq[i..i+k]) {
                println!("{}", str::from_utf8(&seq[start..i+k]).unwrap());
                let v = segments.get(&seq[start..k+i]);
                match v {
                    None => {
                        segments.insert(seq[start..i+k].to_owned(), segid);
                        paths.last_mut().unwrap().push(segid);
                        segid += 1;
                    },
                    Some(v) => { paths.last_mut().unwrap().push(*v); }
                }
                start = i;
            }
        }
        println!("{:?}", segments);
        println!("{:?}", paths);

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
