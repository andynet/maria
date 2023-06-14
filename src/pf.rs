use bio::data_structures::suffix_array::lcp as lcp_array;
use gfa::gfa::GFA;
use std::io::Write;
use std::io;
use std::fs;
use std::collections::HashMap;
use gfa::gfa::Path;
use gfa::gfa::Segment;
use gfa::parser::GFAParser;
use crate::arrays::SuffixArray;
use std::ops::Add;
use std::str;
use crate::reverse_complement;

/// Data required to iterate through suffix array in sublinear space
pub struct PFData {
    segment_join: Vec<u8>,
    /// suffix array of segment join
     sa: Vec<usize>, 
    /// longest common prefix of segment join
    lcp: Vec<isize>,
    /// id of segments associated with suffix array
     id: Vec<usize>,
    /// position in the segment associated with suffix array
    pos: Vec<usize>,

    path_join: Vec<usize>,
    /// length of segments
    seg_len: Vec<usize>,
    /// starting positions of a segment in expanded path join
    seq_pos: Vec<Vec<usize>>,
    /// rank of right contexts of path join
    rc_rank: Vec<Vec<usize>>,

    /// segment overlap
    overlap: usize
}

pub struct PFDataIterator<'a> {
    data: &'a PFData,
    block: Block<'a>,
}

impl PFData {
    pub fn new(segments: &[Vec<u8>], paths: &[Vec<usize>], overlap: usize) -> Self {
        let segment_join = join(segments);

        let sa  = SuffixArray::create(&*segment_join);
        let lcp = lcp_array(&segment_join, &sa).decompress();
        let isa = permutation_invert(&sa);
        let id  = get_node_ids(&segment_join, &isa);
        let pos = get_node_pos(&segment_join, &isa);

        let path_join = join(paths);

        let seg_len = get_lengths(segments);
        let mut rc_rank = get_right_context_rank(&path_join, segments.len());
        let mut seq_pos = get_sequence_position(&path_join, &seg_len, overlap);

        for i in 0..segments.len() {
            let iperm = argsort(&rc_rank[i]);
            rc_rank[i] = permutation_apply(&iperm, &rc_rank[i]);
            seq_pos[i] = permutation_apply(&iperm, &seq_pos[i]);
        }

        Self{
            segment_join, sa, lcp, id, pos,
            path_join, seg_len, seq_pos, rc_rank,
            overlap 
        }
    }

    pub fn from_pfgraph(filename: &str) -> Self {
        let parser: GFAParser<usize, ()> = GFAParser::new();
        let gfa = parser.parse_file(filename)
            .expect("Error parsing GFA file.");

        let overlap = 1; // TODO determine_overlap();
        let segments: Vec<Vec<u8>> = parse_segments(&gfa.segments);
        let paths: Vec<Vec<usize>> = parse_paths(&gfa.paths);

        Self::new(&segments, &paths, overlap)
    }

    pub fn iter(&self) -> PFDataIterator {
        let block = Block::get_block_at(self, 0)
            .expect("No block found.");

        return PFDataIterator { 
            data: self,
            block
        }
    }

    pub fn print(&self) {
        let mut segment_join_repr = self.segment_join.clone();
        for c in segment_join_repr.iter_mut() {
            match *c {
                0   => { *c = b'$'; },
                1   => { *c = b'#'; },
                2.. => { *c -= 2; }
            }
        }
        for i in 0..self.sa.len() {
            println!("{}\t{}\t{}\t{}\t{}\t{}", 
                i, self.sa[i], self.lcp[i], self.id[i], self.pos[i],
                str::from_utf8(&segment_join_repr[self.sa[i]..]).unwrap()
            );
        }
        println!();

        for id in 0..self.seg_len.len() {
            println!("{}\t{}\t{:?}\t{:?}", 
                id, self.seg_len[id], self.seq_pos[id], self.rc_rank[id]
            );
        }
        println!();
    }
}

impl<'a> Iterator for PFDataIterator<'a> {
    type Item = (usize, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let block = &mut self.block;
        if let Some(x) = block.next() { return Some(x); }

        let tmp = Block::get_block_at(self.data, block.end);
        if tmp.is_none() { return None; }
        *block = tmp.unwrap();
        return block.next();
    }
}

struct Block<'a> {
     data: &'a PFData,
    start: usize,       // start of a block (inclusive)
      end: usize,       //   end of a block (exclusive)
      idx: Vec<usize>,  // index of smallest rank at position start + i
     rank: Vec<usize>,  // smalles rank at position start + i
}

impl<'a> Block<'a> {
    fn get_block_at(data: &'a PFData, mut start: usize) -> Option<Block<'a>> {
        if start >= data.sa.len() { return None; }

        let mut remaining = data.seg_len[data.id[start]] - data.pos[start];
        while remaining <= data.overlap { 
            start += 1; 
            remaining = data.seg_len[data.id[start]] - data.pos[start];
        }

        let mut end = start + 1;
        while data.lcp[end] >= remaining as isize {
            end += 1;
        }

        let idx = vec![0; end-start];

        let mut rank = Vec::new();
        for id in &data.id[start..end] {
            rank.push(data.rc_rank[*id][0]);
        }

        return Some(Block{data, start, end, idx, rank});
    }
}

impl<'a> Iterator for Block<'a> {
    type Item = (usize, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let k = argmin(&self.rank);
        if self.rank[k] == usize::MAX { return None; }

        let data = self.data;
        let id = data.id[self.start + k];
        let pos = data.pos[self.start + k];

        let result = (data.seq_pos[id][self.idx[k]] + pos, id, pos);

        self.idx[k] += 1;
        match data.rc_rank[data.id[self.start + k]].get(self.idx[k]) {
            Some(val) => { self.rank[k] = *val; },
            None      => { self.rank[k] = usize::MAX; }
        }
        return Some(result);
    }
}

fn parse_segments(segments: &[Segment<usize, ()>]) -> Vec<Vec<u8>> {
    let mut result = Vec::new();
    for segment in segments {
        result.push(segment.sequence.clone());
    }
    return result;
}

fn parse_paths(paths: &[Path<usize, ()>]) -> Vec<Vec<usize>> {
    let mut result = Vec::new();
    for path in paths {
        let p = str::from_utf8(&path.segment_names).unwrap();
        let to_usize = |x: &str| (x[0..x.len()-1].parse::<usize>()
            .expect("Cannot parse path."));

        let p: Vec<_> = p.split(',').map(to_usize).collect();
        result.push(p);
    }
    return result;
}

fn join<T>(slice: &[Vec<T>]) -> Vec<T> 
where
    T: From<u8> + Add<T, Output = T> + Copy
{
    let mut result = Vec::new();
    for vector in slice {
        let v: Vec<T> = vector.iter().map(|&x| x + 2.into()).collect();
        result.extend_from_slice(&v);
        result.push(1.into());
    }
    result.push(0.into());
    return result;
}

fn get_lengths(segments: &[Vec<u8>]) -> Vec<usize> { 
    let mut len = Vec::new();
    for segment in segments {
        len.push(segment.len());
    }
    len.push(1);
    return len;
}

fn get_node_ids(segment_join: &[u8], isa: &[usize]) -> Vec<usize> {
    let mut result = vec![0; segment_join.len()];
    let mut id = 0;
    for i in 0..segment_join.len() {
        result[isa[i]] = id;
        if segment_join[i] == 1 { id += 1; }
    }
    return result;
}

fn get_node_pos(segment_join: &[u8], isa: &[usize]) -> Vec<usize> {
    let mut result = vec![0; segment_join.len()];
    let mut pos = 0;
    for i in 0..segment_join.len() {
        result[isa[i]] = pos;
        pos += 1;
        if segment_join[i] == 1 { pos = 0; }
    }
    return result;
}

fn get_sequence_position(path_join: &[usize], len: &[usize], overlap: usize) -> Vec<Vec<usize>> {
    let mut result = vec![Vec::new(); len.len()];
    let mut start = 0;
    for id in path_join {
        if *id >= 2 {
            result[id-2].push(start);
            start += len[id-2] - overlap;
        }
    }
    return result;
}

fn get_right_context_rank(path_join: &[usize], size: usize) -> Vec<Vec<usize>> {
    let mut result = vec![Vec::new(); size];

    let sa = SuffixArray::create(&path_join[..]);
    let isa = permutation_invert(&sa);

    for (i, id) in path_join.iter().enumerate() {
        if *id >= 2 {
            result[id-2].push(isa[i+1]);
        }
    }
    return result;
}

fn permutation_invert(perm: &[usize]) -> Vec<usize> {
    let mut inverse = vec![0; perm.len()];
    for i in 0..perm.len() { inverse[perm[i]] = i; }
    return inverse;
}

fn permutation_apply(perm: &[usize], data: &[usize]) -> Vec<usize> {
    let mut result = vec![0; perm.len()];
    for i in 0..perm.len() {
        result[perm[i]] = data[i];
    }
    return result;
}

fn argsort(data: &[usize]) -> Vec<usize> {
    let mut v: Vec<(usize, usize)> = data.iter().enumerate().map(|(i, x)| (*x, i)).collect();
    v.sort_unstable();
    let v: Vec<usize> = v.iter().map(|(_, i)| *i).collect();
    let v = permutation_invert(&v);
    return v;
}

fn argmin(data: &[usize]) -> usize {
    let mut min_pos = 0;
    let mut min_val = usize::MAX;
    for (p, &v) in data.iter().enumerate() {
        if v < min_val { min_pos = p; min_val = v; }
    }
    return min_pos;
}

pub fn get_triggers(trigs: &[u8], size: usize) -> Vec<&[u8]> {
    let mut result = Vec::new();
    for i in (0..trigs.len()).step_by(size+1) {
        result.push(&trigs[i..i+size]);
    }
    return result;
}

pub fn load_trigs(filename: &str) -> (Vec<u8>, usize) {
    let trigs = fs::read_to_string(filename)
        .expect("Unable to read the triggers file")
        .trim().as_bytes().to_owned();

    let trigs_size;
    match trigs.iter().position(|&x| x == b'\n') {
        None    => { trigs_size = trigs.len(); },
        Some(x) => { trigs_size = x; }
    }
    return (trigs, trigs_size);
}

pub fn normalize(segments: HashMap<Vec<u8>, usize>, mut paths: Vec<Vec<usize>>)
    -> (Vec<Vec<u8>>, Vec<Vec<usize>>) 
{
    let mut segments: Vec<_> = segments.into_iter().collect();
    segments.sort_unstable();

    let mut mapping = vec![0; segments.len()];
    for (i, (_, id)) in segments.iter().enumerate() { mapping[*id] = i; }

    for (_, id) in segments.iter_mut() { *id = mapping[*id]; }
    for path in paths.iter_mut() {
        for id in path { *id = mapping[*id]; }
    }

    let segments = segments.iter().map(|x| x.0.clone() ).collect();

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

pub fn split_prefix_free(
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

fn into_path_step(step: &str) -> (usize, u8) {
    let id = step[0..step.len()-1].parse::<usize>().expect("Cannot parse path.");
    let sign = step.bytes().last().unwrap();
    return (id, sign);
}

pub fn reconstruct_path(path: &Path<usize, ()>, gfa: &GFA<usize, ()>) -> Vec<u8> {
    let mut result = Vec::new();

    let mut map = HashMap::new();
    for segment in &gfa.segments {
        map.insert(segment.name, segment.sequence.clone());
    }

    let p = str::from_utf8(&path.segment_names).unwrap();
    let p: Vec<_> = p.split(',').map(into_path_step).collect();

    for (id, sign) in p {
        let sequence = map.get(&id).expect("Segment not found");
        match sign {
            b'+' => { result.extend_from_slice(&sequence) },
            b'-' => { result.extend_from_slice(&reverse_complement(&sequence)) },
            x    => { panic!("Unexpected direction {x}") }
        }
    } 
    return result;
}

pub fn print_gfa<T: Write>(
      segments: &Vec<Vec<u8>>,
         paths: &Vec<Vec<usize>>,
             k: usize,  // size of the trigger words
    mut output: T
) -> io::Result<()> {

    writeln!(output, "H\tVN:Z:1.1")?;
    for (id, seq) in segments.iter().enumerate() {
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
    use super::into_path_step;
    #[test]
    fn into_iter_works() {
        let s = "2+";
        let res = into_path_step(s);
        assert_eq!(res, (2, b'+'));
    }

}
