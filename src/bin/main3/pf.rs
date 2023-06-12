use bio::data_structures::suffix_array::lcp as lcp_array;
use gfa::gfa::Path;
use gfa::gfa::Segment;
use gfa::parser::GFAParser;
use maria::arrays::SuffixArray;
use std::ops::Add;
use std::str;

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
    pub fn from_pfgraph(filename: &str) -> Self {
        let parser: GFAParser<usize, ()> = GFAParser::new();
        let gfa = parser.parse_file(filename)
            .expect("Error parsing GFA file.");

        let overlap = 1; // TODO
        let segments: Vec<Vec<u8>> = parse_segments(&gfa.segments);
        let paths: Vec<Vec<usize>> = parse_paths(&gfa.paths);

        let segment_join = join(&segments);

        let sa  = SuffixArray::create(&*segment_join);
        let lcp = lcp_array(&segment_join, &sa).decompress();
        let isa = permutation_invert(&sa);
        let id  = get_node_ids(&segment_join, &isa);
        let pos = get_node_pos(&segment_join, &isa);

        let path_join = join(&paths);

        let len = get_lengths(&segments);
        let mut rc_rank = get_right_context_rank(&path_join, segments.len());
        let mut seq_pos = get_sequence_position(&path_join, &len, overlap);

        for i in 0..segments.len() {
            let iperm = argsort(&rc_rank[i]);
            rc_rank[i] = permutation_apply(&iperm, &rc_rank[i]);
            seq_pos[i] = permutation_apply(&iperm, &seq_pos[i]);
        }

        Self{
            segment_join, sa, lcp, id, pos,
            path_join, seg_len: len, seq_pos, rc_rank,
            overlap 
        }
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

