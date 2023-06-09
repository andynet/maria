use bio::data_structures::suffix_array::lcp as lcp_array;
use gfa::gfa::Path;
use gfa::gfa::Segment;
use gfa::parser::GFAParser;
use maria::arrays::SuffixArray;
use std::ops::Add;
use std::str;

/// Data required to iterate through suffix array in sublinear space
pub struct PFData {
    /// suffix array of segment join
     sa: Vec<usize>, 
    /// longest common prefix of segment join
    lcp: Vec<isize>,
    /// id of segments associated with suffix array
     id: Vec<usize>,
    /// position in the segment associated with suffix array
    pos: Vec<usize>,

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
    i: usize,
    j: usize,
    idx: Option<Vec<usize>>,
    rcr: Vec<usize>,
    block: Option<Block<'a>>,
}

impl PFData {
    pub fn from_pfgraph(filename: &str) -> Self {
        let parser: GFAParser<usize, ()> = GFAParser::new();
        let gfa = parser.parse_file(filename)
            .expect("Error parsing GFA file.");

        let overlap = 1;
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

        Self{ sa, lcp, id, pos, seg_len: len, seq_pos, rc_rank, overlap } 
    }

    pub fn iter(&self) -> PFDataIterator {
        return PFDataIterator { 
            data: self,
            i: 0,
            j: self.seg_len.len() + 1, // #separators + sentinel
            idx: None,
            rcr: Vec::new(),
            block: None,
        }
    }
}

impl<'a> Iterator for PFDataIterator<'a> {
    type Item = (usize, usize, usize);

    // fn next(&mut self) -> Option<Self::Item> {
    //     if let Some(x) = self.block.next() { return x; }
    //     if let Some(self.block) = get_next_block() { return self.block.next(); }
    //     return None;
    // }
    fn next(&mut self) -> Option<Self::Item> {
        let data = self.data;

        if self.idx == None {
            // find new block
            //      find new start
            self.i = self.j;
            if self.i >= data.sa.len() { return None; }
            let mut remaining = data.seg_len[data.id[self.i]] - data.pos[self.i];
            while remaining <= data.overlap { 
                self.i += 1; 
                remaining = data.seg_len[data.id[self.i]] - data.pos[self.i];
            }
            //      find new end
            self.j = self.i+1;
            while data.lcp[self.j] >= remaining as isize {
                self.j += 1;
            }
            //      construct new block
            println!("block: {}..{}", self.i, self.j);
            self.block = Some(
                Block::new(
                    &data.id[self.i..self.j], 
                    &data.pos[self.i..self.j], 
                    &data.seq_pos, 
                    &data.rc_rank)
            );
            // self.idx = Some(vec![0; self.j - self.i]);
            // self.rcr = Vec::new();
            // for k in self.i..self.j { 
            //     self.rcr.push(data.rc_rank[data.id[k]][0]) 
            // }
        }

        let block = self.block.unwrap();
        let result = block.next();
        return result;
        // let k = argmin(&self.rcr);
        // // if self.rcr[k] == usize::MAX { return None; }
        // let Some(idx) = self.idx;
        // let result = (
        //     data.seq_pos[data.id[k]][idx[k]] + data.pos[k],
        //     data.id[k],
        //     data.pos[k]
        // );
        //
        // // self.idx[k] += 1;
        // let val = data.rc_rank[data.id[k]].get(idx[k]);
        // match val {
        //     None => { self.rcr[k] = usize::MAX; },
        //     Some(v) => { self.rcr[k] = *v; }
        // }
        //
        // return Some(result);

        // if self.processing == None {
        //     let remaining = data.seg_len[data.id[i]] - data.pos[i];
        //     if remaining <= data.overlap { i += 1; continue; }
        //     let mut j = i+1;
        //     while data.lcp[j] >= remaining as isize { j += 1; }

        //     let block = Block::new(&data.id[i..j], &data.pos[i..j], &data.seq_pos, &data.rc_rank);

        // }

        // return block.next(); 
        // i = j;
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

struct Block<'a> {
         id: &'a[usize],
        pos: &'a[usize],
    seq_pos: &'a[Vec<usize>],
    rc_rank: &'a[Vec<usize>],
        idx: Vec<usize>,
        rcr: Vec<usize>
}

impl<'a> Block<'a> {
    fn new(
        id: &'a[usize], pos: &'a[usize], 
        seq_pos: &'a[Vec<usize>], rc_rank: &'a[Vec<usize>]
    ) -> Block<'a> {

        let mut rcr = Vec::new();
        for i in id {
            rcr.push(rc_rank[*i][0]);
        }
        Block { 
            id, pos, seq_pos, rc_rank, 
            idx: vec![0; id.len()],
            rcr
        }
    }
}

impl<'a> Iterator for Block<'a> {
    type Item = (usize, usize, usize);

    fn next(&mut self) -> Option<Self::Item> {
        let k = argmin(&self.rcr);
        if self.rcr[k] == usize::MAX { return None; }
        let result = (
            self.seq_pos[self.id[k]][self.idx[k]] + self.pos[k],
            self.id[k],
            self.pos[k]
        );

        self.idx[k] += 1;
        let val = self.rc_rank[self.id[k]].get(self.idx[k]);
        match val {
            None => { self.rcr[k] = usize::MAX; },
            Some(v) => { self.rcr[k] = *v; }
        }

        return Some(result);
    }
}

fn argmin(data: &[usize]) -> usize {
    let mut min_pos = 0;
    let mut min_val = usize::MAX;
    for (p, &v) in data.iter().enumerate() {
        if v < min_val { min_pos = p; min_val = v; }
    }
    return min_pos;
}

