pub struct Block<'a> {
         id: &'a[usize],
        pos: &'a[usize],
    seq_pos: &'a[Vec<usize>],
    rc_rank: &'a[Vec<usize>],
        idx: Vec<usize>,
        rcr: Vec<usize>
}

impl<'a> Block<'a> {
    pub fn new(
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

