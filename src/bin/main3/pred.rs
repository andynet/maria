pub trait Predecessor {
    fn argpred(&self, item: usize) -> usize;
}

impl Predecessor for Vec<usize> {
    fn argpred(&self, item: usize) -> usize {
        let mut i = self.len() - 1;
        while self[i] > item { i -= 1; }
        return i;
    }
}

