// https://github.com/rust-bio/rust-bio/issues/3
// SAIS is generic, but suffix_array is not... WTF
// https://stackoverflow.com/a/41597323
use bio::data_structures::suffix_array::suffix_array;

pub trait SuffixArray {
    fn create(&self) -> Vec<usize>;
}

impl SuffixArray for [u8] {
    fn create(&self) -> Vec<usize> {
        return suffix_array(self);
    }
}

impl SuffixArray for [usize] {
    fn create(&self) -> Vec<usize> {
        let n = self.len();
        let mut v: Vec<(&[usize], usize)> = Vec::with_capacity(n);
        for i in 0..n {
            v.push((&self[i..], i));
        }
        v.sort_unstable();
        let result = v.iter().map(|x| x.1).collect();
        return result;
    }
}

#[cfg(test)]
mod tests {
    use crate::arrays::SuffixArray;

    #[test]
    fn test_suffix_array() {
        let v: Vec<u8> = b"GATCAG$".to_vec();
        let sa = SuffixArray::create(&v[..]);
        // 6: $GATCAG
        // 4: AG$GATC
        // 1: ATCAG$G
        // 3: CAG$GAT
        // 5: G$GATCA
        // 0: GATCAG$
        // 2: TCAG$GA
        assert_eq!(vec![6, 4, 1, 3, 5, 0, 2], sa);

        let v: Vec<usize> = vec![3, 1, 4, 2, 1, 3, 0];
        let sa = SuffixArray::create(&v[..]);
        assert_eq!(vec![6, 4, 1, 3, 5, 0, 2], sa);
    }
}


