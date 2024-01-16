use std::cmp::Ordering;

pub trait Predecessor {
    fn argpred(&self, item: usize) -> usize;
}

impl Predecessor for Vec<usize> {
    /// Returns index of an element smaller or equal to item
    /// Assumes self is monotonically increasing and item is larger than at
    /// least one element of self
    fn argpred(&self, item: usize) -> usize {
        let mut l = 0;
        let mut r = self.len();

        while l < r-1 {
            let m = (l + r) / 2;
            match item.cmp(&self[m]) {
                Ordering::Greater => { l = m; },
                Ordering::Less    => { r = m; },
                Ordering::Equal   => { return m; },
            }
        }
        return l;
    }
}

#[test]
fn test_predecessor() {
    let v = vec![0, 3, 5, 9];

    assert!(v.argpred(4) == 1);
    assert!(v.argpred(5) == 2);
    assert!(v.argpred(6) == 2);
    assert!(v.argpred(10) == 3);
}

