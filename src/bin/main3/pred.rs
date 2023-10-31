pub trait Predecessor {
    fn argpred(&self, item: usize) -> usize;
}

impl Predecessor for Vec<usize> {
    fn argpred(&self, item: usize) -> usize {
        let mut l = 0;
        let mut r = self.len();

        while l < r-1 {
            let m = (l + r) / 2;
            if item > self[m] { l = m; }
            else { r = m; }
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

