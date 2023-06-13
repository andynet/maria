pub mod arrays;
pub mod naive;
pub mod pf;

fn reverse_complement(s: &[u8]) -> Vec<u8> {
    let n = s.len();
    let mut result = vec![0; n];
    for (i, c) in s.iter().enumerate() {
        match *c {
            b'A' => { result[n-1-i] = b'T'; },
            b'C' => { result[n-1-i] = b'G'; },
            b'G' => { result[n-1-i] = b'C'; },
            b'T' => { result[n-1-i] = b'A'; },
            b'N' => { result[n-1-i] = b'N'; },
            _ => { panic!("Unexpected letter in s."); }
        }
    }
    return result;
}

#[cfg(test)]
mod tests;
