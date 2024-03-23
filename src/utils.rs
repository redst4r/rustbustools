//! Utilities
use indicatif::{ProgressBar, ProgressStyle};
use std::collections::HashSet;
use std::hash::Hash;

/// turning a vector into a HashSet
pub fn vec2set<T: Eq + Hash>(x: Vec<T>) -> HashSet<T> {
    x.into_iter().collect::<HashSet<T>>()
}

/// Encoding a base sequence into int
pub fn seq_to_int(seq: String) -> u64 {
    assert!(seq.len() <= 32); // cant handle longer sequences in a single 64bit integer!
    let s: String = seq
        .chars()
        .map(|x| match x {
            'A' => '0',
            'C' => '1',
            'G' => '2',
            'T' => '3',
            _ => panic!("unkown seq character"),
        })
        .collect();
    u64::from_str_radix(&s, 4).unwrap()
}

/// Decoding an int into a base sequence. len of the sequence must be specified
pub fn int_to_seq(i: u64, seq_len: usize) -> String {
    let mut q = i;
    let mut result: Vec<u64> = Vec::with_capacity(seq_len);
    while q >= 4 {
        let quotient = q / 4;
        let remainder = q % 4;
        result.push(remainder);
        q = quotient;
    }
    result.push(q);

    while result.len() < seq_len {
        result.push(0);
    }
    result.reverse();

    let s: String = result
        .iter()
        .map(|x| match x {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => panic!("unkown seq character"),
        })
        // .collect::<String>();
        .collect();
    s
}

/// returns a progress bar instance with standard formatting
pub fn get_progressbar(total: u64) -> ProgressBar {
    let bar = ProgressBar::new(total);
    bar.set_style(
        ProgressStyle::default_bar()
            .template("[{elapsed_precise} ETA {eta}] {bar:40.cyan/blue} {pos}/{len} {per_sec}")
            .unwrap()
            .progress_chars("##-"),
    );
    bar
}

pub mod argsort {
    //! A workaround for the unsortable `Vec<f64>` (due to Nan)
    //! # Example
    //! allows something like
    //! ```rust
    //! # use bustools::utils::argsort::{argsort_float,argmax_float};
    //! argsort_float(&vec![1.1_f64, -0.1_f64], true);
    //! argmax_float(&vec![1.0_f64, 10_f64]);
    //! ```
    //!
    //! # References  
    //! <https://stackoverflow.com/questions/69764050/how-to-get-the-indices-that-would-sort-a-vector-in-rust>
    //! <https://stackoverflow.com/questions/28247990/how-to-do-a-binary-search-on-a-vec-of-floats>

    use std::cmp::Ordering;

    #[derive(PartialEq, PartialOrd)]
    struct NonNan(f64);

    impl NonNan {
        fn new(val: f64) -> Option<NonNan> {
            if val.is_nan() {
                None
            } else {
                Some(NonNan(val))
            }
        }
    }
    impl Eq for NonNan {}

    impl Ord for NonNan {
        fn cmp(&self, other: &NonNan) -> Ordering {
            self.partial_cmp(other).unwrap()
        }
    }

    /// Argsort a `slice[T]`
    pub fn argsort<T: Ord>(slice: &[T]) -> Vec<usize> {
        let n = slice.len();
        let mut keys: Vec<_> = (0..n).collect();
        keys.sort_by_key(|x| &slice[*x]);
        keys
    }

    /// argsort of a f64 vector assuming no NAN (will panic otherwise)
    pub fn argsort_float(fvec: &Vec<f64>, ascending: bool) -> Vec<usize> {
        let _fvec: Vec<NonNan> = fvec
            .iter()
            .map(|x| NonNan::new(*x).unwrap_or_else(|| panic!("Nan values in {fvec:?}")))
            .collect();

        let mut fvec_sorted_ix = argsort(&_fvec);
        if ascending {
            fvec_sorted_ix.reverse();
        }
        fvec_sorted_ix
    }

    /// argmax of a f64 vector assuming no NAN (will panic otherwise)
    pub fn argmax_float(fvec: &Vec<f64>) -> (usize, f64) {
        let ix = argsort_float(fvec, true);
        let i = ix[0];
        let value = fvec[i];
        (i, value)
    }

    /// argmin of a f64 vector assuming no NAN (will panic otherwise)
    pub fn argmin_float(fvec: &Vec<f64>) -> (usize, f64) {
        let ix = argsort_float(fvec, false);
        let i = ix[0];
        let value = fvec[i];
        (i, value)
    }

    #[test]
    fn test_argmin() {
        let v = vec![1.0, -1.0, 10.0, 0.0];
        assert_eq!(argmin_float(&v), (1, -1.0));
    }

    #[test]
    fn test_argmax() {
        let v = vec![1.0, -1.0, 10.0, 0.0];
        assert_eq!(argmax_float(&v), (2, 10.0));
    }
}

#[cfg(test)]
mod tests {
    #[test]
    fn encode_seq() {
        use crate::utils::seq_to_int;
        assert_eq!(seq_to_int("A".to_string()), 0);
        assert_eq!(seq_to_int("C".to_string()), 1);
        assert_eq!(seq_to_int("G".to_string()), 2);
        assert_eq!(seq_to_int("T".to_string()), 3);
        assert_eq!(seq_to_int("GCCA".to_string()), 148);
    }

    #[test]
    fn decode_seq() {
        use crate::utils::int_to_seq;

        //  base order
        assert_eq!(int_to_seq(0, 1), "A");
        assert_eq!(int_to_seq(1, 1), "C");
        assert_eq!(int_to_seq(2, 1), "G");
        assert_eq!(int_to_seq(3, 1), "T");
        //
        // # padding leading A's
        assert_eq!(int_to_seq(0, 3), "AAA");
        assert_eq!(int_to_seq(1, 3), "AAC");
        assert_eq!(int_to_seq(2, 3), "AAG");
        assert_eq!(int_to_seq(3, 3), "AAT");
        //
        assert_eq!(int_to_seq(4, 2), "CA");
        assert_eq!(int_to_seq(5, 2), "CC");
        assert_eq!(int_to_seq(6, 2), "CG");
        assert_eq!(int_to_seq(7, 2), "CT");
        //
        assert_eq!(int_to_seq(148, 4), "GCCA");
        //
        // # make sure to raise an error when the decoded string is actually longer
        // # then requested (since its probably a bug in the code calling _decode_int_to_ACGT)
        // with pytest.raises(AssertionError):
        // busio._decode_int_to_ACGT(148, seq_len=2)
        // with pytest.raises(AssertionError):
        // busio._decode_int_to_ACGT(-1, seq_len=1)
    }
}
