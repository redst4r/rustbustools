use indicatif::{ProgressBar, ProgressStyle};


pub fn seq_to_int(seq: String) -> u64{

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
    let int_seq = u64::from_str_radix(&s, 4).unwrap();
    int_seq
}


pub fn int_to_seq(i: u64, seq_len:u64) -> String{

    let mut q = i;
    let mut result: Vec<u64> = Vec::with_capacity(seq_len as usize);
    while q>= 4{
        let quotient = q / 4;
        let remainder = q % 4;
        result.push(remainder);
        q = quotient;
    }
    result.push(q);

    while result.len() < seq_len as usize{
        result.push(0);
    }

    result.reverse();

    let s: String = result.iter()
    .map(|x| match x {
        0 => 'A',
        1 => 'C',
        2 => 'G',
        3 => 'T',
        _ => panic!("unkown seq character"),
    })
    .collect::<String>();

    // println!("{:?}", s);
    s
}   

pub fn get_progressbar(total: u64) -> ProgressBar{
    let bar = ProgressBar::new(total);
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise} ETA {eta}] {bar:40.cyan/blue} {pos}/{len} {per_sec}")
        .progress_chars("##-"));
    bar
}


mod tests {
    #[test]
    fn encode_seq(){
        use crate::utils::seq_to_int;
        assert_eq!(seq_to_int("A".to_string()), 0);
        assert_eq!(seq_to_int("C".to_string()), 1);
        assert_eq!(seq_to_int("G".to_string()), 2);
        assert_eq!(seq_to_int("T".to_string()), 3);
        assert_eq!(seq_to_int("GCCA".to_string()), 148);
    }

    #[test]
    fn decode_seq(){
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
