//! Code for `bustools correct` to correct sequencing errors in Cell Barcodes using a whitelist
//! 
//! Pretty straight forward: Operates via a `BKTree`, which allows for quick 
//! "approximate" matching
//! 
#![deny(missing_docs)]
use std::{
    collections::{HashMap, HashSet},
    fs::File,
    io::{BufRead, BufReader},
};

use bktree::BkTree;

use crate::{
    io::{BusReader, BusWriter},
    utils::{get_progressbar, int_to_seq, seq_to_int},
};

const MAX_DIST: isize = 1;

fn my_hamming(a: &String, b: &String) -> isize {
    // hamming distance for two strings of the same size
    assert_eq!(a.len(), b.len());
    let mut counter: isize = 0;
    // for (c1, c2) in  std::iter::zip((*a).chars(), (*b).chars()){  // todo: change to bytes, might be faster
    for (c1, c2) in std::iter::zip((*a).bytes(), (*b).bytes()) {
        // todo: change to bytes, might be faster
        if c1 != c2 {
            counter += 1;
        }
    }
    counter
}

#[derive(Debug, Eq, PartialEq)]
enum CorrectionResult {
    SingleHit(String),
    NoHit,
    Ambigous(Vec<String>),
}

/// Correct a single barcode using the whitelist (represented as a BKTree)
/// Checks if any whitelisted barcode is <= 1 away from the query
fn correct_single_cb(cb: String, bk: &BkTree<String>) -> CorrectionResult {
    let matches = bk.find(cb, MAX_DIST);
    match matches.len() {
        0 => CorrectionResult::NoHit,
        1 => {
            let (new_cb, _distance) = matches[0];
            CorrectionResult::SingleHit(new_cb.to_owned())
        }
        _ => {
            // more complicated there
            // bktree find also returns EXACT matches!
            let perfect_match: Vec<String> = matches
                .iter()
                .filter_map(|(cb, dist)| {
                    if *dist == 0 {
                        Some((*cb).clone())
                    } else {
                        None
                    }
                })
                .collect();
            if perfect_match.len() == 1 {
                let cb_correct = perfect_match.get(0).unwrap().clone();
                CorrectionResult::SingleHit(cb_correct)
            } else {
                // panic!("Shouldnt happen. Whitelist shouldnt have two hits 2BP appart: {:?}", matches),
                // actually it does happen: the query can fall exactly between two whitelisted CBs
                // just remove it
                let multi: Vec<String> = matches
                    .into_iter()
                    .map(|(cb_whitelist, _dist)| cb_whitelist.clone())
                    .collect();
                CorrectionResult::Ambigous(multi)
            }
        }
    }
}

/// Corrects observed barcodes in the busfile using a whitelist of barcodes and writes the results to disk
/// 
/// # Parameters
/// * `busfile`: filename of the busfile to be corrected
/// * `busfile_out`: file where the corrected records are written
/// * `whitelist_filename` : the file with the whitelisted barcodes (one per line)
/// 
/// # Overview/Performance tricks
/// The CBs are highly repetitive; would be slow to query the BKtree for each CB (they'll repeat ALOt)
/// 1. gather all the unique CBs in the busfile
/// 2. correct them and create a HashMap<uncorrected, corrected>
/// 3. iterate over the bus file, correct the individual entries and write to disk
/// 
pub fn correct(busfile: &str, busfile_out: &str, whitelist_filename: &str) {
    println!("Loading whitelist");
    let whitelist = load_whitelist(whitelist_filename);
    println!("Loaded whitelist");

    println!("Building BKTree");
    let mut bk: BkTree<String> = BkTree::new(my_hamming);
    bk.insert_all(whitelist.clone().into_iter());
    println!("Built BKTree");

    let breader = BusReader::new(busfile);
    let cb_len = breader.bus_header.cb_len as usize;

    // note the file might be unsorted, so cant realy on groupby_cb
    println!("collecting CBs");
    let unique_cbs: HashSet<String> = breader.map(|r| int_to_seq(r.CB, cb_len)).collect();
    println!("collected CBs");

    println!("correcting unique CBs");
    // mapping on the int represnetation of the barcodes! saves some time
    let mut corrector: HashMap<u64, u64> = HashMap::with_capacity(unique_cbs.len());
    let bar = get_progressbar(unique_cbs.len() as u64);
    let mut cb_correct = 0;
    let mut cb_total = 0;
    for (counter, cb) in unique_cbs.into_iter().enumerate() {
        cb_total += 1;

        // to save time (BKtree is slow) check if we have a direct match
        if whitelist.contains(&cb) {
            let cbint = seq_to_int(cb);
            corrector.insert(cbint, cbint);
            cb_correct += 1
        // if its not a direct match, check the BKTree for 1 error
        } else if let CorrectionResult::SingleHit(corrected_cb) = correct_single_cb(cb.clone(), &bk)
        {
            corrector.insert(seq_to_int(cb), seq_to_int(corrected_cb));
            cb_correct += 1
        } else {
            // simply dont do anything. Later if we look up a query-CB and cant find it in the map
            // it cant be corrected!
        }

        if counter % 10_000 == 0 {
            bar.inc(10_000)
        }
    }
    println!("corrected unique CBs: {cb_correct}/{cb_total}");

    println!("writing corrected busfile");

    // now with a map of uncorrected->corrected fix the busfile
    let breader_tmp = BusReader::new(busfile);
    let total_records = breader_tmp.count();
    let bar = get_progressbar(total_records as u64);

    let breader = BusReader::new(busfile);
    let mut bwriter = BusWriter::new(busfile_out, breader.bus_header.clone());

    for (counter, record) in breader.enumerate() {
        if let Some(corrected_cb) = corrector.get(&record.CB) {
            let mut new_record = record.clone();
            new_record.CB = *corrected_cb;
            bwriter.write_record(&new_record);
        }

        if counter % 10_000 == 0 {
            bar.inc(10_000)
        }
    }
    println!("wrote corrected busfile");
}

/// Parse the whitelist-file (one whitelisted barcode per line) into a HashSet
fn load_whitelist(whitelist_filename: &str) -> HashSet<String> {
    let whitelist_reader = BufReader::new(File::open(whitelist_filename).unwrap());
    let whitelist_header: HashSet<String> = whitelist_reader.lines().map(|f| f.unwrap()).collect();
    whitelist_header
}

#[cfg(test)]
mod testing {
    use bktree::BkTree;

    use crate::correct::{correct, correct_single_cb, load_whitelist, CorrectionResult};

    use super::my_hamming;
    #[test]
    fn test_correct() {
        let whitelist = vec!["AAAA".to_string(), "BBBB".to_string()];
        let mut bk: BkTree<String> = BkTree::new(my_hamming);
        bk.insert_all(whitelist.into_iter());

        // perfect match
        assert_eq!(
            correct_single_cb("AAAA".to_string(), &bk),
            CorrectionResult::SingleHit("AAAA".to_string())
        );

        // one mismatch match
        assert_eq!(
            correct_single_cb("AAAB".to_string(), &bk),
            CorrectionResult::SingleHit("AAAA".to_string())
        );

        // too far away match
        assert_eq!(
            correct_single_cb("BBAA".to_string(), &bk),
            CorrectionResult::NoHit
        );

        let whitelist = vec!["AAAA".to_string(), "AABB".to_string()];
        let mut bk: BkTree<String> = BkTree::new(my_hamming);
        bk.insert_all(whitelist.into_iter());
        // two hits, not clear which one
        assert_eq!(
            correct_single_cb("AABA".to_string(), &bk),
            CorrectionResult::Ambigous(vec!["AAAA".to_string(), "AABB".to_string()])
        );

        // make sure that a perfect match is respected too
        let whitelist = vec!["AAAA".to_string(), "AAAB".to_string()];
        let mut bk: BkTree<String> = BkTree::new(my_hamming);
        bk.insert_all(whitelist.into_iter());
        // two hits, not clear which one
        assert_eq!(
            correct_single_cb("AAAA".to_string(), &bk),
            CorrectionResult::SingleHit("AAAA".to_string())
        );
    }

    #[test]
    fn test_real_file() {
        // pub const TEST_BUSFILE: &str = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
        // pub const TEST_BUSFOLDER: &str = "/home/michi/bus_testing/bus_output/";
        pub const TEST_BUSFILE: &str =
            "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        // pub const TEST_BUSFOLDER: &str = "/home/michi/bus_testing/bus_output_short/";
        // let the_whitelist = "/home/michi/mounts/TB4drive/kallisto_resources/3M-february-2018.txt";
        let the_whitelist = "/home/michi/bus_testing/3M-february-2018.txt";
        correct(TEST_BUSFILE, "/tmp/corrected.bus", the_whitelist)
    }

    #[test]
    fn test_bk() {
        let the_whitelist = "/home/michi/bus_testing/3M-february-2018.txt";
        let wl = load_whitelist(the_whitelist);
        let mut bk: BkTree<String> = BkTree::new(my_hamming);
        bk.insert_all(wl.into_iter());

        let r = bk.find("TAACCTGAGACTCGGA".to_string(), 1);
        println!("{:?}", r);
    }
}


/*
cargo run --release -- --output /tmp/corr.bus correct\
  --ifile ~/bus_testing/bus_output_shorter/output.corrected.sort.bus \
  --whitelist ~/bus_testing/3M-february-2018.txt
 */