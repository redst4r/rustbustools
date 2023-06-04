use std::{fs, time::Instant, io::{self, Write}};

use rustbustools::{count::{count}, io::{BusFolder, BusReader, write_partial_busfile}, count2, iterators::CellGroupIterator};
use rustbustools::countmatrix::CountMatrix;

// pub const TEST_T2G: &str = "/home/michi/transcripts_to_genes.txt";
// pub const TEST_BUSFILE: &str = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/DSP1/kallisto/sort_bus/bus_output/output.corrected.sort.bus";
// pub const TEST_BUSFOLDER: &str = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/DSP1/kallisto/sort_bus/bus_output/";

pub const TEST_T2G: &str = "/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt";
// pub const TEST_BUSFILE: &str = "/home/michi/ISB_data/LT_pilot/LT_pilot/kallisto_quant/DSP1/kallisto/sort_bus/bus_output/output.corrected.sort.bus";
// pub const TEST_BUSFOLDER: &str = "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/DSP1/kallisto/sort_bus/bus_output/";


pub const TEST_BUSFILE: &str = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
pub const TEST_BUSFOLDER: &str = "/home/michi/bus_testing/bus_output/";



#[test]
fn test_vs_bustools() {
    // comparing our implementation vs kallisto-bustools on a real bus file
    // run:
    // RUST_BACKTRACE=1 cargo test --release --package rustbustools --lib -- --nocapture count::test::test_vs_bustools
    use std::process::Command;

    let outfolder = "/tmp/bustest_rust";
    let outfolder_kallisto = "/tmp/bustest_kallisto";


    let bfolder = BusFolder::new(TEST_BUSFOLDER, TEST_T2G);
    let tfile = bfolder.get_transcript_file();
    let ecfile = bfolder.get_ecmatrix_file();
    let bfile = bfolder.get_busfile();

    fs::create_dir_all(outfolder).expect("Failed to create outfolder");
    fs::create_dir_all(outfolder_kallisto).expect("Failed to create outfolder_kallisto");

    // -------------------
    // Doing our own count
    // -------------------
    println!("Doing count::count");
    let now = Instant::now();
    let c = count(&bfolder, false);
    let elapsed_time = now.elapsed();
    println!("count::count in in {:?}", elapsed_time);
    c.write(outfolder);

    println!("Doing count::count2");
    let now = Instant::now();
    let c2 = count2::count(&bfolder, false);
    let elapsed_time = now.elapsed();
    println!("count2::count in in {:?}", elapsed_time);
    assert!(c2.is_equal(&c));

    // -------------------
    // Bustools count
    // -------------------
    // Command

    let bustools_binary = "/home/michi/miniconda3_newtry/envs/nextflow_bioinformatics/bin/bustools";
    println!("Doing kallisto::count");
    let now = Instant::now();
    let output = Command::new(bustools_binary)
        .arg("count")
        .arg("-o")
        .arg(format!("{outfolder_kallisto}/gene"))
        .arg("-e")
        .arg(ecfile)
        .arg("-g")
        .arg(TEST_T2G)
        .arg("-t")
        .arg(tfile)
        .arg("--genecounts")
        .arg(bfile)
        .output().unwrap();
    let elapsed_time = now.elapsed();
    println!("kallisto::count in in {:?}", elapsed_time);

    println!("status: {}", output.status);

    // -------------------
    // Comparing results
    // -------------------
    let cmat_kallisto = CountMatrix::from_disk(
        &format!("{outfolder_kallisto}/gene.mtx"),
        &format!("{outfolder_kallisto}/gene.barcodes.txt"),
        &format!("{outfolder_kallisto}/gene.genes.txt"),
    );

    // no need to load from disk!!
    // let cmat_rust = CountMatrix::from_disk(
    //     &format!("{outfolder}/gene.mtx"),
    //     &format!("{outfolder}/gene.barcodes.txt"),
    //     &format!("{outfolder}/gene.genes.txt"),
    // );
    let cmat_rust = c;


    let sum1: usize = cmat_kallisto.matrix.iter().map(|(v, _s)| *v).sum();
    let sum2: usize = cmat_rust.matrix.iter().map(|(v, _s)| *v).sum();
    assert_eq!(sum1, sum2);

    // let h1 = cmat_kallisto.to_map();
    // let h2 = cmat_rust.to_map();
    // for (k, v) in h1.iter() {
    //     match h2.get(k) {
    //         Some(v2) => {
    //             if *v != *v2 {
    //                 println!("h2 wrong value {k:?}, {v}  <-> {k:?}, {v2}")
    //             }
    //         }
    //         None => {
    //             println!("h1 unique: {k:?}, {v}")
    //         }
    //     }
    // }
    // for (k, v) in h2.iter() {
    //     match h1.get(k) {
    //         Some(v2) => {
    //             if *v != *v2 {
    //                 println!("h1 wrong value {k:?}, {v}  <-> {k:?}, {v2}")
    //             }
    //         }
    //         None => {
    //             println!("h2 unique: {k:?}, {v}")
    //         }
    //     }
    // }

    assert!(cmat_kallisto.is_equal(&cmat_rust));
}

    // #[test]
    fn test_cb_iter_speed() {
        use std::time::Instant;
        let foldername = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
        let n = 100000;

        let b = BusReader::new(foldername);
        let biter2 = b.groupby_cb();

        let now = Instant::now();
        let _s2: Vec<_> = biter2.take(n).map(|(_a, records)| records).collect();
        let elapsed_time = now.elapsed();
        println!(
            "Running CellIterator({}) took {} seconds.",
            n,
            elapsed_time.as_secs()
        );
    }


    // #[test]
    fn test_write() {
        // let fname = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
        let fname = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
        // let outname = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        // write_partial_busfile(fname, outname, 10_000_000);
    
        let outname = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";
        write_partial_busfile(fname, outname, 1_500_000);
    }