#[test]
fn test_vs_bustools() {
    // comparing our implementation vs kallisto-bustools on a real bus file
    // run:
    // RUST_BACKTRACE=1 cargo test --release --package rustbustools --lib -- --nocapture count::test::test_vs_bustools
    use std::process::Command;

    let outfolder = "/tmp/bustest_rust";
    let outfolder_kallisto = "/tmp/bustest_kallisto";

    let busfolder = "/home/michi/bus_testing/bus_output_shorter";
    // let t2g = "/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt";
    let t2g = "/home/michi/bus_testing/transcripts_to_genes.txt";

    let bfolder = BusFolder::new(busfolder, t2g);
    let tfile = bfolder.get_transcript_file();
    let ecfile = bfolder.get_ecmatrix_file();
    let bfile = bfolder.get_busfile();

    fs::create_dir_all(outfolder).expect("Failed to create outfolder");
    fs::create_dir_all(outfolder_kallisto).expect("Failed to create outfolder_kallisto");

    // -------------------
    // Doing our own count
    // -------------------
    let now = Instant::now();
    let c = count(&bfolder, false);
    let elapsed_time = now.elapsed();
    println!("count::count in in {:?}", elapsed_time);
    c.write(outfolder);

    let now = Instant::now();
    let c2 = count2::count(&bfolder, false);
    let elapsed_time = now.elapsed();
    println!("count2::count in in {:?}", elapsed_time);
    assert!(c2.is_equal(&c));

    // -------------------
    // Bustools count
    // -------------------
    // Command
    // println!("source ~/miniconda3_newtry/bin/activate June2021");
    println!("bustools count -o {outfolder_kallisto}/gene -e {ecfile} -g {t2g} -t {tfile} --genecounts  {bfile}");

    let cmat_kallisto = CountMatrix::from_disk(
        &format!("{outfolder_kallisto}/gene.mtx"),
        &format!("{outfolder_kallisto}/gene.barcodes.txt"),
        &format!("{outfolder_kallisto}/gene.genes.txt"),
    );

    let cmat_rust = CountMatrix::from_disk(
        &format!("{outfolder}/gene.mtx"),
        &format!("{outfolder}/gene.barcodes.txt"),
        &format!("{outfolder}/gene.genes.txt"),
    );
    // println!("kallisto \n{:?}", cmat_kallisto.to_map());
    // println!("rust\n{:?}", cmat_rust.to_map());

    let sum1: usize = cmat_kallisto.matrix.iter().map(|(v, _s)| *v).sum();
    let sum2: usize = cmat_rust.matrix.iter().map(|(v, _s)| *v).sum();
    println!("{sum1} {sum2}");

    let h1 = cmat_kallisto.to_map();
    let h2 = cmat_rust.to_map();

    for (k, v) in h1.iter() {
        match h2.get(k) {
            Some(v2) => {
                if *v != *v2 {
                    println!("h2 wrong value {k:?}, {v}  <-> {k:?}, {v2}")
                }
            }
            None => {
                println!("h1 unique: {k:?}, {v}")
            }
        }
    }
    for (k, v) in h2.iter() {
        match h1.get(k) {
            Some(v2) => {
                if *v != *v2 {
                    println!("h1 wrong value {k:?}, {v}  <-> {k:?}, {v2}")
                }
            }
            None => {
                println!("h2 unique: {k:?}, {v}")
            }
        }
    }

    assert!(cmat_kallisto.is_equal(&cmat_rust))
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