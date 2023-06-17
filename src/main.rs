use clap::{self, Args, Parser, Subcommand};
use itertools::Itertools;
// use rand::distributions::Uniform;
use rustbustools::butterfly::make_ecs;
use rustbustools::consistent_genes::{GeneId, Genename, EC};
use rustbustools::correct;
use rustbustools::inspect::inspect;
use rustbustools::io::{BusHeader, BusReader};
use rustbustools::iterators::CellGroupIterator;
use rustbustools::sort::sort_on_disk;
use rustbustools::utils::int_to_seq;
use rustbustools::{busmerger, io::BusFolder};
// use std::collections::{HashMap, HashSet};
use std::fs::{self, File};
use std::io::{BufWriter, Write};
// use std::time::Instant;

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {
    /// Path to output file
    #[clap(short = 'o', long = "output")]
    output: String,

    #[clap(subcommand)]
    command: MyCommand,
}

#[allow(non_camel_case_types)]
#[derive(Subcommand)]
enum MyCommand {
    busmerge(BusMergeArgs),
    count(CountArgs),
    count2(CountArgs),
    resolve_ec(ResolveArgs),
    inspect(InspectArgs),
    sort(SortArgs),
    getcb(GetCBArgs),
    butterfly(ButterflyArgs),
    correct(CorrectArgs),
}

/// correct CBs with whitelist
#[derive(Args)]
struct CorrectArgs {


    /// Input busfile
    #[clap(long = "ifile", short = 'i')]
    inbus: String,

    /// Cell Barcode Whitelist
    #[clap(long = "whitelist")]
    whitelist: String,
}

/// Buttefly/ amplification profile
#[derive(Args)]
struct ButterflyArgs {
    /// input busfolder
    #[clap(long = "ifile", short = 'i')]
    inbus: String,
    /// Transcript-to-gene file
    #[clap(long = "t2g")]
    t2g: String,
    /// CB-UMI entries with multiple ECs will be collapsed into a single record (if they are consistent with a single gene)
    #[clap(long = "collapse")]
    collapse_ec: bool,
}

/// Sort busfile by CB/UMI/EC
#[derive(Args)]
struct SortArgs {
    /// input busfolder
    #[clap(long = "ifile", short = 'i')]
    inbus: String,
}

/// count the mRNAs  per cell and write to file
#[derive(Args)]
struct GetCBArgs {
    /// input busfolder
    #[clap(long = "ifile", short = 'i')]
    inbus: String,
}

/// countmatrix from busfile
#[derive(Args)]
struct CountArgs {
    /// input busfolder
    #[clap(long = "ifolder")]
    inbus: String,

    /// Transcript-to-gene file
    #[clap(long = "t2g")]
    t2g: String,

    /// ignore multimapped busrecords (same CB/UMI but different EC)
    #[clap(long = "ignoremm")]
    ignoremm: bool,
}

/// find overlap between busfiles and write out overlapping molecules
#[derive(Args)]
struct BusMergeArgs {
    /// 1st Input busfile
    #[clap(long = "i1")]
    inbus1: String,
    /// 2nd Input busfile
    #[clap(long = "i2")]
    inbus2: String,

    /// 1st output busfile
    #[clap(long = "o1")]
    outbus1: String,
    /// 2nd output busfile
    #[clap(long = "o2")]
    outbus2: String,
}

/// resovle an EC into gene names
#[derive(Args)]
struct ResolveArgs {
    /// input busfolder
    #[clap(long = "ifolder")]
    inbus: String,
    #[clap(long = "t2g")]
    /// Transcript-to-gene file
    t2g: String,

    /// Equivalence class to query genes for
    #[clap(long = "ec")]
    ec: u32,
}

/// Inspect busfile for stats
#[derive(Args)]
struct InspectArgs {
    /// input busfolder
    #[clap(short = 'i', long = "input")]
    inbus: String,
}

// mod count;
use rustbustools::count;
use rustbustools::count2;

fn main() {
    // use rustbustools::multinomial::{test_multinomial_stats, test_multinomial, multinomial_sample};

    // use std::time::Instant;

    // let dim = 10_000_000;

    // let n = 1_000_000;
    // let p: Vec<f64> = (1..dim).map(|x| x as f64).rev().collect();

    // let now = Instant::now();
    // multinomial_sample(n, p);
    // let elapsed_time = now.elapsed();
    // println!("Running multinomial_sample({}) took {} seconds.", dim, elapsed_time.as_secs());

    // let now = Instant::now();
    // test_multinomial_stats(dim);
    // let elapsed_time = now.elapsed();
    // println!("Running test_multinomial_stats({}) took {} seconds.", dim, elapsed_time.as_secs());

    // let now = Instant::now();
    // test_multinomial(dim);
    // let elapsed_time = now.elapsed();
    // println!("Running test_multinomial({}) took {} seconds.", dim, elapsed_time.as_secs());

    let cli = Cli::parse();
    match cli.command {
        MyCommand::busmerge(args) => {
            println!("Doing bus merging");
            busmerger::merge_busfiles_on_overlap(
                &args.inbus1,
                &args.inbus2,
                &args.outbus1,
                &args.outbus2,
            )
        }
        MyCommand::count(args) => {
            println!("Doing count");

            fs::create_dir(&cli.output).unwrap();

            let bfolder = BusFolder::new(&args.inbus, &args.t2g);
            let c = count::count(&bfolder, args.ignoremm);

            c.write(&cli.output);
            // write_matrix_market(&cli.output, &c.matrix).unwrap();
        }
        MyCommand::count2(args) => {
            println!("Doing count");
            fs::create_dir(&cli.output).unwrap();

            let bfolder = BusFolder::new(&args.inbus, &args.t2g);

            let c = count2::count(&bfolder, args.ignoremm);
            c.write(&cli.output);
            // write_matrix_market(&cli.output, &c.matrix).unwrap();
        }

        MyCommand::resolve_ec(args) => {
            println!("Doing resolve");
            let bfolder = BusFolder::new(&args.inbus, &args.t2g);
            let mut genes: Vec<&GeneId> =
                bfolder.ec2gene.get_genes(EC(args.ec)).into_iter().collect();
            genes.sort();
            println!("EC {} -> {:?}", args.ec, genes);

            let mut genenames: Vec<Genename> = bfolder
                .ec2gene
                .get_genenames(EC(args.ec))
                .into_iter()
                .collect();
            genenames.sort();

            println!("EC {} -> {:?}", args.ec, genenames);
        }
        MyCommand::inspect(args) => {
            inspect(&args.inbus);
        }

        MyCommand::getcb(args) => {
            let fh = File::create(cli.output).unwrap();
            let mut writer = BufWriter::new(fh);
            // let cb_len = 16;

            let header = BusHeader::from_file(&args.inbus);
            let cb_len = header.get_cb_len() as usize;
            let bus_cb = BusReader::new(&args.inbus)
                .groupby_cb()
                .map(|(cb, records)| {
                    (
                        // CB,decoded
                        int_to_seq(cb, cb_len),
                        // number of UMIs
                        records.iter().map(|r| r.UMI).unique().count(),
                    )
                });

            for (cb, nrecords) in bus_cb {
                writeln!(writer, "{},{}", cb, nrecords).unwrap();
            }
        }
        MyCommand::sort(args) => {
            let chunksize = 10_000_000; // roughly 300MB size
            sort_on_disk(&args.inbus, &cli.output, chunksize)
        }
        MyCommand::butterfly(args) => {
            let bfolder = BusFolder::new(&args.inbus, &args.t2g);
            let cuhist = make_ecs(&bfolder, args.collapse_ec);
            cuhist.to_disk(&cli.output);
        }
        MyCommand::correct(args) => {
            correct::correct(&args.inbus, &cli.output, &args.whitelist);
        }
    }
}

/*
flamegraph --flamechart  -- ~/rust_target/release/rustbustools --output /dev/null count --ifolder /home/michi/bus_testing/bus_output_shorter --t2g /home/michi/bus_testing/transcripts_to_genes.txt
 */
