use rustbustools::busmerger;

use clap::{self, Parser, Subcommand, Args};

#[derive(Parser)]
#[clap(author, version, about, long_about = None)]
struct Cli {  
    /// Path to output file
    #[clap(short ='o', long = "output")] 
    output: String,    

    #[clap(subcommand)]
    command: MyCommand
}

#[allow(non_camel_case_types)]
#[derive(Subcommand)]
enum MyCommand {
    busmerge(BusMergeArgs),
}

#[derive(Args)]
struct BusMergeArgs{
    /// 1st Input busfile
    #[clap(long= "i1")] 
    inbus1: String,
    /// 2nd Input busfile
    #[clap(long= "i2")] 
    inbus2: String,

    /// 1st output busfile
    #[clap(long= "o1")] 
    outbus1: String,
    /// 2nd output busfile
    #[clap(long= "o2")] 
    outbus2: String,  
}


// mod count;
use rustbustools::count::{test_multinomial_stats, test_multinomial, multinomial_sample};
use std::time::Instant;
fn main() {

    let dim = 10_000_000;

    let n = 1_000_000;
    let p: Vec<f64> = (1..dim).map(|x| x as f64).rev().collect();

    let now = Instant::now();
    multinomial_sample(n, p);
    let elapsed_time = now.elapsed();
    println!("Running multinomial_sample({}) took {} seconds.", dim, elapsed_time.as_secs());

    // let now = Instant::now();
    // test_multinomial_stats(dim);
    // let elapsed_time = now.elapsed();
    // println!("Running test_multinomial_stats({}) took {} seconds.", dim, elapsed_time.as_secs());


    // let now = Instant::now();
    // test_multinomial(dim);
    // let elapsed_time = now.elapsed();
    // println!("Running test_multinomial({}) took {} seconds.", dim, elapsed_time.as_secs());



    // let cli = Cli::parse();
    // match cli.command{
    //     MyCommand::busmerge(args) => {
    //         println!("Doing bus merging");
    //         busmerger::merge_busfiles_on_overlap(args.inbus1, args.inbus2, args.outbus1, args.outbus2)      
    //     }
    // }


}


