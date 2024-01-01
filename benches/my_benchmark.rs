use std::collections::{HashMap, HashSet};
use std::fs::File;

use bitvec::field::BitField;
use bustools::busz::{BuszWriter, BuszReader};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use itertools::Itertools;
use bustools::bus_multi::CellIteratorMulti;
use bustools::consistent_genes::{Ec2GeneMapper, Genename, EC};
use bustools::io::{BusReader, BusRecord, BusWriter, BusParams};
use bustools::iterators::{CellGroup, CbUmiGroup, CellGroupIterator, CbUmiGroupIterator};
use bustools::merger::MultiIterator;
use bustools::multinomial::{multinomial_sample, multinomial_sample_binary_search};

// use bustools::iterators_clone::CellIteratorClone;


/*
fn criterion_benchmark_multinomial(c: &mut Criterion) {
    // use statrs::distribution::Binomial as Binomial_statrs;
    // use probability::distribution::Binomial as Binomial_prob;
    use bustools::multinomial::{ multinomial_sample_statrs};
    // c.bench_function("fib 20", |b| b.iter(|| fibonacci(black_box(20))));

    // use bustools::io::parse_ecmatrix;

    // let fname = "/tmp/matrix.ec";
    // c.bench_function("ec1", |b| b.iter(|| parse_ecmatrix(black_box(fname.to_string()))));
    // c.bench_function("ec2", |b| b.iter(|| parse_ecmatrix2(black_box(fname.to_string()))));


    // multinomial benchmarks
    // c.bench_function("MN1", |b| b.iter(|| test_multinomial(black_box( 1000))));
    // c.bench_function("Mn2", |b| b.iter(|| test_multinomial_stats(black_box(1000))));

    // multinomial benchmarks, my implementation
    let n = 100000;
    let dim = 10000;
    let p: Vec<f64> = (1..dim).map(|x| x as f64).collect();

    c.bench_function("MN1", |b| b.iter(|| multinomial_sample(black_box( n),  black_box( p.clone()))));
    c.bench_function("MN2", |b| b.iter(|| multinomial_sample_statrs(black_box( n),  black_box( p.clone()))));
    // c.bench_function("Mn2", |b| b.iter(|| test_multinomial_stats(black_box(1000))));

    /*
    comparing the Binomial sampling algorthims in 
    statrs and probability crates
    1. statrs uses sum of bernoulli, which is very slow for large N
    2. probability uses inverse cdf sampling, much faster
    */
    // use probability::prelude::*;
    // use rand::distributions::Distribution;
    // let N = 10000;
    // let p = 0.001;

    // let mut source = source::default(42);
    // let mut r = rand::thread_rng();

    // let b1 = Binomial_prob::new(N, p);
    // let b2 = Binomial_statrs::new(p, N as u64).unwrap();

    // c.bench_function("Binomial prob", |b| b.iter(|| b1.sample(&mut source)));
    // c.bench_function("Binomial statrs", |b| b.iter(|| b2.sample(&mut r)));

}
 */
#[allow(dead_code)]
fn multinomial_speed(c: &mut Criterion){

    use probability::prelude::*;

    fn binary_search_dummy(N: u64, d: u64){

        let p = (1..d).map(|x| x as f64).collect();

        let mut random_source = source::default(4);   
        multinomial_sample_binary_search(N, &p, &mut random_source);
    }

    fn binomial_dummy(N: u64, d: u64){
        let p = (1..d).map(|x| x as f64).collect();

        let mut random_source = source::default(4);   
        multinomial_sample(N, &p, &mut random_source);
    }

    let dims = vec![10_000, 100_000, 1_000_000];
    let N = 1_000_000;
    for d in dims{

        let name = format!("Binary, dim {}", d);
        c.bench_function(&name, |b| b.iter(|| 
            binary_search_dummy(black_box(N), 
                            black_box(d), 
            )));

        let name = format!("Binomial, dim {}", d);
        c.bench_function(&name, |b| b.iter(|| 
            binomial_dummy(black_box(N), 
                            black_box(d), 
            )));
    }

}

#[allow(dead_code)]
fn plain_iterator_speed(c: &mut Criterion){
    
    fn _dummy(n: usize) -> u32{
        let busname = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
        let biter= BusReader::new(busname);
        let total = biter.take(n).map(|r|r.COUNT).sum::<u32>();
        total
    }

    c.bench_function("plain_iterator", |b| b.iter(|| _dummy(black_box(100000))));
}


/* testing the speed of my iterators (CB/UMI) against the grouped version
 */
#[allow(dead_code)]
fn iterator_speed(c: &mut Criterion){

     pub const FOLDERNAME:&str = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";

    fn dummy_CB_group(n: usize) ->Vec<Vec<BusRecord>>{
        let biter= BusReader::new(FOLDERNAME);
        let s:Vec<Vec<BusRecord>> = biter
            .group_by(|r| r.CB)
            .into_iter().map(|(_a, records)|records.collect())
            .take(n)
            .collect();
        s
    }
    
    fn dummy_CB_mine(n: usize) ->Vec<Vec<BusRecord>>{
        let biter2= BusReader::new(FOLDERNAME).groupby_cb();
        let s2: Vec<Vec<BusRecord>> = biter2.take(n).map(|(_a, records)|records).collect();
        s2
    }

    fn dummy_CB_mine_clone(n: usize) ->Vec<Vec<BusRecord>>{
        let biter2= BusReader::new(FOLDERNAME).groupby_cb();
        let s2: Vec<Vec<BusRecord>> = biter2.take(n).map(|(_a, records)|records).collect();
        s2
    }
    
    fn dummy_CBUMI_group(n: usize) ->Vec<Vec<BusRecord>>{
        let biter= BusReader::new(FOLDERNAME);
        let s:Vec<Vec<BusRecord>> = biter.group_by(|r| (r.CB, r.UMI)).into_iter().map(|(_a, records)|records.collect()).take(n).collect();
        s
    }
    
    fn dummy_CBUMI_mine(n: usize) ->Vec<Vec<BusRecord>>{
        let biter2= BusReader::new(FOLDERNAME).groupby_cbumi();
        let s2: Vec<Vec<BusRecord>> = biter2.take(n).map(|(_a, records)|records).collect();
        s2
    }


    c.bench_function("CB group", |b| b.iter(|| dummy_CB_group(black_box(100000))));
    c.bench_function("CB mine", |b| b.iter(|| dummy_CB_mine(black_box(100000))));
    c.bench_function("CB mine clone", |b| b.iter(|| dummy_CB_mine_clone(black_box(100000))));


    // c.bench_function("S3", |b| b.iter(|| dummy_CBUMI_group(black_box(100000))));
    // c.bench_function("S4", |b| b.iter(|| dummy_CBUMI_mine(black_box(100000))));

}

#[allow(dead_code)]
fn speed_with_capacity(c: &mut Criterion){

    fn dummy_no_capacity(n: usize) -> usize{
        let mut v = Vec::new();
        for i in 1..n {
            v.push(i);
        }
        v.len()
    }
    
    fn dummy_with_capacity(n: usize) -> usize{
        let mut v = Vec::with_capacity(n);
        for i in 1..n {
            v.push(i);
        }
        v.len()
    }

    c.bench_function("no capacity", |b| b.iter(|| dummy_no_capacity(black_box(10000))));
    c.bench_function("with capacity", |b| b.iter(|| dummy_with_capacity(black_box(10000))));

}


// fn clone_vs_ref_CB_iterators(c: &mut Criterion){

//     fn dummy_clone(n: usize) -> Vec<(u64, Vec<BusRecord>)> {
//         let foldername = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
//         let biter2= CellIteratorClone::new(foldername);
//         let v: Vec<_> = biter2.take(n).collect();
//         v
//     }

//     fn dummy_ref(n: usize) -> Vec<(u64, Vec<BusRecord>)> {
//         let foldername = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
//         let biter2= CellIterator::new(foldername);
//         let v: Vec<_> = biter2.take(n).collect();
//         v
//     }

//     c.bench_function("CB iter with clone", |b| b.iter(|| dummy_clone(black_box(10000))));
//     c.bench_function("CB iter with ref", |b| b.iter(|| dummy_ref(black_box(10000))));

// }


#[allow(dead_code)]
fn create_dummy_ec() ->Ec2GeneMapper{
    let ec0: HashSet<Genename> = vec![Genename("A".to_string())].into_iter().collect();
    let ec1: HashSet<Genename> = vec![Genename("B".to_string())].into_iter().collect();
    let ec2: HashSet<Genename> = vec![Genename("A".to_string()), Genename("B".to_string())].into_iter().collect();
    let ec3: HashSet<Genename> = vec![Genename("C".to_string()), Genename("D".to_string())].into_iter().collect();

    let ec_dict: HashMap<EC, HashSet<Genename>> = HashMap::from([
        (EC(0), ec0), // A
        (EC(1), ec1), // B
        (EC(2), ec2), // A,B
        (EC(3), ec3), // C,D
        ]);

    Ec2GeneMapper::new(ec_dict)
}

// fn setup_groupby() -> (HashMap<String, BusRecord> ,HashMap<String, &'static Ec2GeneMapper>){
//     let es1 = create_dummy_ec();
//     let es2 = create_dummy_ec();
//     let es_dict: HashMap<String, &Ec2GeneMapper> = vec![
//         ("s1".to_string(), &es1),
//         ("s2".to_string(), &es2),
//         ]
//     .into_iter().collect();

//     let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 1, FLAG: 0};
//     let s1 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 1, FLAG: 0};
//     let record_dict = vec![
//         ("s1".to_string(), r1),
//         ("s2".to_string(), s1),
//     ].into_iter().collect();

//     (record_dict, es_dict)
// }

// fn bench_groupby_genes(c: &mut Criterion){
//     // use bustools::phantompurger::groupby_gene_simple;
//     // fn dummy_simple(n: usize)->Vec<HashMap<String, BusRecord>>{
//     //     let (record_dict, es_dict) = setup_groupby();
//     //     let res = groupby_gene_simple(record_dict, &es_dict);
//     //     res
//     // }
//     // c.bench_function("sinple", |b| b.iter(|| dummy_simple(black_box(10000))));

//     use rustphantompurger::phantompurger::groupby_gene_even_simpler;
//     fn dummy_even(n: usize)->Vec<HashMap<String, BusRecord>>{
//         let (record_dict, es_dict) = setup_groupby();
//         let res = groupby_gene_even_simpler(record_dict, &es_dict);
//         res
//     }

//     c.bench_function("even simpler", |b| b.iter(|| dummy_even(black_box(10000))));
// }

#[allow(dead_code)]
fn bench_busreader_buffersize(c: &mut Criterion){

    let buffersize_vector = vec![800, 8_000, 80_000, 800_000, 8_000_000];
    fn dummy(buffersize: usize) -> usize{
        let bfile = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        return BusReader::new_with_capacity(bfile, buffersize).count();
    }
    for bsize in buffersize_vector{
        c.bench_function(&format!("Buffersize: {bsize}"), |b| b.iter(|| dummy(black_box(bsize))));
    }
}

#[allow(dead_code)]
fn bench_buswriter_buffersize(c: &mut Criterion){
    let buffersize_vector = vec![800, 8_000, 80_000, 800_000, 8_000_000];
    fn dummy(buffersize: usize) -> usize{

        let fname_nozip = "/tmp/test.bus";
        let params = BusParams{ cb_len:16, umi_len: 12 };
        let mut bw =  BusWriter::new_with_capacity( File::create(fname_nozip).unwrap(), params, buffersize);
        for i in 0..1000{
            for j in 0..10000{
                let r = BusRecord {CB: i, UMI: j*i, EC: 0, COUNT:1, FLAG: 0};
                // let r = BusRecord {CB: i, UMI: j, EC:0, COUNT:1, FLAG: 0};
                bw.write_record(&r);
            }
        }
        return 1;
    }
    for bsize in buffersize_vector{
        c.bench_function(&format!("Buffersize: {bsize}"), |b| b.iter(|| dummy(black_box(bsize))));
    }
}

fn bench_busz_compression_write(c: &mut Criterion){

    fn dummy_uncompressed() -> usize{
        let fname_nozip = "/tmp/test.bus";
        let params = BusParams{ cb_len:16, umi_len: 12 };
        let mut bw =  BusWriter::new( fname_nozip, params);
        for i in 1..501{
            for j in 1..10001{
                let r = BusRecord {CB: i, UMI: j*i, EC: 0, COUNT:1, FLAG: 0};
                // let r = BusRecord {CB: i, UMI: j, EC: (j % i) as u32, COUNT:1, FLAG: 0};
                bw.write_record(&r);
            }
        }
        return 1;
    }

    fn dummy_compressed() -> usize{
        let fname_zip = "/tmp/test.busz";
        let buszblocksize = 1000;
        let params = BusParams{ cb_len:16, umi_len: 12 };
        let mut bw =  BuszWriter::new( fname_zip, params, buszblocksize);
        
        for i in 1..501{
            for j in 1..10001{
                let r = BusRecord {CB: i, UMI: j, EC:(j % i) as u32, COUNT:1, FLAG: 0};
                bw.write_record(r);
            }
        }
        return 1;
    }

    c.bench_function(&format!("Uncompressed writing"), |b| b.iter(|| dummy_uncompressed()));
    c.bench_function(&format!("Compressed writing"), |b| b.iter(|| dummy_compressed()));
}

fn bench_busz_compression_read(c: &mut Criterion){
    
    fn dummy_uncompressed(fname :&str) -> usize{

        let r = BusReader::new(fname);
        let mut counter = 0;
        for record in r {
            counter += record.COUNT as usize;
        }
        counter
    }

    fn dummy_compressed(fname :&str) -> usize{

        let r = BuszReader::new(fname);
        let mut counter = 0;
        for record in r {
            counter += record.COUNT as usize;
        }
        counter
    }
    let fname_uncomp = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";
    let fname_comp = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.busz";

    c.bench_function(&format!("Uncompressed reading"), |b| b.iter(|| dummy_uncompressed(black_box(fname_uncomp))));
    c.bench_function(&format!("Compressed reading"), |b| b.iter(|| dummy_compressed(black_box(fname_comp))));

}

// fn bench_busreader_gzip(c: &mut Criterion){
//     let fname_zip = "/tmp/test.bus.gz";

//     let header = BusHeader::new(16, 12, 20);
//     let mut bw =  BusWriter::from_gzip(fname_zip, header);
//     println!("Writing big record");
//     for i in 0..1000{
//         for j in 0..1000{
//             let r = BusRecord {CB: i, UMI: j, EC:0, COUNT:1, FLAG: 0};
//             bw.write_record(&r);
//         }
//     }
//     drop(bw);

//     let fname_nozip = "/tmp/test.bus";
//     let header = BusHeader::new(16, 12, 20);
//     let mut bw =  BusWriter::new( fname_nozip, header);
//     println!("Writing big record");
//     for i in 0..1000{
//         for j in 0..1000{
//             let r = BusRecord {CB: i, UMI: j, EC:0, COUNT:1, FLAG: 0};
//             bw.write_record(&r);
//         }
//     }
//     drop(bw);

//     fn dummy_gzip(fname: &str) -> usize{
//         let reader = BusReader::from_gzip(fname);
//         let n = reader.count();
//         n
//     }
//     fn dummy_no_gzip(fname: &str) -> usize{
//         let reader = BusReader::new(fname);
//         let n = reader.count();
//         n
//     }
//     c.bench_function(&format!("Read CB gzipped"), |b| b.iter(|| dummy_gzip(black_box(fname_zip))));
//     c.bench_function(&format!("Read CB unzipped"), |b| b.iter(|| dummy_no_gzip(black_box(fname_nozip))));

// }

#[allow(dead_code)]
fn busmulti_vs_merger(c: &mut Criterion){
    let fname1 = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
    let fname2 = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
    let n = 100000;

    fn dummy_busmulti(fname1: &str,fname2: &str, n:usize) -> usize{

        // let b1 = BusReader::new(fname1).groupby_cb();
        // let b2 = BusReader::new(fname2).groupby_cb();
        let hashmap = HashMap::from([
            ("test1".to_string(), fname1.to_string()),
            ("test2".to_string(), fname2.to_string()),
        ]);
        let iii = CellIteratorMulti::new(&hashmap);
        let v: Vec<_> = iii.take(n).collect();
        v.len()
    }

    fn dummy_merger(fname1: &str,fname2: &str, n:usize) -> usize{

        let b1 = BusReader::new(fname1).groupby_cb();
        let b2 = BusReader::new(fname2).groupby_cb();
        let hashmap = HashMap::from([
            ("test1".to_string(), b1),
            ("test2".to_string(), b2),
        ]);
        let iii = MultiIterator::new(hashmap);
        let v: Vec<_> = iii.take(n).collect();
        v.len()
    }

    c.bench_function("busmulti", |b| b.iter(|| dummy_busmulti(black_box(fname1), black_box(fname2), n )));
    c.bench_function("merger"  , |b| b.iter(|| dummy_merger(black_box(fname1), black_box(fname2), n )));

}

use bitvec::prelude::{self as bv, Msb0};
use rand::distributions::Uniform;
use rand::prelude::Distribution;
pub fn bitslice_to_bytes(bits: &bv::BitSlice<u8, bv::Msb0>) -> Vec<u8>{

    assert_eq!(bits.len() % 8,  0, "cant covnert to bytes if Bitsclie is not a multiple of 8");

    let nbytes = bits.len() / 8;
    let mut bytes = Vec::with_capacity(nbytes);
    for i in 0..nbytes {
        let b = &bits[i*8.. (i+1)*8];
        let a: u8 = b.load_be(); // doesnt matter if be/le here since its only 1byte anyway
        bytes.push(a);
    }
    bytes
}

pub fn swap_endian(bytes: &[u8], wordsize: usize) -> Vec<u8>{
    let mut swapped_endian: Vec<u8> = Vec::with_capacity(bytes.len());
    for bytes in bytes.chunks(wordsize){
        swapped_endian.extend(bytes.iter().rev());
    }
    swapped_endian
}
fn swap_vecu8_vs_bistream(c: &mut Criterion){
    let n_elements = 1000;
    let data_dist = Uniform::from(0..255);
    let mut rng = rand::thread_rng();
    let mut data: Vec<u8> = Vec::with_capacity(n_elements);
    for _ in 0..n_elements {
        data.push(data_dist.sample(&mut rng) as u8);
    }

    let bits = bv::BitVec::from_slice(&data);

    fn convert_and_swap(bits: bv::BitVec<u8, bv::Msb0>) -> bv::BitVec<u8, bv::Msb0>  {
        let x = bitslice_to_bytes(&bits);
        let swapped = swap_endian(&x, 8);
        let swapped_bits = bv::BitVec::from_slice(&swapped);
        swapped_bits
    }

    fn direct_swap(bits: bv::BitVec<u8, Msb0>) -> bv::BitVec<u8, bv::Msb0>  {

        let mut res: bv::BitVec<u8, bv::Msb0> = bv::BitVec::with_capacity(bits.len());
        for c in bits.chunks(8*8) {
            // let mut owned = c.to_bitvec();
            // owned.reverse();
            res.extend(c.iter().rev());
        }
        res
    }

    c.bench_function(&format!("convert and swap"), |b| b.iter(|| convert_and_swap(black_box(bits.clone()))));
    c.bench_function(&format!("direct swap"), |b| b.iter(|| direct_swap(black_box(bits.clone()))));

    // let x_swapped = 
}


// criterion_group!(benches, criterion_benchmark);
// criterion_group!(benches, bench_buswriter_buffersize);
// criterion_group!(benches, busmulti_vs_merger);
// criterion_group!(benches, plain_iterator_speed);
// criterion_group!(benches, bench_busz_compression_write);
criterion_group!(benches, bench_busz_compression_read);
// criterion_group!(benches, swap_vecu8_vs_bistream);

criterion_main!(benches);