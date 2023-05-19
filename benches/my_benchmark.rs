use std::collections::{HashMap, HashSet};
use std::fs::File;

use criterion::{black_box, criterion_group, criterion_main, Criterion};
use itertools::Itertools;
use rustbustools::bus_multi::CellIteratorMulti;
use rustbustools::consistent_genes::{Ec2GeneMapper, Genename, EC};
use rustbustools::io::{BusReader, BusRecord, BusWriter, BusHeader};
use rustbustools::iterators::{CellGroup, CbUmiGroup, CellGroupIterator, CbUmiGroupIterator};
use rustbustools::merger::MultiIterator;
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
        let s:Vec<Vec<BusRecord>> = biter.group_by(|r| (r.CB, r.UMI)).into_iter().map(|(a, records)|records.collect()).take(n).collect();
        s
    }
    
    fn dummy_CBUMI_mine(n: usize) ->Vec<Vec<BusRecord>>{
        let biter2= BusReader::new(FOLDERNAME).groupby_cbumi();
        let s2: Vec<Vec<BusRecord>> = biter2.take(n).map(|(a, records)|records).collect();
        s2
    }


    c.bench_function("CB group", |b| b.iter(|| dummy_CB_group(black_box(100000))));
    c.bench_function("CB mine", |b| b.iter(|| dummy_CB_mine(black_box(100000))));
    c.bench_function("CB mine clone", |b| b.iter(|| dummy_CB_mine_clone(black_box(100000))));


    // c.bench_function("S3", |b| b.iter(|| dummy_CBUMI_group(black_box(100000))));
    // c.bench_function("S4", |b| b.iter(|| dummy_CBUMI_mine(black_box(100000))));

}
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


fn bench_buswriter_buffersize(c: &mut Criterion){
    let buffersize_vector = vec![800, 8_000, 80_000, 800_000, 8_000_000];
    fn dummy(buffersize: usize) -> usize{

        let fname_nozip = "/tmp/test.bus";
        let header = BusHeader::new(16, 12, 20);
        let mut bw =  BusWriter::new_with_capacity( File::create(fname_nozip).unwrap(), header, buffersize);
        for i in 0..1000{
            for j in 0..10000{
                let r = BusRecord {CB: i, UMI: j, EC:0, COUNT:1, FLAG: 0};
                bw.write_record(&r);
            }
        }
        return 1;

    }
    for bsize in buffersize_vector{
        c.bench_function(&format!("Buffersize: {bsize}"), |b| b.iter(|| dummy(black_box(bsize))));
    }
}



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


// criterion_group!(benches, criterion_benchmark);
// criterion_group!(benches, bench_buswriter_buffersize);
criterion_group!(benches, busmulti_vs_merger);
criterion_main!(benches);