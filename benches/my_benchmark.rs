use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufReader, Read};

use bitvec::field::BitField;
use bustools::busz::{BuszWriter, BuszReader};
use criterion::{black_box, criterion_group, criterion_main, Criterion};
use itertools::{izip, Itertools};
use bustools::bus_multi::CellIteratorMulti;
use bustools::consistent_genes::{Ec2GeneMapper, Genename, EC};
use bustools::io::{BusReaderPlain, BusRecord, BusWriterPlain, BusParams, setup_busfile};
use bustools::iterators::{CbUmiGroup, CbUmiGroupIterator, CellGroup, CellGroupIterator};
use bustools::merger::{MultiIteratorSlow, MultiIterator};


/// puls a busfile into memory; careful for large files!
fn bus_to_mem(busfile: &str) -> Vec<u8>{
    let mut buffer = Vec::new();
    let mut f= File::open(busfile).unwrap();
    f.read_to_end(&mut buffer).unwrap();
    buffer
}
#[allow(dead_code)]
/// profiling the deserializing of the busfile more or less
/// based on IN MEMORY busfiles to get around I/O
fn plain_and_compressed_iterator_speed(c: &mut Criterion){
    
    // let busname = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
    let busname = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
    let busname_compressed = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.busz";

    // pull the file into memory, to alleviate disk access
    let buffer = bus_to_mem(busname);
    let n = 1_000_000;

    c.bench_function("plain_iterator in memory",
     |b| b.iter(|| {
        BusReaderPlain::from_read(buffer.as_slice()).take(n).map(|r|r.COUNT).sum::<u32>();
     }
    ));

    // // do it on disk
    // c.bench_function("plain_iterator on disk",
    //  |b| b.iter(|| {
    //     BusReaderPlain::new(busname).take(n).map(|r|r.COUNT).sum::<u32>();
    //  }
    // ));

    // Compressed iterators
    let buffer_compressed = bus_to_mem(busname_compressed);

    c.bench_function("compressed_iterator in memory",
    |b| b.iter(|| {
       BuszReader::from_read(buffer_compressed.as_slice()).take(n).map(|r|r.COUNT).sum::<u32>();
    }
    ));

    // // do it on disk
    // c.bench_function("compressed_iterator on disk",
    //  |b| b.iter(|| {
    //     BuszReader::new(busname_compressed).take(n).map(|r|r.COUNT).sum::<u32>();
    //  }
    // ));


    c.bench_function("u64 compressed_iterator in memory",
    |b| b.iter(|| {
       bustools::busz::decode2::BuszReader::from_read(buffer_compressed.as_slice()).take(n).map(|r|r.COUNT).sum::<u32>();
    }
    ));

}


/* testing the speed of my iterators (CB/UMI) against the grouped version
 */
#[allow(dead_code)]
fn iterator_speed(c: &mut Criterion){

     pub const FOLDERNAME:&str = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
    
    fn dummy_cb_mine(n: usize) ->Vec<Vec<BusRecord>>{
        let biter2= BusReaderPlain::new(FOLDERNAME).groupby_cb();
        let s2: Vec<Vec<BusRecord>> = biter2.take(n).map(|(_a, records)|records).collect();
        s2
    }
    
    fn dummy_cbumi_mine(n: usize) ->Vec<Vec<BusRecord>>{
        let biter2= BusReaderPlain::new(FOLDERNAME).groupby_cbumi();
        let s2: Vec<Vec<BusRecord>> = biter2.take(n).map(|(_a, records)|records).collect();
        s2
    }

    c.bench_function("CB mine", |b| b.iter(|| dummy_cb_mine(black_box(100000))));
    c.bench_function("CBUMI mine", |b| b.iter(|| dummy_cbumi_mine(black_box(100000))));
}

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
        // return BusReader::new_with_capacity(bfile, buffersize).count();
        let reader = BufReader::with_capacity(buffersize, File::open(bfile).unwrap());
        BusReaderPlain::from_read(reader).count()

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
        let mut bw =  BusWriterPlain::new_with_capacity( File::create(fname_nozip).unwrap(), params, buffersize);
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

#[allow(dead_code)]
fn bench_busz_compression_write(c: &mut Criterion){

    fn dummy_uncompressed() -> usize{
        let fname_nozip = "/tmp/test.bus";
        let params = BusParams{ cb_len:16, umi_len: 12 };
        let mut bw =  BusWriterPlain::new( fname_nozip, params);
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

#[allow(dead_code)]
fn bench_busz_compression_read(c: &mut Criterion){
    
    fn dummy_uncompressed(fname :&str) -> usize{

        let r = BusReaderPlain::new(fname);
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


#[allow(dead_code)]
fn busmulti_vs_merger(c: &mut Criterion){
    let fname1 = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
    let fname2 = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
    let n = 100000;

    fn dummy_busmulti(fname1: &str,fname2: &str, n:usize) -> usize{

        let b1 = BusReaderPlain::new(fname1).groupby_cb();
        let b2 = BusReaderPlain::new(fname2).groupby_cb();
        let hashmap = HashMap::from([
            ("test1".to_string(), b1),
            ("test2".to_string(), b2),
        ]);
        let iii = CellIteratorMulti::new(hashmap);
        let v: Vec<_> = iii.take(n).collect();
        v.len()
    }

    fn dummy_merger(fname1: &str,fname2: &str, n:usize) -> usize{

        let b1 = BusReaderPlain::new(fname1).groupby_cb();
        let b2 = BusReaderPlain::new(fname2).groupby_cb();
        let hashmap = HashMap::from([
            ("test1".to_string(), b1),
            ("test2".to_string(), b2),
        ]);
        let iii = MultiIteratorSlow::new(hashmap);
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

#[allow(dead_code)]
fn bench_merger(c: &mut Criterion){

    let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0 };
    let r2 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
    let r3 = BusRecord { CB: 1, UMI: 3, EC: 0, COUNT: 2, FLAG: 0 };
    let r4 = BusRecord { CB: 3, UMI: 0, EC: 0, COUNT: 2, FLAG: 0 };

    let v1 = vec![r1.clone(), r2.clone(), r3.clone(), r4.clone()];

    let s1 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
    let s2 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 12, FLAG: 0 };
    let s3 = BusRecord { CB: 1, UMI: 3, EC: 1, COUNT: 12, FLAG: 0 };
    let s4 = BusRecord { CB: 2, UMI: 3, EC: 1, COUNT: 2, FLAG: 0 };
    let s5 = BusRecord { CB: 2, UMI: 3, EC: 2, COUNT: 2, FLAG: 0 };
    let v2 = vec![s1.clone(), s2.clone(), s3.clone(), s4.clone(), s5.clone()];

    // write the records to file
    let (busname1, _dir1) = setup_busfile(&v1);
    let (busname2, _dir2) = setup_busfile(&v2);

    fn dummy_merger(fname1: &str,fname2: &str, n:usize) -> usize{

        let b1 = BusReaderPlain::new(fname1).groupby_cb();
        let b2 = BusReaderPlain::new(fname2).groupby_cb();
        let hashmap = HashMap::from([
            ("test1".to_string(), b1),
            ("test2".to_string(), b2),
        ]);
        let iii = MultiIteratorSlow::new(hashmap);
        let v: Vec<_> = iii.take(n).collect();
        v.len()
    }

    fn dummy_merger_fast(fname1: &str,fname2: &str, n:usize) -> usize{

        let b1 = BusReaderPlain::new(fname1).groupby_cb();
        let b2 = BusReaderPlain::new(fname2).groupby_cb();
        let hashmap = HashMap::from([
            ("test1".to_string(), b1),
            ("test2".to_string(), b2),
        ]);
        let iii = MultiIterator::new(hashmap);
        let v: Vec<_> = iii.take(n).collect();
        v.len()
    }

    let n = 1000;
    c.bench_function(
        "merger",
        |b| b.iter(
            || dummy_merger(
                black_box(&busname1), 
                black_box(&busname2), 
                n)
        ));

        c.bench_function(
            "merger fast",
            |b| b.iter(
                || dummy_merger_fast(
                    black_box(&busname1), 
                    black_box(&busname2), 
                    n)
            ));
}

#[allow(dead_code)]
fn bench_merger_real(c: &mut Criterion){
    let busname1 = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
    let busname2 = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";

    let n = 100000;

    // lets pull those into mem, so that its not just disk access we measure
    let b1 = BusReaderPlain::new(busname1).groupby_cbumi().take(n);
    let b2 = BusReaderPlain::new(busname2).groupby_cbumi().take(n);

    let r1 = b1.collect::<Vec<_>>();
    let r2 = b2.collect::<Vec<_>>();

    fn dummy_merger_mem(r1: Vec<((u64, u64), Vec<BusRecord>)>,r2:Vec<((u64, u64), Vec<BusRecord>)>) -> usize{

        let hashmap = HashMap::from([
            ("test1".to_string(), r1.into_iter()),
            ("test2".to_string(), r2.into_iter()),
        ]);
        let iii = MultiIteratorSlow::new(hashmap);
        let v: Vec<_> = iii.collect();
        v.len()
    }

    fn dummy_merger_fast_mem(r1: Vec<((u64, u64), Vec<BusRecord>)>,r2:Vec<((u64, u64), Vec<BusRecord>)>) -> usize{

        let hashmap = HashMap::from([
            ("test1".to_string(), r1.into_iter()),
            ("test2".to_string(), r2.into_iter()),
        ]);
        let iii = MultiIterator::new(hashmap);
        let v: Vec<_> = iii.collect();
        v.len()
    }

    fn dummy_merger(fname1: &str,fname2: &str, n:usize) -> usize{

        let b1 = BusReaderPlain::new(fname1).groupby_cbumi();
        let b2 = BusReaderPlain::new(fname2).groupby_cbumi();
        let hashmap = HashMap::from([
            ("test1".to_string(), b1),
            ("test2".to_string(), b2),
        ]);
        let iii = MultiIteratorSlow::new(hashmap);
        let v: Vec<_> = iii.take(n).collect();
        v.len()
    }

    fn dummy_merger_fast(fname1: &str,fname2: &str, n:usize) -> usize{

        let b1 = BusReaderPlain::new(fname1).groupby_cbumi();
        let b2 = BusReaderPlain::new(fname2).groupby_cbumi();
        let hashmap = HashMap::from([
            ("test1".to_string(), b1),
            ("test2".to_string(), b2),
        ]);
        let iii = MultiIterator::new(hashmap);
        let v: Vec<_> = iii.take(n).collect();
        v.len()
    }

    // c.bench_function(
    //     "merger",
    //     |b| b.iter(
    //         || dummy_merger(
    //             black_box(&busname1), 
    //             black_box(&busname2), 
    //             n)
    //     ));

    // c.bench_function(
    //     "merger fast",
    //     |b| b.iter(
    //         || dummy_merger_fast(
    //             black_box(&busname1), 
    //             black_box(&busname2), 
    //             n)
    //     ));    
    c.bench_function(
        "merger",
        |b| b.iter(
            || dummy_merger_mem(
                black_box(r1.clone()), 
                black_box(r2.clone()), 
                )
        ));

    c.bench_function(
        "merger fast",
        |b| b.iter(
            || dummy_merger_fast_mem(
                black_box(r1.clone()), 
                black_box(r2.clone()), 
                )
        ));    
}


use ahash::{AHasher, RandomState};

fn bench_hashmaps(c: &mut Criterion){

    let keys = vec!["ada", "gfsg", "sgrs","hsdfs", "dsf","asd","favm"];
    let values = vec![1,2,3,4,5,6,7];

    c.bench_function(
        "ahash fast",
        |b| b.iter(
            || {
                let mut h : HashMap<String, usize, RandomState> = HashMap::default();

                for (k, v) in izip!(keys.iter(), values.iter()) {
                    h.insert(k.to_string(), *v);
                }
            }
        ));  
    c.bench_function(
        "regular hashmap",
        |b| b.iter(
            || {
                let mut h : HashMap<String, usize> = HashMap::new();

                for (k, v) in izip!(keys.iter(), values.iter()) {
                    h.insert(k.to_string(), *v);
                }
            }
        ));  

}


/// grouping on plain and compressed iterators
fn grouped_plain_and_compressed_iterator_speed(c: &mut Criterion){
    
    // let busname = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
    let busname = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
    let busname_compressed = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.busz";

    // pull the file into memory, to alleviate disk access
    let buffer = bus_to_mem(busname);
    let n = 100_000;

    c.bench_function("plain_iterator groupby_cb in memory",
     |b| b.iter(|| {
        BusReaderPlain::from_read(buffer.as_slice()).groupby_cb().take(n).map(|(_cb, rlist)|rlist.len()).sum::<usize>();
     }
    ));


    // Compressed iterators
    let buffer_compressed = bus_to_mem(busname_compressed);

    c.bench_function("compressed_iterator in memory",
    |b| b.iter(|| {
       BuszReader::from_read(buffer_compressed.as_slice()).groupby_cb().take(n).map(|(_cb, rlist)|rlist.len()).sum::<usize>();;
    }
    ));

}


// criterion_group!(benches, criterion_benchmark);
// criterion_group!(benches, bench_buswriter_buffersize);
// criterion_group!(benches, busmulti_vs_merger);
// criterion_group!(benches, iterator_speed);
criterion_group!(benches, plain_and_compressed_iterator_speed);
// criterion_group!(benches, grouped_plain_and_compressed_iterator_speed);

// criterion_group!(benches, bench_io_vs_inmem);


// criterion_group!(benches, bench_busz_compression_write);
// criterion_group!(benches, bench_busz_compression_read);
// criterion_group!(benches, swap_vecu8_vs_bistream);

// criterion_group!(benches, bench_merger_real);
// criterion_group!(benches, bench_hashmaps);

criterion_main!(benches);