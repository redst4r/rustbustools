mod bus;

use crate::bus::{BusRecord, BusHeader, read_bus_header, read_bus_records, BusIteratorBuffered};

use std::collections::HashMap;

fn main() {
    // trying out serializing and deserialzing
    let r = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};

    // serialize
    let rbin = bincode::serialize(&r).unwrap();
    println!("struct Point serializes into byte array {:?}", rbin);

    // deserialize
    let deser: BusRecord = bincode::deserialize(&rbin).unwrap();
    println!("deser into {:?}", deser);


    // reading a busfile header
    let fname = "/home/michi/ms_python_packages/rustbustools/some.bus";
    let header_struct = read_bus_header(fname);
    println!("Header {:?}", header_struct);

    // read the busrecords
    let records = read_bus_records(fname);
    for r in records{
        println!("{:?}", r)
    }
    let mut bi = BusIteratorBuffered::new(&fname.to_string());
    println!("BI {:?}", bi);
    let record = bi.next();
    println!("{:?}", record);
    let record = bi.next();
    println!("{:?}", record);
    let record = bi.next();
    println!("{:?}", record);

    // assert_eq!(BI.next(), None);
    // let fh = std::fs::File::open(fname).expect("FAIL");
    // BufReader::new(fh)

    println!("===================");
    let mut counter: HashMap<u64, u32> = HashMap::new();

    // scp scp -r michi@baker:~/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/Fresh1/kallisto/sort_bus/bus_output/output.corrected.sort.bus /tmp/test.bus
    let fname = "/home/michi/output.corrected.sort.bus";
    let BI = BusIteratorBuffered::new(&fname.to_string());
    // let BI = BusIterator::new(&fname.to_string());
    println!("{:?}", BI.bus_header);

    for (i, r) in BI.enumerate(){
        // println!("{:?}", r);
        // break;
        let count = counter.entry(r.CB).or_insert(0);
        *count += 1;

        if i % 1_000_000 == 0{
            println!("Iteration {} Mio", i/ 1_000_000)
        }
    }
    println!("{:?}", counter.len());

    // let mut fh = std::fs::File::open("/home/michi/ms_python_packages/rustbustools/test.txt").expect("FAIL");
    // // let mut buffer = Vec::with_capacity(32);
    // let mut buffer = [0;32];
    //
    // println!("XXX Buffer {} {:?}",buffer.len(), buffer);
    //
    // let n = fh.read(&mut buffer).expect("fail");
    // println!("Buffer {} {:?}",n, buffer);


}
