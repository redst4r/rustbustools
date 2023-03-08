use crate::{io::BusIteratorBuffered, iterators::{CellGroupIterator, CbUmiGroupIterator}};
use crate::io::{BusRecord, setup_busfile};


pub fn inspect(busfile: &str){

    let n_cbumi = BusIteratorBuffered::new(busfile).groupby_cbumi().count();
    let n_cells =  BusIteratorBuffered::new(busfile).groupby_cb().count();

    let mut nreads = 0;
    let mut nrecords = 0;

    let bus = BusIteratorBuffered::new(busfile);
    for r in bus{
        nrecords+=1;
        nreads+= r.COUNT
    }

    println!("{nrecords} BUS records");
    println!("{nreads} reads");
    println!("{n_cells} cell-barcodes");
    println!("{n_cbumi} CB-UMIs");
}

#[test]
fn test_inspect(){
    let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
    let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
    let r3 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
    let r4 = BusRecord{CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
    let r5 = BusRecord{CB: 2, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
    let r6 = BusRecord{CB: 3, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

    let records = vec![r1.clone(),r2.clone(),r3.clone(),r4.clone(),r5.clone(), r6.clone()];
    // let records = vec![r1,r2,r3,r4,r5, r6].to_vec();

    let (busname, _dir) = setup_busfile(&records);

    inspect(&busname);

}