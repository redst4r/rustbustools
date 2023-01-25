use crate::io::{BusIteratorBuffered, BusRecord};

pub struct CbUmiIterator {
    pub(crate) busiter: BusIteratorBuffered,
    pub(crate) last_record: Option<BusRecord>  //option needed to mark the final element of the iteration
}

impl CbUmiIterator {

    pub fn new(fname: &str) ->CbUmiIterator{
        let mut busiter = BusIteratorBuffered::new(fname);
        let last_record = busiter.next(); //initilize with the first record in the file
        CbUmiIterator {busiter, last_record}
    }
}

impl Iterator for CbUmiIterator {
    type Item = ((u64, u64), Vec<BusRecord>);

    fn next(&mut self) -> Option<Self::Item> {

        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted
    
        loop {
            if let Some(last_record) = self.last_record {  //if we're not done with the iteration
                // try to get a new record
                if let Some(new_record) = self.busiter.next(){

                    // determine if we encounter a new cell in this iteration
                    let (current_cb, current_umi) = (last_record.CB, last_record.UMI);

                    // last_record is overhauled now (since we got a "new last record")
                    // and we add it 
                    busrecords.push(last_record); // the stored element from the previous iteration
                    self.last_record = Some(new_record);

                    // now we just need to decide if we want to emit, or continue growing
                    if new_record.CB > current_cb || (new_record.CB == current_cb &&  new_record.UMI > current_umi){  
                        // we ran into a new CB/UMI and its records
                        // println!("\tyielding {:?}", (current_cb, &busrecords));

                        return Some(((current_cb, current_umi), busrecords));
                    }
                    else if (new_record.CB == current_cb) && new_record.UMI == current_umi {
                        // nothing happens, just keep growing busrecords
                    }
                    else{  // the new cb is smaller then the current state: this is a bug due to an UNOSORTED busfile
                        panic!("Unsorted busfile: {}/{} -> {}/{}", current_cb, current_umi, new_record.CB, new_record.UMI)
                    }
                }
                else {
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    busrecords.push(last_record);
                    let current_cb = last_record.CB;
                    let current_umi = last_record.UMI;
                    // to mark the end of iteration and all items emitted, set last_item to None
                    self.last_record = None;
                    return Some(((current_cb, current_umi), busrecords));  
                }
            }
            else{  // last_record == None
                // we are done
                return None
            }
        }
    }
}

//=================================================================================
pub struct CellIterator {
    pub(crate) busiter: BusIteratorBuffered,
    pub(crate) last_record: Option<BusRecord>  //option needed to mark the final element of the iteration
}

impl CellIterator {

    pub fn new(fname: &str) ->CellIterator{
        let mut busiter = BusIteratorBuffered::new(fname);
        let last_record = busiter.next(); //initilize with the first record in the file
        CellIterator {busiter, last_record}
    }
}

impl Iterator for CellIterator {
    type Item = (u64, Vec<BusRecord>);

    fn next(&mut self) -> Option<Self::Item> {

        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted
    
        loop {
            if let Some(last_record) = self.last_record{  //if we're not done with the iteration
                // try to get a new record
                if let Some(new_record) = self.busiter.next(){

                    // determine if we encounter a new cell in this iteration
                    let current_cb = last_record.CB;

                    // last_record is overhauled now (since we got a "new last record")
                    // and we add it 
                    busrecords.push(last_record); // the stored element from the previous iteration
                    self.last_record = Some(new_record);


                    if new_record.CB > current_cb {  
                        // we ran into a new CB and its records
                        // println!("\tyielding {:?}", (current_cb, &busrecords));
                        return Some((current_cb, busrecords));
                    }
                    else if new_record.CB == current_cb {
                        // nothing
                    }
                    else{  // the new cb is smaller then the current state: this is a bug due to an UNOSORTED busfile
                        panic!("Unsorted busfile: {} -> {}", current_cb, new_record.CB)
                    }
                }
                else {
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    busrecords.push(last_record);
                    let current_cb = last_record.CB;
                    // to mark the end of iteration and all items emitted, set last_item to None
                    self.last_record = None;
                    return Some((current_cb, busrecords));  
                }
            }
            else{  // last_record == None
                // we are done
                return None
            }
        }
    }
}



 mod tests{
    use crate::io::{BusRecord, setup_busfile};
    use crate::iterators::{CbUmiIterator, CellIterator};
    #[test]
    fn test_cb_iter(){   
        let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r4 = BusRecord{CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r5 = BusRecord{CB: 2, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r6 = BusRecord{CB: 3, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
    
        let records = vec![r1,r2,r3,r4,r5, r6];
    
        let busname = "/tmp/test_iter.bus";
        setup_busfile(&records, &busname);
    
    
        let cb_iter = CellIterator::new(busname);
        // let n: Vec<Vec<BusRecord>> = cb_iter.map(|(_cb, rec)| rec).collect();
        let n: Vec<(u64, Vec<BusRecord>)> = cb_iter.collect();
        // println!("{:?}", n);
    
        assert_eq!(n.len(), 4);
        // println!("{:?}", n);
        // println!("First");
        let c1 = &n[0];
        assert_eq!(*c1, (0, vec![r1,r2]));
    
        // println!("Second");
        let c2 = &n[1];
        assert_eq!(*c2, (1,vec![r3]));
    
        // println!("Third");
        let c3 = &n[2];
        assert_eq!(*c3, (2, vec![r4,r5]));
    
        let c4 = &n[3];
        assert_eq!(*c4, (3, vec![r6]));
    
        // assert_eq!(n, vec![vec![r1,r2], vec![r3], vec![r4,r5]])
    
     }
    #[test]
    fn test_cbumi_iter(){   
        let r1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r4 = BusRecord{CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r5 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0};
        let r6 = BusRecord{CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

        let records = vec![r1,r2,r3,r4,r5, r6];

        let busname = "/tmp/test_cbumi_iter.bus";
        setup_busfile(&records, &busname);


        let cb_iter = CbUmiIterator::new(busname);
        // let n: Vec<Vec<BusRecord>> = cb_iter.map(|(_cb, rec)| rec).collect();
        let n: Vec<((u64, u64), Vec<BusRecord>)> = cb_iter.collect();
        println!("{:?}", n);

        assert_eq!(n.len(), 5);
        // println!("{:?}", n);
        // println!("First");
        let c1 = &n[0];
        assert_eq!(*c1, ((0,1 ), vec![r1,r2]));

        // println!("Second");
        let c2 = &n[1];
        assert_eq!(*c2, ((0, 2), vec![r3]));

        // println!("Third");
        let c3 = &n[2];
        assert_eq!(*c3, ((1,1), vec![r4]));

        let c4 = &n[3];
        assert_eq!(*c4, ((1,2), vec![r5]));

        let c5 = &n[4];
        assert_eq!(*c5, ((2,1), vec![r6]));
        // assert_eq!(n, vec![vec![r1,r2], vec![r3], vec![r4,r5]])

    }
    #[test]
    #[should_panic(expected = "Unsorted busfile: 2 -> 0")]
    fn test_panic_on_unsorted(){  
        let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 2, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r4 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

        let records = vec![r1,r2,r3,r4];

        let busname = "/tmp/test_panic_on_unsorted.bus";
        setup_busfile(&records, busname);

        let cb_iter = CellIterator::new(busname);
        let _n: Vec<(u64, Vec<BusRecord>)> = cb_iter.collect();
    }
 }