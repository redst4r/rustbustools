use crate::io::{BusRecord, BusIteratorBuffered};

pub struct CbUmiIterator_clone {
    pub(crate) busiter: BusIteratorBuffered,
    pub(crate) last_record: Option<BusRecord>  //option needed to mark the final element of the iteration
}

impl CbUmiIterator_clone {

    pub fn new(fname: &str) ->CbUmiIterator_clone{
        let mut busiter = BusIteratorBuffered::new(fname);
        let last_record = busiter.next(); //initilize with the first record in the file
        CbUmiIterator_clone {busiter, last_record}
    }
}

impl Iterator for CbUmiIterator_clone {
    type Item = ((u64, u64), Vec<BusRecord>);


    fn next(&mut self) -> Option<Self::Item> {

        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted
    
        loop {
            
            if let Some(last_record) = self.last_record.clone() {  //if we're not done with the iteration

                // we need to move last_record out of the struct and get ownership of it (since we're emitting it within this next call)

                
                // try to get a new record
                if let Some(new_record) = self.busiter.next(){

                    // determine if we encounter a new cell in this iteration
                    let (current_cb, current_umi) = (last_record.CB, last_record.UMI);

                    // last_record is overhauled now (since we got a "new last record")
                    // and we add it 
                    busrecords.push(last_record); // the stored element from the previous iteration

                    let (newcb, newumi) = (new_record.CB, new_record.UMI);
                    self.last_record = Some(new_record);

                    // now we just need to decide if we want to emit, or continue growing
                    if newcb > current_cb || (newcb == current_cb &&  newumi > current_umi){  
                        // we ran into a new CB/UMI and its records
                        // println!("\tyielding {:?}", (current_cb, &busrecords));

                        return Some(((current_cb, current_umi), busrecords));
                    }
                    else if (newcb == current_cb) && newumi == current_umi {
                        // nothing happens, just keep growing busrecords
                    }
                    else{  // the new cb is smaller then the current state: this is a bug due to an UNOSORTED busfile
                        // panic!("Unsorted busfile: {}/{} -> {}/{}", current_cb, current_umi, new_record.CB, new_record.UMI)
                        panic!("Unsorted busfile")
                    }
                }
                else {
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    let current_cb = last_record.CB;
                    let current_umi = last_record.UMI;
                    busrecords.push(last_record);

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


pub struct CellIterator_clone {
    pub(crate) busiter: BusIteratorBuffered,
    pub(crate) last_record: Option<BusRecord>,  //option needed to mark the final element of the iteration
    // buffersize: usize
}

impl CellIterator_clone {

    pub fn new(fname: &str) ->CellIterator_clone{
        let mut busiter = BusIteratorBuffered::new(fname);
        let last_record = busiter.next(); //initilize with the first record in the file
        CellIterator_clone {busiter, last_record}
        // CellIterator_clone {busiter, last_record, buffersize:1}
    }
}

impl Iterator for CellIterator_clone {
    type Item = (u64, Vec<BusRecord>);

    fn next(&mut self) -> Option<Self::Item> {

        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted
    
        loop {
            if let Some(last_record) = self.last_record.clone(){  //if we're not done with the iteration
                // try to get a new record
                if let Some(new_record) = self.busiter.next(){

                    // determine if we encounter a new cell in this iteration
                    let current_cb = last_record.CB;

                    // last_record is overhauled now (since we got a "new last record")
                    // and we add it 
                    busrecords.push(last_record); // the stored element from the previous iteration
                    let newcb = new_record.CB;

                    self.last_record = Some(new_record);


                    if newcb > current_cb {  
                        // we ran into a new CB and its records
                        // println!("\tyielding {:?}", (current_cb, &busrecords));
                        return Some((current_cb, busrecords));
                    }
                    else if newcb == current_cb {
                        // nothing
                    }
                    else{  // the new cb is smaller then the current state: this is a bug due to an UNOSORTED busfile
                        panic!("Unsorted busfile: {} -> {}", current_cb, newcb)
                    }
                }
                else {
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    let current_cb = last_record.CB;
                    busrecords.push(last_record);
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
