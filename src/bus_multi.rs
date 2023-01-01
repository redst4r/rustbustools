use std::collections::HashMap;
use crate::io::{CellIterator, BusRecord, CbUmiIterator};


pub struct CellIteratorMulti {
    pub iterators: HashMap<String, CellIterator>,
    // last_record: BusRecord
    pub current_items: HashMap<String, (u64, Vec<BusRecord>)> // filename - > (CB, ListOfRecords)
}
impl CellIteratorMulti {

    pub fn new(fnames: &HashMap<String, String>) ->CellIteratorMulti{

        // let iterators: HashMap<String, CellIterator> = fnames.iter()
        //                                              .map(|(name, fname)| (name.clone(), CellIterator::new(fname)))
        //                                              .collect();

        let mut iterators: HashMap<String, CellIterator> = HashMap::new();
        let mut current_items: HashMap<String, (u64, Vec<BusRecord>)> = HashMap::new();

        for (name, fname) in fnames{
            // create new cell iterator
            let mut the_iter = CellIterator::new(fname);

            // populate first elements from that iterator
            let item = the_iter.next();
            match item{
                Some(record_list) => current_items.insert(name.clone(), record_list),
                None => None  //TODO we should probably drop this iterator here already
            };

            // store for later
            iterators.insert(name.clone(), the_iter);
        }

        CellIteratorMulti { iterators, current_items}
    }

    fn advance_iter(&mut self, itername: &String) -> Option<(u64, Vec<BusRecord>)>{
        // advance the specified iterator and return the result
        if let Some(iter) = self.iterators.get_mut(itername){
            if let Some(item) = iter.next(){
                println!("Iterator {} yielding item {:?}", itername, item);
                Some(item)
            }
            else{
                println!("Iterator {} empty", itername);
                None
            }
        }
        else{
            // the iterators already got removed
            // technically this shouldnt even happen
            panic!("not supposed to happen");
        }
    }

    fn get_current_min_items(&mut self) -> u64{
        // gets the smallest current CB across the iterators
        let min = 
        self.current_items.values()
                          .map(|(cb, _records)| cb)
                          .reduce(|cb1, cb2|{
                              if cb1<cb2{
                                  cb1
                              }
                              else{
                                  cb2
                              }
                          }).unwrap();
        min.clone() 
    }
}

impl Iterator for CellIteratorMulti {
    // for each file, returns a vector of records for that cell
    type Item = (u64, HashMap<String, Vec<BusRecord>>);

    fn next(&mut self) -> Option<Self::Item> {

        if self.iterators.len() == 0{
            return None
        }

        let current_min_cb = self.get_current_min_items();
        println!("Current min {}", current_min_cb);
        // let current_min_cb = &1;

        // all iterators/names that are at the minimum item will emit
        let names_to_emit: Vec<String>;

        {   //scope to limit the life of current_items
            // not sure if that has anything to do with the problem actually
            //
            // this issue seemd to be that names_to_emit was containing REFERENCES to things (name) in current_items
            // those REFERENCES causes alot of issues down the road
            // if instead we clone `name`, we dont refer to things inside current_items any more
            let current_items = &self.current_items;  // we need a reference here, cant move self.current_item's ownership out of self
            names_to_emit = 
                    current_items.iter().by_ref()
                                .filter(|(_name, (cb, _r))| *cb == current_min_cb)
                                .map(|(name, (_cb, _r))| name.clone())
                                .collect();
        }
        // let items_to_emit: HashMap<&String, &Vec<BusRecord>> = 
        //     self.current_items.iter()
        //                       .filter(|(name, (cb, r))| cb == current_min_cb)
        //                       .map(|(name, (cb, r))| (name, r))
        //                       .collect();
        
        // advance all the iterators that did emit
        // for name in items_to_emit.iter().map(|(name, _rec| name)){
        let mut the_emission: HashMap<String, Vec<BusRecord>> = HashMap::new();

        // lets pop out current items out of the struct: we need to modify it!
        for name in names_to_emit{

            // first pop that item out of current
            let the_item = self.current_items.get(&name).unwrap().clone();  //todo bad clone: to get around the .insert in the next line
            // and add to emission
            the_emission.insert(name.clone(), the_item.1);

            // advance the iterator once more
            // store the result accordingly

            match self.advance_iter(&name){

                Some(cb_rlist) => {
                println!("Advancing {} --> {:?}", name, cb_rlist);

                    self.current_items.insert(name.clone(), cb_rlist); //overwrite the already emitted item
                }, 
                None => {
                    println!("Advancing {} --> EMPTY", name);

                    self.current_items.remove(&name); // clean up the emitted item, remove the iterator itself
                    self.iterators.remove(&name);
                }
            };

            // different plan:
            // we could advance the iterator first: .insert pops out an existing element (thats to be emitted)
        };
        Some((current_min_cb, the_emission))
    }
}

// =================================================================
pub struct CellUmiIteratorMulti {
    pub iterators: HashMap<String, CbUmiIterator>,
    // last_record: BusRecord
    pub current_items: HashMap<String, ((u64, u64), Vec<BusRecord>)> // filename - > (CB, ListOfRecords)
}
impl CellUmiIteratorMulti {

    pub fn new(fnames: &HashMap<String, String>) ->CellUmiIteratorMulti{

        // let iterators: HashMap<String, CellIterator> = fnames.iter()
        //                                              .map(|(name, fname)| (name.clone(), CellIterator::new(fname)))
        //                                              .collect();

        let mut iterators: HashMap<String, CbUmiIterator> = HashMap::new();
        let mut current_items: HashMap<String, ((u64, u64), Vec<BusRecord>)> = HashMap::new();

        for (name, fname) in fnames{
            // create new cell iterator
            let mut the_iter = CbUmiIterator::new(fname);

            // populate first elements from that iterator
            let item = the_iter.next();
            match item{
                Some(record_list) => current_items.insert(name.clone(), record_list),
                None => None  //TODO we should probably drop this iterator here already
            };

            // store for later
            iterators.insert(name.clone(), the_iter);
        }

        CellUmiIteratorMulti { iterators, current_items}
    }

    fn advance_iter(&mut self, itername: &String) -> Option<((u64, u64), Vec<BusRecord>)>{
        // advance the specified iterator and return the result
        if let Some(iter) = self.iterators.get_mut(itername){
            if let Some(item) = iter.next(){
                // println!("Iterator {} yielding item {:?}", itername, item);
                Some(item)
            }
            else{
                // println!("Iterator {} empty", itername);
                None
            }
        }
        else{
            // the iterators already got removed
            // technically this shouldnt even happen
            panic!("not supposed to happen");
        }
    }

    fn get_current_min_items(&mut self) -> (u64, u64){
        // gets the smallest current CB across the iterators
        let min= 
        self.current_items.values()
                          .map(|(cbumi, _records)| cbumi)
                          .min().unwrap();

        *min
    }
}

impl Iterator for CellUmiIteratorMulti {
    // for each file, returns a vector of records for that cell
    type Item = ((u64, u64), HashMap<String, Vec<BusRecord>>);

    fn next(&mut self) -> Option<Self::Item> {

        if self.iterators.len() == 0{
            return None
        }

        let current_min_cb = self.get_current_min_items();
        // println!("Current min {:?}", current_min_cb);
        // let current_min_cb = &1;

        // all iterators/names that are at the minimum item will emit
        let names_to_emit: Vec<String>;

        {   //scope to limit the life of current_items
            // not sure if that has anything to do with the problem actually
            //
            // this issue seemd to be that names_to_emit was containing REFERENCES to things (name) in current_items
            // those REFERENCES causes alot of issues down the road
            // if instead we clone `name`, we dont refer to things inside current_items any more
            let current_items = &self.current_items;  // we need a reference here, cant move self.current_item's ownership out of self
            names_to_emit = 
                    current_items.iter().by_ref()
                                .filter(|(_name, (cbumi, _r))| *cbumi == current_min_cb)
                                .map(|(name, (_cb, _r))| name.clone())
                                .collect();
        }
        // let items_to_emit: HashMap<&String, &Vec<BusRecord>> = 
        //     self.current_items.iter()
        //                       .filter(|(name, (cb, r))| cb == current_min_cb)
        //                       .map(|(name, (cb, r))| (name, r))
        //                       .collect();
        
        // advance all the iterators that did emit
        // for name in items_to_emit.iter().map(|(name, _rec| name)){
        let mut the_emission: HashMap<String, Vec<BusRecord>> = HashMap::new();

        // lets pop out current items out of the struct: we need to modify it!
        for name in names_to_emit{

            // first pop that item out of current
            let the_item = self.current_items.get(&name).unwrap().clone();  //todo bad clone: to get around the .insert in the next line
            // and add to emission
            the_emission.insert(name.clone(), the_item.1);

            // advance the iterator once more
            // store the result accordingly

            match self.advance_iter(&name){

                Some(cb_rlist) => {
                // println!("Advancing {} --> {:?}", name, cb_rlist);

                    self.current_items.insert(name.clone(), cb_rlist); //overwrite the already emitted item
                }, 
                None => {
                    // println!("Advancing {} --> EMPTY", name);

                    self.current_items.remove(&name); // clean up the emitted item, remove the iterator itself
                    self.iterators.remove(&name);
                }
            };

            // different plan:
            // we could advance the iterator first: .insert pops out an existing element (thats to be emitted)
        };
        Some((current_min_cb, the_emission))
    }
}




#[cfg(test)]
mod tests {
    use crate::io::{BusRecord, BusHeader, BusWriter};
    use std::collections::HashMap;
    use super::*;
    #[test]
    fn test_read_write(){

        let r1 =BusRecord{CB: 0, UMI: 21, EC: 0, COUNT: 2, FLAG: 0};
        let r2 =BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0}; 
        let r3 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
        let r4 =BusRecord{CB: 3, UMI: 0, EC: 0, COUNT:  2, FLAG: 0}; 

        let v1 = vec![r1,r2,r3, r4];

        let s1 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let s2 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 12, FLAG: 0};
        let s3 = BusRecord{CB: 2, UMI: 3, EC: 1, COUNT:  2, FLAG: 0}; 
        let v2 = vec![s1,s2,s3];


        let header = BusHeader::new(16, 12, 20);
        let busname1 = "/tmp/test1.bus".to_string();
        let mut writer = BusWriter::new(&busname1, header);
        writer.write_records(&v1);

        let header = BusHeader::new(16, 12, 20);
        let busname2 = "/tmp/test2.bus".to_string();
        let mut writer = BusWriter::new(&busname2, header);
        writer.write_records(&v2);

        let hashmap = HashMap::from([
            ("test1".to_string(), busname1),
            ("test2".to_string(), busname2)
        ]);


        let mut iii = CellIteratorMulti::new(&hashmap); //warning: this triggers the .next() method for both ierators once, consuming the cell 0
        
        println!("========================================");
        println!("Before iter1 {:?}", iii.current_items);
        let e1 = HashMap::from([
            ("test1".to_string(), vec![r1]),
            ("test2".to_string(), vec![s1]),
        ]);
        let o1 = iii.next().unwrap();
        assert_eq!(o1, (0 ,e1));

        println!("========================================");
        println!("After iter1 {:?}", iii.current_items);


        let e2 = HashMap::from([
            ("test1".to_string(), vec![r2, r3]),
            ("test2".to_string(), vec![s2]),
        ]);
        let o2 = iii.next().unwrap();
        assert_eq!(o2, (1 ,e2));

        println!("========================================");
        println!("After iter2 {:?}", iii.current_items);

        let e3 = HashMap::from([
            ("test2".to_string(), vec![s3]),
        ]);
        let o3 = iii.next().expect("something");
        assert_eq!(o3, (2 ,e3));

        println!("========================================");
        println!("After iter3 {:?}", iii.current_items);

        let e4 = HashMap::from([
            ("test1".to_string(), vec![r4]),
        ]);
        let o4 = iii.next().unwrap();
        assert_eq!(o4, (3 ,e4));
    }
}