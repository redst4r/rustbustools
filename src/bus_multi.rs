//! A module that allows iteration of multiple busfiles simulatniously
//! 
//! One can iterate e.g. two busfiles, going through it cell-by-cell, i.e. 
//! if the same cell is present in both files, we get their BusRecords in a single emission
//! 
//! # Example
//! ```rust, no_run
//! # use std::collections::HashMap;
//! # use rustbustools::bus_multi::CellUmiIteratorMulti;
//! // two busfiles, named
//! let hashmap = HashMap::from([
//!     ("test1".to_string(), "/tmp/some1.bus".to_string()),
//!     ("test2".to_string(), "/tmp/some2.bus".to_string())]);
//! let iii = CellUmiIteratorMulti::new(&hashmap); //warning: this triggers the .next() method for both ierators once, consuming the cell 0
//! 
//! for (cb, emission_dict) in iii{
//!     // get the records wrt to cb
//!     // if the file doesnt contain that cell 
//!     // emission_dict wont contain an entry
//!     let records1 = emission_dict.get("test1");
//!     let records2 = emission_dict.get("test2");
//! }
//! ```
//! # Structure
//! Taking multiple iterators over busfiles, we merge them into a new iterator
//! which emits the cell (cell/umi) and a dictionary.
//! The dict contains the file identifiers a keys, and values are a Vec of BusRecords with that cell (cell/umi) in the respective file
//! There's two main iterators in here
//! - CellIteratorMulti  -> iterate, grouping by cell
//! - CellUmiIteratorMulti -> iterate, grouping bt cell/umi
//! 
use crate::io::{BusReader, BusRecord};
use crate::iterators::{CbUmiGroup, CellGroup};
use std::collections::HashMap;

pub struct CellIteratorMulti {
    pub iterators: HashMap<String, CellGroup<BusReader>>,
    pub current_items: HashMap<String, (u64, Vec<BusRecord>)>, // filename - > (CB, ListOfRecords)
}

impl CellIteratorMulti {
    pub fn new(fnames: &HashMap<String, String>) -> CellIteratorMulti {
        let mut iterators: HashMap<String, CellGroup<BusReader>> = HashMap::new();
        let mut current_items: HashMap<String, (u64, Vec<BusRecord>)> = HashMap::new();

        for (name, fname) in fnames {
            // create new cell iterator
            let mut the_iter = CellGroup::new(BusReader::new(fname));

            // populate first elements from that iterator
            let item = the_iter.next();
            match item {
                Some(record_list) => current_items.insert(name.clone(), record_list),
                None => None, //TODO we should probably drop this iterator here already
            };

            // store for later
            iterators.insert(name.clone(), the_iter);
        }

        CellIteratorMulti { iterators, current_items }
    }

    fn advance_iter(&mut self, itername: &String) -> Option<(u64, Vec<BusRecord>)> {
        // advance the specified iterator and return the result
        if let Some(iter) = self.iterators.get_mut(itername) {
            iter.next()
        } else {
            // the iterators already got removed
            // technically this shouldnt even happen
            panic!("not supposed to happen");
        }
    }

    fn get_current_min_items(&mut self) -> u64 {
        // gets the smallest current CB across the iterators
        let min = self
            .current_items
            .values()
            .map(|(cb, _records)| cb)
            .reduce(|cb1, cb2| if cb1 < cb2 { cb1 } else { cb2 })
            .unwrap();
        *min
    }
}

impl Iterator for CellIteratorMulti {
    // for each file, returns a vector of records for that cell
    type Item = (u64, HashMap<String, Vec<BusRecord>>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.iterators.is_empty() {
            return None;
        }

        let current_min_cb = self.get_current_min_items();
        // println!("Current min {}", current_min_cb);

        // all iterators/names that are at the minimum item will emit
        let names_to_emit: Vec<String> = self
            .current_items
            .iter()
            .filter_map(|(sname, (cb, _r))| {
                if *cb == current_min_cb {
                    Some(sname.clone())
                } else {
                    None
                }
            })
            .collect();

        // advance all the iterators that did emit
        // for name in items_to_emit.iter().map(|(name, _rec| name)){
        let mut the_emission: HashMap<String, Vec<BusRecord>> = HashMap::new();

        // lets pop out current items out of the struct: we need to modify it!
        for name in names_to_emit {
            // first pop that item out of current
            let (_cb, the_item) = self.current_items.remove(&name).unwrap(); //todo bad clone: to get around the .insert in the next line

            // and add to emission
            the_emission.insert(name.clone(), the_item);

            // advance the iterator once more
            // store the result accordingly

            match self.advance_iter(&name) {
                Some(cb_rlist) => {
                    // println!("Advancing {} --> {:?}", name, cb_rlist);
                    self.current_items.insert(name.clone(), cb_rlist); //overwrite the already emitted item
                }
                None => {
                    // println!("Advancing {} --> EMPTY", name);
                    self.iterators.remove(&name);
                }
            };

            // different plan:
            // we could advance the iterator first: .insert pops out an existing element (thats to be emitted)
        }
        Some((current_min_cb, the_emission))
    }
}

// =================================================================
pub struct CellUmiIteratorMulti {
    pub iterators: HashMap<String, CbUmiGroup<BusReader>>,
    // last_record: BusRecord
    pub current_items: HashMap<String, ((u64, u64), Vec<BusRecord>)>, // filename - > (CB, ListOfRecords)
}
impl CellUmiIteratorMulti {
    pub fn new(fnames: &HashMap<String, String>) -> CellUmiIteratorMulti {
        let mut iterators: HashMap<String, CbUmiGroup<BusReader>> = HashMap::new();
        let mut current_items: HashMap<String, ((u64, u64), Vec<BusRecord>)> = HashMap::new();

        for (name, fname) in fnames {
            // create new cell iterator
            let mut the_iter = CbUmiGroup::new(BusReader::new(fname));

            // populate first elements from that iterator
            let item = the_iter.next();
            match item {
                Some(record_list) => current_items.insert(name.clone(), record_list),
                None => None, //TODO we should probably drop this iterator here already
            };

            // store for later
            iterators.insert(name.clone(), the_iter);
        }
        CellUmiIteratorMulti { iterators, current_items }
    }

    fn advance_iter(&mut self, itername: &String) -> Option<((u64, u64), Vec<BusRecord>)> {
        // advance the specified iterator and return the result
        if let Some(iter) = self.iterators.get_mut(itername) {
            iter.next()
        } else {
            // the iterators already got removed
            // technically this shouldnt even happen
            panic!("not supposed to happen");
        }
    }

    fn get_current_min_items(&mut self) -> (u64, u64) {
        // gets the smallest current CB across the iterators
        let min = self
            .current_items
            .values()
            .map(|(cbumi, _records)| cbumi)
            .min()
            .unwrap();
        //TODO is this min over tuples really working?
        *min
    }
}

impl Iterator for CellUmiIteratorMulti {
    // for each file, returns a vector of records for that cell
    type Item = ((u64, u64), HashMap<String, Vec<BusRecord>>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.iterators.is_empty() {
            return None;
        }

        let current_min_cb = self.get_current_min_items();
        // println!("Current min {:?}", current_min_cb);

        // all iterators/names that are at the minimum item will emit
        let names_to_emit: Vec<String> = self
            .current_items
            .iter()
            .filter_map(|(sname, (cbumi, _r))| {
                if *cbumi == current_min_cb {
                    Some(sname.clone())
                } else {
                    None
                }
            })
            .collect();

        // advance all the iterators that did emit
        // for name in items_to_emit.iter().map(|(name, _rec| name)){
        let mut the_emission: HashMap<String, Vec<BusRecord>> = HashMap::new();

        // lets pop out current items out of the struct: we need to modify it!
        for name in names_to_emit {
            // first pop that item out of current
            let (_cbumi, the_item) = self.current_items.remove(&name).unwrap(); //todo bad clone: to get around the .insert in the next line
            // and add to emission
            the_emission.insert(name.clone(), the_item);

            // advance the iterator once more
            // store the result accordingly

            match self.advance_iter(&name) {
                Some(cb_rlist) => {
                    // println!("Advancing {} --> {:?}", name, cb_rlist);
                    self.current_items.insert(name.clone(), cb_rlist); //overwrite the already emitted item
                }
                None => {
                    // println!("Advancing {} --> EMPTY", name);
                    self.iterators.remove(&name);
                }
            }
        }
        Some((current_min_cb, the_emission))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::io::{setup_busfile, BusRecord};
    use std::collections::HashMap;

    #[test]
    fn test_read_write() {
        /*
        iterate two busfiles in parallel across cells
        it should yield the following:
        r1, s1
        r2,r3, s2
        s3
        r4
         */
        use std::iter::zip;
        let r1 = BusRecord { CB: 0, UMI: 21, EC: 0, COUNT: 2, FLAG: 0 };
        let r2 = BusRecord { CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r3 = BusRecord { CB: 1, UMI: 3, EC: 0, COUNT: 2, FLAG: 0 };
        let r4 = BusRecord { CB: 3, UMI: 0, EC: 0, COUNT: 2, FLAG: 0 };

        let v1 = vec![r1.clone(), r2.clone(), r3.clone(), r4.clone()];

        let s1 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
        let s2 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 12, FLAG: 0 };
        let s3 = BusRecord { CB: 2, UMI: 3, EC: 1, COUNT: 2, FLAG: 0 };
        let v2 = vec![s1.clone(), s2.clone(), s3.clone()];

        // write the records to file
        let (busname1, _dir1) = setup_busfile(&v1);
        let (busname2, _dir2) = setup_busfile(&v2);

        // create_busfile(busname1, v1);
        // create_busfile(busname2, v2);

        let hashmap = HashMap::from([
            ("test1".to_string(), busname1.to_string()),
            ("test2".to_string(), busname2.to_string()),
        ]);

        // what we expect to get
        let expected_pairs = vec![
            HashMap::from([
                ("test1".to_string(), vec![r1]),
                ("test2".to_string(), vec![s1]),
            ]),
            HashMap::from([
                ("test1".to_string(), vec![r2, r3]),
                ("test2".to_string(), vec![s2]),
            ]),
            HashMap::from([("test2".to_string(), vec![s3])]),
            HashMap::from([("test1".to_string(), vec![r4])]),
        ];

        // iterate, see if it meets the expectations
        let iii = CellIteratorMulti::new(&hashmap); //warning: this triggers the .next() method for both ierators once, consuming the cell 0

        for (r_expect, (_cb, r_obs)) in zip(expected_pairs, iii) {
            assert_eq!(r_expect, r_obs);
        }
    }

    #[test]
    fn test_cbumi_multi() {
        use std::iter::zip;
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

        let hashmap = HashMap::from([
            ("test1".to_string(), busname1.to_string()),
            ("test2".to_string(), busname2.to_string()),
        ]);

        // what we expect to get
        let expected_pairs = vec![
            HashMap::from([
                ("test1".to_string(), vec![r1]),
                ("test2".to_string(), vec![s1]),
            ]),
            HashMap::from([("test1".to_string(), vec![r2])]),
            HashMap::from([("test2".to_string(), vec![s2])]),
            HashMap::from([
                ("test1".to_string(), vec![r3]),
                ("test2".to_string(), vec![s3]),
            ]),
            HashMap::from([("test2".to_string(), vec![s4, s5])]),
            HashMap::from([("test1".to_string(), vec![r4])]),
        ];

        // iterate, see if it meets the expectations
        let iii = CellUmiIteratorMulti::new(&hashmap); //warning: this triggers the .next() method for both ierators once, consuming the cell 0

        for (r_expect, (_cb, r_obs)) in zip(expected_pairs, iii) {
            assert_eq!(r_expect, r_obs);
        }
    }
}
