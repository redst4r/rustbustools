//! Module for more advanced iterators over busrecords.
//! Allows to iterate over 
//! * cells (CB)
//! * mRNAs (CB+UMI)
//! * CB+UMI+gene: in case some CB/UMI collision
//! 
//! To iterate over a **sorted** busfile, grouping all records by CB:
//! ```rust, no_run
//! # use bustools::io::BusReader;
//! use bustools::iterators::CellGroupIterator; //need to bring that trait into scope
//! 
//! let breader = BusReader::new("/path/to/some.bus");
//! for (cb, vector_of_records) in breader.groupby_cb() {
//!     // Example: the number of records in that cell
//!     let n_molecules: usize = vector_of_records.len();
//!     
//! }
//! ```
//! 
//! To iterate over a `sorted` busfile, grouping all records by CB+UMI:
//! ```rust, no_run
//! # use bustools::io::BusReader; 
//! use bustools::iterators::CbUmiGroupIterator; //need to bring that trait into scope
//! 
//! let breader = BusReader::new("/path/to/some.bus");
//! for ((cb, umi), vector_of_records) in breader.groupby_cbumi() {
//!     // Example: the number of reads of that molecule (CB/UMI)
//!     let n_reads: u32 = vector_of_records.iter().map(|r| r.COUNT).sum();
//! }
//! ```
use crate::{
    consistent_genes::{groubygene, CUGset, Ec2GeneMapper},
    io::{BusRecord, CUGIterator},
};

/// groups an iterator over BusRecords by cell+umi
pub struct CbUmiGroup<I: CUGIterator> {
    iter: I,
    last_record: Option<BusRecord>, //option needed to mark the final element of the iteration
}

impl<I> Iterator for CbUmiGroup<I>
where
    I: CUGIterator,
{
    type Item = ((u64, u64), Vec<BusRecord>);
    fn next(&mut self) -> Option<Self::Item> {
        // let mut busrecords: Vec<BusRecord> = Vec::with_capacity(10); // storing the result to be emitted; capacity is small, since we dont expect many records for the same CB/UMI
        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted; capacity is small, since we dont expect many records for the same CB/UMI

        loop {
            // anything left in the basic iterator (if not, just emit whatever is stored in self.last_element)
            if let Some(new_record) = self.iter.next() {
                // the newly seen record
                let (new_cb, new_umi) = (new_record.CB, new_record.UMI);

                // take ownership of self.last_record, which we're goign to emit now (since we have a new item)
                // replace by the new item
                let last_record =
                    std::mem::replace(&mut self.last_record, Some(new_record)).unwrap();

                let (current_cb, current_umi) = (last_record.CB, last_record.UMI);

                busrecords.push(last_record); // the stored element from the previous iteration

                // now we just need to decide if we want to emit, or continue growing
                if new_cb > current_cb || (new_cb == current_cb && new_umi > current_umi) {
                    // we ran into a new CB/UMI and its records
                    // println!("\tyielding {:?}", (current_cb, &busrecords));

                    return Some(((current_cb, current_umi), busrecords));
                } else if new_cb == current_cb && new_umi == current_umi {
                    // nothing happens, just keep growing busrecords
                } else {
                    // the new cb is smaller then the current state: this is a bug due to an UNOSORTED busfile
                    panic!(
                        "Unsorted busfile: {}/{} -> {}/{}",
                        current_cb, current_umi, new_cb, new_umi
                    )
                }
            } else {
                // emit whatever is left in last element
                // we get the ownership and replace with None (which in the next iterator will trigger the end of the entire iterator)
                // to mark the end of iteration and all items emitted, set last_item to None
                // let last_record = std::mem::replace(&mut self.last_record, None);
                let last_record = self.last_record.take(); // swaps the current value for None; clippy suggestion to std::mem::replace
                                                           // get the last element and emit
                if let Some(r) = last_record {
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    let current_cb = r.CB;
                    let current_umi = r.UMI;
                    busrecords.push(r);
                    return Some(((current_cb, current_umi), busrecords));
                } else {
                    return None;
                }
            }
        }
    }
}

impl<I> CbUmiGroup<I>
where
    I: CUGIterator,
{
    pub fn new(mut iter: I) -> Self {
        let last_record = iter.next(); //initilize with the first record in the file
        Self { iter, last_record }
    }
}

/// gets iterator chaining working! Just a wrapper around CbUmiGroupIterator::new()
pub trait CbUmiGroupIterator: CUGIterator + Sized {
    fn groupby_cbumi(self) -> CbUmiGroup<Self> {
        CbUmiGroup::new(self)
    }
}
// implements the .groupby_cbumi() synthax for any `CUGIterator`
impl<I: CUGIterator> CbUmiGroupIterator for I {}

//=================================================================================
/// groups an iterator over BusRecords by cell
pub struct CellGroup<I: CUGIterator> {
    iter: I,
    last_record: Option<BusRecord>, //option needed to mark the final element of the iteration
}

impl<I> Iterator for CellGroup<I>
where
    I: CUGIterator,
{
    type Item = (u64, Vec<BusRecord>);
    fn next(&mut self) -> Option<Self::Item> {
        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted
                                                         // let mut busrecords: Vec<BusRecord> = Vec::with_capacity(self.buffersize); // storing the result to be emitted

        loop {
            if let Some(new_record) = self.iter.next() {
                // the newly seen record
                let new_cb = new_record.CB;

                // take ownership of self.last_record, which we're goign to emit now (since we have a new item)
                // replace by the new item
                // note that before this line, self.last_record CANNOT logically be None, hence the unwrap.
                // If it is somethings wrong with my logic
                let last_record =
                    std::mem::replace(&mut self.last_record, Some(new_record)).unwrap();

                let current_cb = last_record.CB;

                busrecords.push(last_record); // the stored element from the previous iteration

                // now we just need to decide if we want to emit, or continue growing
                match new_cb.cmp(&current_cb) {
                    std::cmp::Ordering::Equal => {} //nothing happens, just keep growing busrecords
                    std::cmp::Ordering::Greater => {
                        return Some((current_cb, busrecords));
                    }
                    std::cmp::Ordering::Less => {
                        panic!("Unsorted busfile: {} -> {}", current_cb, new_cb)
                    }
                }
            } else {
                // get the last element and emit

                // we get the ownership and replace with None (which in the next iterator will trigger the end of the entire iterator)
                // to mark the end of iteration and all items emitted, set last_item to None
                // let last_record = std::mem::replace(&mut self.last_record, None);
                let last_record = self.last_record.take(); // swaps the current value for None; clippy suggestion to std::mem::replace

                if let Some(r) = last_record {
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    let current_cb = r.CB;
                    busrecords.push(r);
                    return Some((current_cb, busrecords));
                } else {
                    return None;
                }
            }
        }
    }
}

impl<I> CellGroup<I>
where
    I: CUGIterator,
{
    pub fn new(mut iter: I) -> Self {
        let last_record = iter.next(); //initilize with the first record in the file
        Self { iter, last_record }
    }
}

/// gets iterator chaining working! Just a wrapper around CellGroupIterator::new()
pub trait CellGroupIterator: CUGIterator + Sized {
    fn groupby_cb(self) -> CellGroup<Self> {
        CellGroup::new(self)
    }
}
impl<I: CUGIterator> CellGroupIterator for I {}

/// enables to group things by gene
/// `CbUmiIterator.iter().group_by_gene()`
/// works for any iterator yielding `Vec<BusRecord>`
pub struct GroupbyGene<I> {
    iter: I,
    ecmapper: Ec2GeneMapper,
}

impl<I> Iterator for GroupbyGene<I>
where
    I: Iterator<Item = Vec<BusRecord>>,
{
    type Item = Vec<CUGset>;
    fn next(&mut self) -> Option<Self::Item> {
        self.iter.next().map(|v| groubygene(v, &self.ecmapper))
    }
}

impl<I> GroupbyGene<I> {
    pub fn new(iter: I, ecmapper: Ec2GeneMapper) -> Self {
        Self { iter, ecmapper }
    }
}

/// gets iterator chaining working! Just a wrapper around GroupbyGene::new() 
pub trait GroupbyGeneIterator<T>: Iterator<Item = T> + Sized {
    fn group_by_gene(self, ecmapper: Ec2GeneMapper) -> GroupbyGene<Self> {
        GroupbyGene::new(self, ecmapper)
    }
}
impl<T, I: Iterator<Item = T>> GroupbyGeneIterator<T> for I {}

#[cfg(test)]
mod tests {
    use std::collections::{HashMap, HashSet};

    use crate::consistent_genes::{Ec2GeneMapper, Genename, EC};
    use crate::io::BusReader;
    use crate::io::{setup_busfile, BusRecord};
    use crate::iterators::{CbUmiGroupIterator, CellGroup, CellGroupIterator};
    use crate::utils::vec2set;

    #[test]
    fn test_cb_iter_last_element1() {
        // make sure everything is emitted, even if the last element is its own group
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let r3 = BusRecord { CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };

        let records = vec![r1.clone(), r2.clone(), r3.clone()];
        // let n: Vec<_> = records.into_iter().groupby_cb().map(|(_cb, records)| records).collect();
        let (busname, _dir) = setup_busfile(&records);

        let b = BusReader::new(&busname);
        let n: Vec<_> = b.groupby_cb().map(|(_cb, records)| records).collect();

        // println!("{:?}", n);
        // assert_eq!(n, vec![vec![r1, r2], vec![r3]]);
        assert_eq!(n.len(), 2);

        let rlist = &n[1];
        assert_eq!(rlist.len(), 1);

        // another wayto initialize,  no chaining
        let (busname, _dir) = setup_busfile(&records);
        let b = BusReader::new(&busname);
        let n: Vec<_> = CellGroup::new(b).map(|(_cb, records)| records).collect();
        assert_eq!(n.len(), 2);

        let rlist = &n[1];
        assert_eq!(rlist.len(), 1);
    }

    #[test]
    fn test_cb_iter_last_element2() {
        // make sure everything is emitted, even if the last element has mutiple elements
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let r3 = BusRecord { CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r4 = BusRecord { CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let records = vec![r1.clone(), r2.clone(), r3.clone(), r4.clone()];
        let (busname, _dir) = setup_busfile(&records);

        let b = BusReader::new(&busname);
        let n: Vec<_> = b.groupby_cb().collect();

        // println!("{:?}", n);

        assert_eq!(n.len(), 2);

        let (_cb, rlist) = &n[1];
        assert_eq!(rlist.len(), 2);
    }

    #[test]
    fn test_cb_iter() {
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let r3 = BusRecord { CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r4 = BusRecord { CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
        let r5 = BusRecord { CB: 2, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let r6 = BusRecord { CB: 3, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };

        let records = vec![
            r1.clone(),
            r2.clone(),
            r3.clone(),
            r4.clone(),
            r5.clone(),
            r6.clone(),
        ];

        // let n: Vec<(u64, Vec<BusRecord>)> = records.into_iter().groupby_cb().collect();

        let (busname, _dir) = setup_busfile(&records);

        let b = BusReader::new(&busname);
        let n: Vec<(u64, Vec<BusRecord>)> = b.groupby_cb().collect();

        assert_eq!(n.len(), 4);
        // println!("{:?}", n);
        // println!("First");
        let c1 = &n[0];
        assert_eq!(*c1, (0, vec![r1, r2]));

        // println!("Second");
        let c2 = &n[1];
        assert_eq!(*c2, (1, vec![r3]));

        // println!("Third");
        let c3 = &n[2];
        assert_eq!(*c3, (2, vec![r4, r5]));

        let c4 = &n[3];
        assert_eq!(*c4, (3, vec![r6]));

        // assert_eq!(n, vec![vec![r1,r2], vec![r3], vec![r4,r5]])
    }

    #[test]
    fn test_cbumi_iter() {
        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
        let r3 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r4 = BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
        let r5 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0 };
        let r6 = BusRecord { CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };

        let records = vec![
            r1.clone(),
            r2.clone(),
            r3.clone(),
            r4.clone(),
            r5.clone(),
            r6.clone(),
        ];
        // let n: Vec<((u64, u64), Vec<BusRecord>)> = records.into_iter().groupby_cbumi().collect();

        let (busname, _dir) = setup_busfile(&records);

        let b = BusReader::new(&busname);
        let cb_iter = b.groupby_cbumi();
        // let n: Vec<Vec<BusRecord>> = cb_iter.map(|(_cb, rec)| rec).collect();
        let n: Vec<((u64, u64), Vec<BusRecord>)> = cb_iter.collect();
        // println!("{:?}", n);

        assert_eq!(n.len(), 5);
        // println!("{:?}", n);
        // println!("First");
        let c1 = &n[0];
        assert_eq!(*c1, ((0, 1), vec![r1, r2]));

        // println!("Second");
        let c2 = &n[1];
        assert_eq!(*c2, ((0, 2), vec![r3]));

        // println!("Third");
        let c3 = &n[2];
        assert_eq!(*c3, ((1, 1), vec![r4]));

        let c4 = &n[3];
        assert_eq!(*c4, ((1, 2), vec![r5]));

        let c5 = &n[4];
        assert_eq!(*c5, ((2, 1), vec![r6]));
        // assert_eq!(n, vec![vec![r1,r2], vec![r3], vec![r4,r5]])
    }

    #[test]
    #[should_panic(expected = "Unsorted busfile: 2 -> 0")]
    fn test_panic_on_unsorted() {
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let r3 = BusRecord { CB: 2, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r4 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };

        let records = vec![r1, r2, r3, r4];

        // let busname = "/tmp/test_panic_on_unsorted.bus";
        let (busname, _dir) = setup_busfile(&records);
        let b = BusReader::new(&busname);
        b.groupby_cb().count();
    }

    #[test]
    #[should_panic(expected = "Unsorted busfile: 2/2 -> 0/1")]
    fn test_panic_on_unsorted_cbumi() {
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let r3 = BusRecord { CB: 2, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r4 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };

        let records = vec![r1, r2, r3, r4];

        let (busname, _dir) = setup_busfile(&records);
        let b = BusReader::new(&busname);
        b.groupby_cbumi().count();
    }

    use crate::iterators::GroupbyGeneIterator;
    #[test]
    fn test_groupby_genes() {
        let ec0: HashSet<Genename> = vec2set(vec![Genename("A".to_string())]);
        let ec1: HashSet<Genename> = vec2set(vec![Genename("B".to_string())]);
        let ec2: HashSet<Genename> =
            vec2set(vec![Genename("A".to_string()), Genename("B".to_string())]);
        let ec3: HashSet<Genename> =
            vec2set(vec![Genename("C".to_string()), Genename("D".to_string())]);

        let ec_dict: HashMap<EC, HashSet<Genename>> = HashMap::from([
            (EC(0), ec0.clone()),
            (EC(1), ec1.clone()),
            (EC(2), ec2.clone()),
            (EC(3), ec3.clone()),
        ]);

        let es = Ec2GeneMapper::new(ec_dict);

        // first sample: two records, with consistent gene A
        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 1, EC: 2, COUNT: 2, FLAG: 0 };

        // second sample: two records, with consistent gene A the other consistent with gene B
        let s1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 3, FLAG: 0 }; // A
        let s2 = BusRecord { CB: 0, UMI: 2, EC: 1, COUNT: 4, FLAG: 0 }; //B

        let records = vec![r1.clone(), r2.clone(), s1.clone(), s2.clone()];

        let (busname, _dir) = setup_busfile(&records);

        let b = BusReader::new(&busname);

        let cb_iter = b.groupby_cbumi();

        let results: Vec<_> = cb_iter.map(|(_cbumi, r)| r).group_by_gene(es).collect();

        assert_eq!(results.len(), 2);
        assert_eq!(results[0].len(), 1);
        assert_eq!(results[1].len(), 2);

        println!("{:?}", results)
    }
}

// use itertools::{GroupBy, Itertools};
//  pub struct CbUmiGroup_GROUP {
//     grouped_iter: Box<dyn Iterator<Item=BusRecord>>,
// }

// impl Iterator for CbUmiGroup_GROUP
// {
//     type Item = ((u64, u64), Vec<BusRecord>);
//     fn next(&mut self) -> Option<Self::Item> {
//         self.grouped_iter.next()
//     }
// }

// impl CbUmiGroup_GROUP
// {
//     pub fn new(mut iter: dyn Iterator<Item=BusRecord> + Sized) -> Self {
//         let grouped_iter = iter.group_by(|r| (r.CB, *r));
//         Self { grouped_iter }
//     }
// }

// pub trait CbUmiGroup_GROUPIterator: CUGIterator + Sized {
//     fn groupby_cbumi(self) -> CbUmiGroup_GROUP {
//         CbUmiGroup_GROUP::new(self)
//     }
// }
// impl<I: CUGIterator> CbUmiGroup_GROUPIterator for I {}
