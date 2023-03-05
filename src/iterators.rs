use std::collections::{HashSet, HashMap};

use crate::{io::{BusIteratorBuffered, BusRecord, setup_busfile}, consistent_genes::{Ec2GeneMapper, CUGset, groubygene}, disjoint::DisjointSubsets, utils::vec2set};

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
            // anything left in the basic iterator (if not, just emit whatever is stored in self.last_element)
            if let Some(new_record) = self.busiter.next(){

                // the newly seen record
                let (new_cb, new_umi) = (new_record.CB, new_record.UMI);

                // take ownership of self.last_record, which we're goign to emit now (since we have a new item)
                // replace by the new item
                let last_record = std::mem::replace(&mut self.last_record, Some(new_record)).unwrap();

                let (current_cb, current_umi) = (last_record.CB, last_record.UMI);

                busrecords.push(last_record); // the stored element from the previous iteration

                // now we just need to decide if we want to emit, or continue growing
                if new_cb > current_cb || (new_cb == current_cb &&  new_umi > current_umi){  
                    // we ran into a new CB/UMI and its records
                    // println!("\tyielding {:?}", (current_cb, &busrecords));

                    return Some(((current_cb, current_umi), busrecords));
                }
                else if new_cb == current_cb && new_umi == current_umi {
                    // nothing happens, just keep growing busrecords
                }
                else{  // the new cb is smaller then the current state: this is a bug due to an UNOSORTED busfile
                    panic!("Unsorted busfile: {}/{} -> {}/{}", current_cb, current_umi, new_cb, new_umi)
                }
            }
            else{ // emit whatever is left in last element
                // we get the ownership and replace with None (which in the next iterator will trigger the end of the entire iterator)
                // to mark the end of iteration and all items emitted, set last_item to None
                let last_record = std::mem::replace(&mut self.last_record, None);

                // get the last element and emit
                if let Some(r) = last_record {  
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    let current_cb = r.CB;
                    let current_umi = r.UMI;
                    busrecords.push(r);
                    return Some(((current_cb, current_umi), busrecords));    
                }   
                else{
                    return None
                }          
            }
        }
    }
}

//=================================================================================
pub struct CellIterator {
    pub(crate) busiter: BusIteratorBuffered,
    pub(crate) last_record: Option<BusRecord>,  //option needed to mark the final element of the iteration
    // buffersize: usize
}

impl CellIterator {
    pub fn new(fname: &str) ->CellIterator{
        let mut busiter = BusIteratorBuffered::new(fname);
        let last_record = busiter.next(); //initilize with the first record in the file
        CellIterator {busiter, last_record}
        // CellIterator {busiter, last_record, buffersize:1}
    }
}

impl Iterator for CellIterator {
    type Item = (u64, Vec<BusRecord>);

    fn next(&mut self) -> Option<Self::Item> {

        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted
        // let mut busrecords: Vec<BusRecord> = Vec::with_capacity(self.buffersize); // storing the result to be emitted
    
        loop {
            if let Some(new_record) = self.busiter.next(){
                // the newly seen record
                let new_cb = new_record.CB;

                // take ownership of self.last_record, which we're goign to emit now (since we have a new item)
                // replace by the new item
                // note that before this line, self.last_record CANNOT logically be None, hence the unwrap. 
                // If it is somethings wrong with my logic
                let last_record = std::mem::replace(&mut self.last_record, Some(new_record)).unwrap();

                let current_cb = last_record.CB;

                busrecords.push(last_record); // the stored element from the previous iteration

                // now we just need to decide if we want to emit, or continue growing
                match new_cb.cmp(&current_cb){
                    std::cmp::Ordering::Equal => { }, //nothing happens, just keep growing busrecords
                    std::cmp::Ordering::Greater => {return Some((current_cb, busrecords));},
                    std::cmp::Ordering::Less => {panic!("Unsorted busfile: {} -> {}", current_cb, new_cb)}
                }
            }
            else{
                // get the last element and emit

                // we get the ownership and replace with None (which in the next iterator will trigger the end of the entire iterator)
                // to mark the end of iteration and all items emitted, set last_item to None
                let last_record = std::mem::replace(&mut self.last_record, None);

                if let Some(r) = last_record {  
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    let current_cb = r.CB;
                    busrecords.push(r);
                    return Some((current_cb, busrecords));    
                }   
                else{
                    return None
                }          
            }
        }
    }
}


/* adaptor to be able to group things by gene 
CbUmiIterator.iter().group_by_gene()

*/
pub struct GroupbyGene<I> {
    iter: I,
    ecmapper: Ec2GeneMapper
}


impl<I> Iterator for GroupbyGene<I>
where I: Iterator<Item = Vec<BusRecord>>
{
    type Item = Vec<CUGset>;
    fn next(&mut self) -> Option<Self::Item> {

        self.iter.next()
            .map(|v| 
                groubygene(v, &self.ecmapper)
            )
    }
}

impl<I> GroupbyGene<I> {
    pub fn new(iter: I, ecmapper: Ec2GeneMapper) -> Self {
        Self { iter, ecmapper }
    }
}

/* Stuff to get the iterator chaining working! Just a wrapper around GroupbyGene::new() */
pub trait GroupbyGeneIterator<T>: Iterator<Item = T> + Sized {
    fn group_by_gene(self, ecmapper: Ec2GeneMapper) -> GroupbyGene<Self> {
        GroupbyGene::new(self, ecmapper)
    }
}
impl<T, I: Iterator<Item = T>> GroupbyGeneIterator<T> for I {}




#[cfg(test)]
mod tests{
    use std::collections::{HashSet, HashMap};

    use crate::consistent_genes::Ec2GeneMapper;
    use crate::io::{BusRecord, setup_busfile};
    use crate::iterators::{CbUmiIterator, CellIterator, GroupbyGeneIterator};
    use crate::utils::vec2set;

    #[test]
    fn test_cb_iter_last_element1(){   
        // make sure everything is emitted, even if the last element is its own group
        let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};  

        let records = vec![r1.clone(),r2.clone(),r3.clone()];
        let (busname, _dir) = setup_busfile(&records);

        let n: Vec<_> = CellIterator::new(&busname).map(|(_cb, records)| records).collect();
        // println!("{:?}", n);
        // assert_eq!(n, vec![vec![r1, r2], vec![r3]]);
        assert_eq!(n.len(), 2);

        let rlist = &n[1];
        assert_eq!(rlist.len(), 1);
    }

    #[test]
    fn test_cb_iter_last_element2(){   
        // make sure everything is emitted, even if the last element has mutiple elements
        let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};  
        let r4 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};  
        let records = vec![r1.clone(),r2.clone(),r3.clone(), r4.clone()];
        let (busname, _dir) = setup_busfile(&records);

        let n: Vec<_> = CellIterator::new(&busname).collect();
        // println!("{:?}", n);

        assert_eq!(n.len(), 2);

        let (_cb, rlist) = &n[1];
        assert_eq!(rlist.len(), 2);
    }


    #[test]
    fn test_cb_iter(){   
        let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r4 = BusRecord{CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r5 = BusRecord{CB: 2, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r6 = BusRecord{CB: 3, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
    
        let records = vec![r1.clone(),r2.clone(),r3.clone(),r4.clone(),r5.clone(), r6.clone()];
        // let records = vec![r1,r2,r3,r4,r5, r6].to_vec();
    
        let (busname, _dir) = setup_busfile(&records);
    

        let cb_iter = CellIterator::new(&busname);
        let n: Vec<(u64, Vec<BusRecord>)> = cb_iter.collect();
    
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

    // #[test]
    fn test_cb_iter_speed(){
        use std::time::Instant;
        let foldername = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
        let n=100000;
        let biter2= CellIterator::new(foldername);
    
        let now = Instant::now();
        let s2: Vec<_> = biter2.take(n).map(|(_a, records)|records).collect();
        let elapsed_time = now.elapsed();
        println!("Running CellIterator({}) took {} seconds.", n, elapsed_time.as_secs());

    }

    #[test]
    fn test_cbumi_iter(){   
        let r1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r4 = BusRecord{CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r5 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0};
        let r6 = BusRecord{CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

        let records = vec![r1.clone(),r2.clone(),r3.clone(),r4.clone(),r5.clone(), r6.clone()];

        let (busname, _dir) = setup_busfile(&records);

        let cb_iter = CbUmiIterator::new(&busname);
        // let n: Vec<Vec<BusRecord>> = cb_iter.map(|(_cb, rec)| rec).collect();
        let n: Vec<((u64, u64), Vec<BusRecord>)> = cb_iter.collect();
        // println!("{:?}", n);

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

        // let busname = "/tmp/test_panic_on_unsorted.bus";
        let (busname, _dir) = setup_busfile(&records);

        let cb_iter = CellIterator::new(&busname);
        let _n: Vec<(u64, Vec<BusRecord>)> = cb_iter.collect();
    }

    #[test]
    fn test_groupby_genes(){
        let ec0: HashSet<String> = vec2set(vec!["A".to_string()]);
        let ec1: HashSet<String> = vec2set(vec!["B".to_string()]);
        let ec2: HashSet<String> = vec2set(vec!["A".to_string(), "B".to_string()]);
        let ec3: HashSet<String> = vec2set(vec!["C".to_string(), "D".to_string()]);
    
        let ec_dict: HashMap<u32, HashSet<String>> = HashMap::from([
            (0, ec0.clone()), 
            (1, ec1.clone()), 
            (2, ec2.clone()), 
            (3, ec3.clone()), 
            ]);
    
        let es = Ec2GeneMapper::new(ec_dict);
    
        // first sample: two records, with consistent gene A
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let r2 =BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 2, FLAG: 0};
    
         // second sample: two records, with consistent gene A the other consistent with gene B
         let s1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 3, FLAG: 0}; // A
         let s2 = BusRecord{CB: 0, UMI: 2, EC: 1, COUNT: 4, FLAG: 0}; //B
    
         let records = vec![r1.clone(),r2.clone(),s1.clone(),s2.clone()];
    
         let (busname, _dir) = setup_busfile(&records);
    
         let cb_iter = CbUmiIterator::new(&busname);
    
         let results: Vec<_> = cb_iter
            .map(|(cbumi, r)| r).group_by_gene(es).collect();
    
        assert_eq!(results.len(), 2);
        assert_eq!(results[0].len(), 1);
        assert_eq!(results[1].len(), 2);
    
        println!("{:?}", results)
    }    
 }