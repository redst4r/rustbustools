use std::{collections::HashMap, fmt::Debug};

//
//  An iterator that merges mutliple sorted iterators by item
//

pub struct MultiIterator<I,T, K> where I: Iterator<Item=(K, T)> , T: Debug, K: Ord+Eq+Debug+Copy {  // T is the type of elements we iterate (jointly). ie. each iterator to be merged emits T-type items
    pub iterators: HashMap<String, I>,  // todo not sure about this dyn here; compiler asks for it
    pub current_items: HashMap<String, (K,T)> // filename - > (CB, ListOfRecords)
}

impl<I, T, K> MultiIterator<I,T, K> 
where I: Iterator<Item=(K,T)>, T: Debug , K: Ord+Eq+Debug+Copy {

    pub fn new(iterators: HashMap<String, I>) -> Self {

        let mut current_items: HashMap<String, (K, T)> = HashMap::new();

        let mut eee: HashMap<String, I> = HashMap::new(); 

        for (name, mut the_iter) in iterators {

            // populate first elements from that iterator
            let item = the_iter.next();
            match item{
                Some((key, record_list)) => current_items.insert(name.clone(), (key, record_list)),
                None => None  //TODO we should probably drop this iterator here already
            };

            // store for later
            eee.insert(name.clone(), the_iter);
        }

        MultiIterator { iterators: eee, current_items}
    }

    fn advance_iter(&mut self, itername: &String) -> Option<(K, T)>{
        // advance the specified iterator and return the result
        if let Some(iter) = self.iterators.get_mut(itername){
            if let Some((key, item)) = iter.next(){
                // println!("Iterator {} yielding item {:?}", itername, item);
                Some((key, item))
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

    fn get_current_min_items(&mut self) -> K{
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
        *min
    }
}

impl<I, T, K> Iterator for MultiIterator<I, T,K> where I: Iterator<Item=(K,T)> , T: Debug, K: Ord+Eq+Debug+Copy {
    // for each file, returns a vector of records for that cell
    type Item = (K, HashMap<String, T>);

    fn next(&mut self) -> Option<Self::Item> {

        if self.iterators.is_empty() {
            return None
        }

        let current_min_cb = self.get_current_min_items();
        // println!("Current min {:?}", current_min_cb);

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
        
        // advance all the iterators that did emit
        // for name in items_to_emit.iter().map(|(name, _rec| name)){
        let mut the_emission: HashMap<String, T> = HashMap::new();

        // lets pop out current items out of the struct: we need to modify it!
        for name in names_to_emit{

            // first pop that item out of current
            let the_item = self.current_items.remove(&name).unwrap();
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
    use std::collections::HashMap;
    use itertools::izip;
    use crate::{io::{BusRecord, setup_busfile, BusReader}, bus_multi::CellUmiIteratorMulti, iterators::{CbUmiGroupIterator}};
    use super::MultiIterator;


    #[test]
    fn test_basic() {
        let h = HashMap::from([
            ("A".to_owned(), izip!(0..10, 0..10)),
            ("B".to_owned(), izip!(0..10, 10..20) ),
        ]);

        let mut multi = MultiIterator::new(h);
        let x = multi.next().unwrap();
        assert_eq!(
            x, 
            (0, HashMap::from([("A".to_owned(), 0), ("B".to_owned(), 10)])));
        
        let x = multi.next().unwrap();
        assert_eq!(
            x, 
            (1, HashMap::from([("A".to_owned(), 1), ("B".to_owned(), 11)])));

        let x = multi.next().unwrap();
        assert_eq!(
            x, 
            (2, HashMap::from([("A".to_owned(), 2), ("B".to_owned(), 12)])));            
    }
    #[test]
    fn test_compare() {
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let r2 =BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0}; 
        let r3 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
        let r4 =BusRecord{CB: 3, UMI: 0, EC: 0, COUNT:  2, FLAG: 0}; 

        let v1 = vec![r1.clone(), r2.clone(), r3.clone(), r4.clone()];

        let s1 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let s2 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 12, FLAG: 0};
        let s3 = BusRecord{CB: 1, UMI: 3, EC: 1, COUNT: 12, FLAG: 0};
        let s4 = BusRecord{CB: 2, UMI: 3, EC: 1, COUNT:  2, FLAG: 0}; 
        let s5 = BusRecord{CB: 2, UMI: 3, EC: 2, COUNT:  2, FLAG: 0}; 
        let v2 = vec![s1.clone(),s2.clone(),s3.clone(), s4.clone(), s5.clone()];

        // write the records to file
        let (busname1, _dir1) = setup_busfile(&v1);
        let (busname2, _dir2) = setup_busfile(&v2);

        let hashmap = HashMap::from([
            ("test1".to_string(), busname1.to_string()),
            ("test2".to_string(), busname2.to_string()),
        ]);

        // // what we expect to get
        // let expected_pairs = vec![
        //     HashMap::from([
        //         ("test1".to_string(), vec![r1]),
        //         ("test2".to_string(), vec![s1]),
        //     ]),
        //     HashMap::from([("test1".to_string(), vec![r2])]),
        //     HashMap::from([("test2".to_string(), vec![s2])]),
        //     HashMap::from([
        //         ("test1".to_string(), vec![r3]),
        //         ("test2".to_string(), vec![s3]),
        //     ]),
        //     HashMap::from([("test2".to_string(), vec![s4, s5])]),
        //     HashMap::from([("test1".to_string(), vec![r4])]),
        // ];

        // iterate, see if it meets the expectations
        let iii = CellUmiIteratorMulti::new(&hashmap); //warning: this triggers the .next() method for both ierators once, consuming the cell 0

        let hashmap = HashMap::from([
            ("test1".to_string(), BusReader::new(&busname1).groupby_cbumi()),
            ("test2".to_string(), BusReader::new(&busname2).groupby_cbumi()),
        ]);

        let miter = MultiIterator::new(hashmap);
            for ((cbumi1, q),(cbumi2, w)) in izip!(iii, miter){
                assert_eq!(cbumi1, cbumi2);
                assert_eq!(q, w);
            }

    }
}