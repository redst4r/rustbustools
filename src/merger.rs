//! An iterator that merges mutliple sorted iterators by item
//!
use std::{collections::HashMap, fmt::Debug, mem::replace};

/// Iterator the merges mulitple sorted iterators by item
///
/// I: Iterator yielding keys K and values T
pub struct MultiIteratorSlow<I, T, K>
where
    I: Iterator<Item = (K, T)>,
    T: Debug,
    K: Ord + Eq + Debug + Copy,
{
    // T is the type of elements we iterate (jointly). ie. each iterator to be merged emits T-type items
    iterators: HashMap<String, I>, // todo not sure about this dyn here; compiler asks for it
    current_items: HashMap<String, (K, T)>, // filename - > (CB, ListOfRecords)
}

impl<I, T, K> MultiIteratorSlow<I, T, K>
where
    I: Iterator<Item = (K, T)>,
    T: Debug,
    K: Ord + Eq + Debug + Copy,
{
    /// construct a key-merged iterator from a dictionary/hashmap of iterators
    pub fn new(iterators: HashMap<String, I>) -> Self {
        let mut current_items: HashMap<String, (K, T)> = HashMap::new();

        let mut eee: HashMap<String, I> = HashMap::new();

        for (name, mut the_iter) in iterators {
            // populate first elements from that iterator
            let item = the_iter.next();
            match item {
                Some((key, record_list)) => current_items.insert(name.clone(), (key, record_list)),
                None => None, //TODO we should probably drop this iterator here already
            };

            // store for later
            eee.insert(name.clone(), the_iter);
        }

        MultiIteratorSlow { iterators: eee, current_items }
    }

    fn advance_iter(&mut self, itername: &String) -> Option<(K, T)> {
        // advance the specified iterator and return the result
        if let Some(iter) = self.iterators.get_mut(itername) {
            if let Some((key, item)) = iter.next() {
                // println!("Iterator {} yielding item {:?}", itername, item);
                Some((key, item))
            } else {
                // println!("Iterator {} empty", itername);
                None
            }
        } else {
            // the iterators already got removed
            // technically this shouldnt even happen
            panic!("not supposed to happen");
        }
    }

    fn get_current_min_items(&mut self) -> K {
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

impl<I, T, K> Iterator for MultiIteratorSlow<I, T, K>
where
    I: Iterator<Item = (K, T)>,
    T: Debug,
    K: Ord + Eq + Debug + Copy,
{
    // for each file, returns a vector of records for that cell
    type Item = (K, HashMap<String, T>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.iterators.is_empty() {
            return None;
        }

        let current_min_cb = self.get_current_min_items();
        // println!("Current min {:?}", current_min_cb);

        // all iterators/names that are at the minimum item will emit
        let names_to_emit: Vec<String> =
        {
            //scope to limit the life of current_items
            // not sure if that has anything to do with the problem actually
            //
            // this issue seemd to be that names_to_emit was containing REFERENCES to things (name) in current_items
            // those REFERENCES causes alot of issues down the road
            // if instead we clone `name`, we dont refer to things inside current_items any more
            let current_items = &self.current_items; // we need a reference here, cant move self.current_item's ownership out of self
            current_items
                .iter()
                .by_ref()
                .filter(|(_name, (cb, _r))| *cb == current_min_cb)
                .map(|(name, (_cb, _r))| name.clone())
                .collect()
        };
        let mut the_emission: HashMap<String, T> = HashMap::new();

        if false {
            // advance all the iterators that did emit
            // for name in items_to_emit.iter().map(|(name, _rec| name)){

            // lets pop out current items out of the struct: we need to modify it!
            for name in names_to_emit {
                // first pop that item out of current
                let the_item = self.current_items.remove(&name).unwrap();
                // and add to emission
                the_emission.insert(name.clone(), the_item.1);

                // advance the iterator once more
                // store the result accordingly
                match self.advance_iter(&name) {
                    Some(cb_rlist) => {
                        // println!("Advancing {} --> {:?}", name, cb_rlist);
                        self.current_items.insert(name.clone(), cb_rlist); //overwrite the already emitted item
                    }
                    None => {
                        // println!("Advancing {} --> EMPTY", name);
                        self.current_items.remove(&name); // clean up the emitted item, remove the iterator itself
                        self.iterators.remove(&name);
                    }
                };
                // different plan:
                // we could advance the iterator first: .insert pops out an existing element (thats to be emitted)
            }
        } else {
            // pull the iterators, and either insert/pop (if the iterator returned something) 
            // or just pop (in case we're done with that iterator)
            //  the old iterms out of self.current_items
            // 
            for name in names_to_emit {
                // advance iterator
                match self.advance_iter(&name) {
                    Some(cb_rlist) => {
                        // update the new item, emit the old

                        // let dest  = self.current_items.get_mut(&name).unwrap();
                        // let to_emit = mem::replace(dest, cb_rlist);

                        let to_emit = self.current_items.insert(name.clone(), cb_rlist).expect("element must exist"); //overwrite the already emitted item
                        the_emission.insert(name.clone(), to_emit.1);
                    },
                    // iterator done
                    None => {
                        let to_emit = self.current_items.remove(&name).expect("element must exist"); // clean up the emitted item, remove the iterator itself
                        the_emission.insert(name.clone(), to_emit.1);
                        self.iterators.remove(&name);
                    },
                }
            }
        }
        Some((current_min_cb, the_emission))
    }
}


/// A faster version of the [`MultiIteratorSlow`], working around the performance issue of using hashmaps 
/// (each emission generates hashmaps which is expensive).
///
/// Rather, we keep three arrays of size n_iterators and do linear search to find the correct one to poll.
/// 
/// ## Generics:
/// - `I`: an iterator of Item=(K,T)
/// - `K`: the key to group items (e.g. CB,UMI)
/// - `T`: the type of item a single iterator emits ([`crate::io::BusRecord`])
/// 
pub struct MultiIterator<I, T, K>
where
    I: Iterator<Item = (K, T)>,
    T: Debug,
    K: Ord + Eq + Debug + Copy,
{
    // T is the type of elements we iterate (jointly). ie. each iterator to be merged emits T-type items
    names: Vec<String>,
    iterators: Vec<I>,
    current_items: Vec<Option<(K, T)>>, // filename - > (CB, ListOfRecords)
}
impl<I, T, K> MultiIterator<I, T, K>
where
    I: Iterator<Item = (K, T)>,
    T: Debug,
    K: Ord + Eq + Debug + Copy,
{
    /// construct a key-merged iterator from a dictionary/hashmap of iterators
    pub fn new(iterators: HashMap<String, I>) -> Self {

        let mut names = Vec::new();
        let mut iii = Vec::new();
        let mut current_items = Vec::new();

        for (name, mut the_iter) in iterators {
            // populate first elements from that iterator
            let item = the_iter.next();
            current_items.push(item);
            iii.push(the_iter);
            names.push(name);
        }
        MultiIterator { names, iterators: iii, current_items }
    }

    /// advance the specified iterator and return the result
    fn advance_iter(&mut self, ix: usize) -> Option<(K, T)> {
        // TODO looks like this is just a wrap around `self.iterators[ix].next()`
        if let Some((key, item)) = self.iterators[ix].next() {
            Some((key, item))
        } else {
            None
        }
    }

    /// gets the smallest key (K) across the iterators
    fn get_current_min_items(&mut self) -> K {
        let min = self.current_items.iter()
            .flatten()
            .map(|(cb, _records)| cb)
            .reduce(|cb1, cb2| if cb1 < cb2 { cb1 } else { cb2 })
            .unwrap();
        *min
    }

    /// remove an iterator from the structure, cleaning up the arrays.
    /// Typically done when the iterator is empty (this will panic if called on a nonempy iterator)
    fn remove_iterator(&mut self, name: &str) {
        let ix = self.names.iter().position(|x| x==name).expect("name needs to be present in iterators");
        let _ii = self.iterators.remove(ix);
        let _name = self.names.remove(ix);
        let curr = self.current_items.remove(ix);
        assert!(curr.is_none());
    }

}

impl<I, T, K> Iterator for MultiIterator<I, T, K>
where
    I: Iterator<Item = (K, T)>,
    T: Debug,
    K: Ord + Eq + Debug + Copy,
{
    // for each file, returns a vector of records for that cell
    type Item = (K, HashMap<String, T>);

    fn next(&mut self) -> Option<Self::Item> {
        if self.iterators.is_empty() {
            return None;
        }

        let current_min_cb = self.get_current_min_items();
        // println!("Current min {:?}", current_min_cb);
        // all iterators/names that are at the minimum item will emit

        let mut ix_to_emit = Vec::with_capacity(self.names.len());
        for (ix, item) in  self.current_items.iter().enumerate(){
            if let Some((key, _ii)) = item {
                if *key == current_min_cb {
                    ix_to_emit.push(ix);
                }
            }
        }

        // TODO hashbrown with ahash??
        let mut the_emission: HashMap<String, T> = HashMap::new();

        // advance all the iterators that did emit
        // 1. advance the respective iterator (a new element to be stored)
        // 2. pop out the preexisting element in the iterator
        // 3. emit the popped element, insert the new
        let mut names_for_cleanup = Vec::new();
        for ix in ix_to_emit {

            // advance the iterator, swap othe old with the new item
            let name_to_emit = self.names[ix].clone();
            let item_to_emit = match self.advance_iter(ix) {
                Some(new_item) => {
                    // let item_to_emit = replace(&mut self.current_items[ix], Some(new_item));
                    // the_emission.insert(name_to_emit, item_to_emit.expect("should never be none").1);
                    replace(&mut self.current_items[ix], Some(new_item))
                },
                None => {
                    // cleanup that iterator; NOTE: cant do it here directly, the next loop iteration relies on the indices being the same
                    names_for_cleanup.push(name_to_emit.clone());
                    // let item_to_emit = self.current_items[ix].take();  // does the same as above: pop and insert non
                    // the_emission.insert(name_to_emit.clone(), item_to_emit.expect("should never be none").1);
                    self.current_items[ix].take() // does the same as above: pop and insert non
                },
            };
            the_emission.insert(name_to_emit, item_to_emit.expect("should never be none").1);
        }

        // clean up exhausted iterators
        for name in names_for_cleanup {
            self.remove_iterator(&name);
        }

        Some((current_min_cb, the_emission))
    }
}

#[cfg(test)]
mod tests {
    use super::{MultiIteratorSlow, MultiIterator};
    use crate::{
        bus_multi::{CellIteratorMulti, CellUmiIteratorMulti},
        io::{BusReader, BusRecord},
        iterators::{CbUmiGroupIterator, CellGroupIterator},
    };
    use itertools::izip;
    use std::collections::HashMap;

    #[test]
    fn test_basic() {
        let h = HashMap::from([
            ("A".to_owned(), izip!(0..10, 0..10)),
            ("B".to_owned(), izip!(0..10, 10..20)),
        ]);

        let mut multi = MultiIteratorSlow::new(h);
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


        // iterate, see if it meets the expectations
        let hashmap = HashMap::from([
            ("test1".to_string(), v1.clone().into_iter().groupby_cbumi()),
            ("test2".to_string(), v2.clone().into_iter().groupby_cbumi()),
        ]);        
        let iii = CellUmiIteratorMulti::new(hashmap); //warning: this triggers the .next() method for both ierators once, consuming the cell 0
        let hashmap = HashMap::from([
            ("test1".to_string(), v1.clone().into_iter().groupby_cbumi()),
            ("test2".to_string(), v2.clone().into_iter().groupby_cbumi()),
        ]);

        let miter = MultiIteratorSlow::new(hashmap);
            for ((cbumi1, q),(cbumi2, w)) in izip!(iii, miter){
                assert_eq!(cbumi1, cbumi2);
                assert_eq!(q, w);
            }

        let hashmap = HashMap::from([
            ("test1".to_string(), v1.clone().into_iter().groupby_cbumi()),
            ("test2".to_string(), v2.clone().into_iter().groupby_cbumi()),
        ]);
        let iii = CellUmiIteratorMulti::new(hashmap); //warning: this triggers the .next() method for both ierators once, consuming the cell 0

        let hashmap = HashMap::from([
            ("test1".to_string(), v1.clone().into_iter().groupby_cbumi()),
            ("test2".to_string(), v2.clone().into_iter().groupby_cbumi()),
        ]);
        let miter = MultiIterator::new(hashmap);
        for ((cbumi1, q),(cbumi2, w)) in izip!(iii, miter){
            assert_eq!(cbumi1, cbumi2);
            assert_eq!(q, w);
        }
    }


    #[test]
    fn test_real_merger() {
        let input1 = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        let input2 = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";

        let hashmap = HashMap::from([
            ("test1".to_string(), BusReader::new(&input1).groupby_cbumi()),
            ("test2".to_string(), BusReader::new(&input2).groupby_cbumi()),
        ]);
        let m = MultiIteratorSlow::new(hashmap);
        let mut counter = 0;
        let before = std::time::Instant::now();

        for (_cbumi, w) in m{
            counter += w.len();
        }
        let elapsed_time = before.elapsed();
        println!("{counter}, {} ms", elapsed_time.as_millis());
    }

    #[test]
    fn test_real_merger_fast() {
        let input1 = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        let input2 = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";

        let hashmap = HashMap::from([
            ("test1".to_string(), BusReader::new(&input1).groupby_cbumi()),
            ("test2".to_string(), BusReader::new(&input2).groupby_cbumi()),
        ]);

        let m = MultiIterator::new(hashmap);
        let mut counter = 0;

        let before = std::time::Instant::now();
        for (_cbumi, w) in m{
            counter += w.len();
        }
        let elapsed_time = before.elapsed();
        println!("{counter}, {} ms", elapsed_time.as_millis());
    }

    // TODO write test where two iterators are removed in the same iteration (changes indexing)
    #[test]
    fn remove_multiple() {
        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0 };
        let v1 = vec![r1.clone()];

        let s1 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
        let s2 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 12, FLAG: 0 };
        let s3 = BusRecord { CB: 1, UMI: 3, EC: 1, COUNT: 12, FLAG: 0 };
        let s4 = BusRecord { CB: 2, UMI: 3, EC: 1, COUNT: 2, FLAG: 0 };
        let s5 = BusRecord { CB: 2, UMI: 3, EC: 2, COUNT: 2, FLAG: 0 };
        let v2 = vec![s1.clone(), s2.clone(), s3.clone(), s4.clone(), s5.clone()];

        let r3 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0 };
        let v3 = vec![r3.clone()];


        let hashmap = HashMap::from([
            ("test1".to_string(), v1.clone().into_iter().groupby_cb()),
            ("test2".to_string(), v2.clone().into_iter().groupby_cb()),
            ("test3".to_string(), v3.clone().into_iter().groupby_cb()),
        ]);        
        let hashmap2 = HashMap::from([
            ("test1".to_string(), v1.clone().into_iter().groupby_cb()),
            ("test2".to_string(), v2.clone().into_iter().groupby_cb()),
            ("test3".to_string(), v3.clone().into_iter().groupby_cb()),
        ]);        

        let iii = CellIteratorMulti::new(hashmap); //warning: this triggers the .next() method for both ierators once, consuming the cell 0

        let miter = MultiIterator::new(hashmap2);
        for ((cbumi1, q),(cbumi2, w)) in izip!(iii, miter){
            assert_eq!(cbumi1, cbumi2);
            assert_eq!(q, w);
        }
    }

}
