use std::collections::HashMap;
use crate::io::BusRecord;

//
//  An iterator that merges mutliple sorted iterators by item
//

struct Testing2<T,S> 
where S: Ord, T: Iterator<Item=S>
{
    iter: HashMap<usize, T>,
}
impl<T,S> Testing2<T,S> where S: Ord, T: Iterator<Item=S> {
    pub fn new(iterators: HashMap<usize, T>) -> Self {
        Testing2 { iter: iterators }
    }
}
impl<T, S> Iterator for Testing2<T, S>
where S: Ord, T: Iterator<Item=S> {
    type Item = HashMap<usize, S>;

    fn next(&mut self) -> Option<Self::Item> {
        let mut h = HashMap::new();
        for (name, the_iter) in self.iter.iter_mut() {
            if let Some(el) = the_iter.next() {
                h.insert(*name, el);
            }
        }
        Some(h)
    }
}


pub struct MultiIterator<T> {  // T is the type of elements we iterate (jointly). ie. each iterator to be merged emits T-type items
    pub iterators: HashMap<String, dyn Iterator<Item=T>>,  // todo not sure about this dyn here; compiler asks for it
    pub current_items: HashMap<String, T> // filename - > (CB, ListOfRecords)
}
impl<T> MultiIterator<T> {

    pub fn new(iterators: HashMap<String, dyn Iterator<Item=T>>) ->MultiIterator<T>{

        let mut current_items: HashMap<String, T> = HashMap::new();

        for (name, the_iter) in iterators{

            // populate first elements from that iterator
            let item = the_iter.next();
            match item{
                Some(record_list) => current_items.insert(name.clone(), record_list),
                None => None  //TODO we should probably drop this iterator here already
            };

            // store for later
            iterators.insert(name.clone(), the_iter);
        }

        MultiIterator { iterators, current_items}
    }

    fn advance_iter(&mut self, itername: &String) -> Option<T>{
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

    fn get_current_min_items(&mut self) -> T{
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

impl<T> Iterator for MultiIterator<T> {
    // for each file, returns a vector of records for that cell
    type Item = HashMap<String, T>;

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


#[cfg(test)]
mod tests {}