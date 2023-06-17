//! Module for the Intersector struct
//!
//!
//! The [Intersector] represents a list of key/values, such that
//! a new item i with key k is added to a key/value pair if the key k
//! intersects with the key in the k/v pair
//!
//! e.g.
//! ```
//!  Action           Intersector state
//! add {'A'}:1 ->   {'A'}:1
//! add {'B'}:2 ->   {'A'}:1, {'B'}:2  / since keys dont overlap
//! add {'A'}:3 ->   {'A'}:[1,3], {'B'}:2  // new element gets merged into the first group "A"
//! add {'A', 'B'}: 4 ->   {'A'}:[1,3, 4], {'B'}:2 //note, its added to the first key found
//! ```
//!
use std::{collections::HashSet, hash::Hash};

use itertools::izip;

fn set_overlap<T: Hash + Eq>(aset: &HashSet<T>, bset: &HashSet<T>) -> bool {
    for el in aset {
        if bset.contains(el) {
            return true;
        }
    }
    false
}

fn set_intersection<T: Hash + Eq + Clone>(a: &HashSet<T>, b: &HashSet<T>) -> HashSet<T> {
    // returns a new set, containing the intersection (all elements are cloned)
    a.intersection(b).cloned().collect()
}

/// special "HashMap" with keys: Sets of type T and values of `Vec<S>`
/// if keys overlap, their values will be merged
pub struct Intersector<T, S> {
    // pub the_sets: HashMap<HashSet<T>, Vec<S> >
    // emulating a dict of Set -> list of items
    pub keys: Vec<HashSet<T>>,
    pub items: Vec<Vec<S>>,
}
impl<T: Hash + Eq + Clone, S> Intersector<T, S> {
    pub fn new() -> Self {
        Intersector { keys: Vec::new(), items: Vec::new() }
    }

    /// Iterate over the key/value pairs
    pub fn iterate_items(&self) -> impl Iterator<Item = (&HashSet<T>, &Vec<S>)> {
        izip!(&self.keys, &self.items)
    }

    /// return the index in the data structure where set overlaps with an element
    fn find_overlap_ix(&self, set: &HashSet<T>) -> Option<usize> {
        for (i, k) in self.keys.iter().enumerate() {
            if set_overlap(k, set) {
                return Some(i);
            }
        }
        None
    }
    /// add a single key/value pair
    pub fn add(&mut self, set: HashSet<T>, item: S) {
        match self.find_overlap_ix(&set) {
            // it the item overlaps with anything (the first entry is returned)
            // update the overlaping key (intersection) and value (adding the current ites)
            Some(i) => {
                let old_key = self.keys.remove(i);
                let mut old_items = self.items.remove(i);

                let update_key = set_intersection(&old_key, &set);
                old_items.push(item);

                self.keys.push(update_key);
                self.items.push(old_items);
            }
            // if no overlap exists, just insert a new k/v pair
            None => {
                self.keys.push(set);
                self.items.push(vec![item]);
            }
        }
    }
}

#[cfg(test)]
mod testing {
    use crate::{disjoint::Intersector, utils::vec2set};

    #[test]
    fn test_inter() {
        let mut it: Intersector<String, usize> = Intersector::new();

        let k1 = vec2set(vec!["A".to_string()]);
        let k2 = vec2set(vec!["B".to_string()]);
        it.add(k1, 1);
        it.add(k2, 2);
        println!("{:?}", it.items);
        assert_eq!(it.keys.len(), 2);

        let mut it: Intersector<String, usize> = Intersector::new();
        let k1 = vec2set(vec!["A".to_string()]);
        let k2 = vec2set(vec!["A".to_string(), "B".to_string()]);
        it.add(k1, 1);
        it.add(k2, 2);
        println!("{:?}", it.items);

        assert_eq!(it.keys.len(), 1);
        assert_eq!(it.keys, vec![vec2set(vec!["A".to_string()]),]);

        let mut it: Intersector<String, usize> = Intersector::new();
        let k1 = vec2set(vec!["A".to_string(), "B".to_string()]);
        let k2 = vec2set(vec!["A".to_string()]);
        let k3 = vec2set(vec!["C".to_string()]);
        it.add(k1, 1);
        it.add(k2, 2);
        it.add(k3, 3);
        println!("{:?}", it.items);

        assert_eq!(it.keys.len(), 2);
        assert_eq!(
            it.keys,
            vec![
                vec2set(vec!["A".to_string()]),
                vec2set(vec!["C".to_string()]),
            ]
        );
    }
}
