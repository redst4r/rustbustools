use std::{
    collections::{HashMap, HashSet},
    hash::Hash,
};

pub struct DisjointSubsets<T> {
    pub disjoint_sets: HashMap<String, HashSet<T>>,
}

pub const SEPARATOR: &str = "@SEP@";

impl<T: Hash + Eq> DisjointSubsets<T> {
    pub fn new() -> Self {
        DisjointSubsets {
            disjoint_sets: HashMap::new(),
        }
    }

    pub fn add(&mut self, setname: String, aset: HashSet<T>) {
        if self.disjoint_sets.contains_key(&setname) {
            panic!("inserting into existing setname")
        }
        // check for any existing set that shares a member with aset
        let mut candidates: Vec<String> = Vec::new();
        for (name, the_set) in &self.disjoint_sets {
            if set_overlap(&aset, the_set) {
                candidates.push(name.clone());
            }
        }

        // aset links everything in the candidates into a single cell
        // create that new set!
        let mut newset: HashSet<T> = HashSet::new();
        let mut newname = String::new();
        for c in candidates {
            let mut s = self.disjoint_sets.remove(&c).unwrap(); //poing the old set out
            for el in s.drain() {
                newset.insert(el);
            }
            newname.push_str(&c);
            // newname.push('_');
            newname.push_str(SEPARATOR);
        }

        // add the set itself
        for el in aset {
            newset.insert(el);
        }
        newname.push_str(&setname);

        // insert the newly made set
        self.disjoint_sets.insert(newname, newset);
    }

    pub fn get_disjoint_set_ids(&self) -> Vec<Vec<String>> {
        // return the list of disjoint sets, for each set, report it by its element's IDs
        let mut setlist: Vec<Vec<String>> = Vec::with_capacity(self.disjoint_sets.len());
        for id_str in self.disjoint_sets.keys() {
            let s: Vec<String> = id_str.split(SEPARATOR).map(|x| x.to_string()).collect();
            setlist.push(s);
        }
        setlist
    }

    pub fn get_disjoint_set_ids_and_set(&self) -> Vec<(Vec<String>, &HashSet<T>)> {
        // return the list of disjoint sets, for each set, report it by its element's IDs and the set
        let mut setlist: Vec<(Vec<String>, &HashSet<T>)> =
            Vec::with_capacity(self.disjoint_sets.len());
            
        for id_str in self.disjoint_sets.keys() {
            let s: Vec<String> = id_str.split(SEPARATOR).map(|x| x.to_string()).collect();
            let theset = self.disjoint_sets.get(id_str).unwrap();
            setlist.push((s, theset));
        }
        setlist
    }
}

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

/*
 Represents a list of key/values, such that
 a new item i with key k is added to a key/value pair if the key k
 intersects with the key in the k/v pair

 e.g.
 add {'A'}:1 ->   {'A'}:1
 add {'B'}:2 ->   {'A'}:1, {'B'}:2
 add {'A'}:3 ->   {'A'}:[1,3], {'B'}:2
 add {'A', 'B'}: 4 ->   {'A'}:[1,3, 4], {'B'}:2 //note, its added to the first key found

*/
pub struct Intersector<T, S> {
    // pub the_sets: HashMap<HashSet<T>, Vec<S> >
    // emulating a dict of Set -> list of items
    pub keys: Vec<HashSet<T>>,
    pub items: Vec<Vec<S>>,
}
impl<T: Hash + Eq + Clone, S> Intersector<T, S> {
    pub fn new() -> Self {
        Intersector {
            keys: Vec::new(),
            items: Vec::new(),
        }
    }

    fn find_overlap_ix(&self, set: &HashSet<T>) -> Option<usize> {
        // return the index in the data structure where set overlaps with an element
        for (i, k) in self.keys.iter().enumerate() {
            if set_overlap(k, set) {
                return Some(i);
            }
        }
        None
    }
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

    // pub fn drain(&mut self){

    // }
}

#[cfg(test)]
mod testing {
    use crate::{
        disjoint::{DisjointSubsets, Intersector},
        utils::vec2set,
    };
    use std::collections::{HashMap, HashSet};

    #[test]
    fn testing_disjoint() {
        let a = vec2set(vec!["A", "B", "C"]);
        let b = vec2set(vec!["D", "E"]);
        let c = vec2set(vec!["B", "D"]);
        let d = vec2set(vec!["Z"]);

        let mut ds = DisjointSubsets::new();
        ds.add("set1".to_string(), a);
        ds.add("set2".to_string(), b);

        let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        let set1 = vec2set(vec!["vec2set(A", "B", "C"]);
        expected.insert("set1".to_string(), set1);
        let set2 = vec2set(vec!["D", "E"]);
        expected.insert("set2".to_string(), set2);

        assert_eq!(ds.disjoint_sets, expected);

        // println!("{:?}", ds.disjoint_sets);

        ds.add("set3".to_string(), c);

        let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        let set1 = vec2set(vec!["A", "B", "C", "D", "E"]);
        expected.insert("set1@SEP@set2@SEP@set3".to_string(), set1);
        assert_eq!(ds.disjoint_sets, expected);
        // println!("{:?}", ds.disjoint_sets);

        ds.add("set4".to_string(), d);

        let mut expected: HashMap<String, HashSet<&str>> = HashMap::new();
        let set1 = vec2set(vec!["A", "B", "C", "D", "E"]);
        expected.insert("set1@SEP@set2@SEP@set3".to_string(), set1);
        let set2 = vec2set(vec!["Z"]);
        expected.insert("set4".to_string(), set2);
        // println!("{:?}", ds.disjoint_sets);
        assert_eq!(ds.disjoint_sets, expected);
    }

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
