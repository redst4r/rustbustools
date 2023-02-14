use std::{collections::{HashMap, HashSet}, hash::Hash};

pub struct DisjointSubsets<T>{
    pub disjoint_sets: HashMap<String, HashSet<T>>
}

pub const SEPARATOR: &str ="@SEP@";

impl <T:Hash+Eq> DisjointSubsets<T>{
    pub fn new() -> Self{
        DisjointSubsets{
            disjoint_sets: HashMap::new()
        }
    }

    pub fn add(&mut self, setname:String, aset:HashSet<T>){

        if self.disjoint_sets.contains_key(&setname){
            panic!("inserting into existing setname")
        }
        // check for any existing set that shares a member with aset
        let mut candidates: Vec<String> = Vec::new();
        for (name, the_set) in &self.disjoint_sets{
            if set_overlap(&aset, the_set){
                candidates.push(name.clone());
            }
        }

        // aset links everything in the candidates into a single cell
        // create that new set!
        let mut newset:HashSet<T> = HashSet::new();
        let mut newname = String::new();
        for c in candidates{
            let mut s = self.disjoint_sets.remove(&c).unwrap(); //poing the old set out
            for el in s.drain(){
                newset.insert(el);
            }
            newname.push_str(&c);
            // newname.push('_');
            newname.push_str(SEPARATOR);
        }

        // add the set itself
        for el in aset{
            newset.insert(el);
        }
        newname.push_str(&setname);

        // insert the newly made set
        self.disjoint_sets.insert(newname, newset);
    }
}

fn set_overlap<T: Hash+Eq>(aset: &HashSet<T>, bset: &HashSet<T>) -> bool{
    for el in aset{
        if bset.contains(el){
            return true
        }
    }
    false
}

#[test]
fn testing_disjoint(){

    let a: HashSet<_> = vec!["A", "B","C"].into_iter().collect();
    let b: HashSet<_> = vec!["D", "E"].into_iter().collect();
    let c: HashSet<_> = vec!["B", "D"].into_iter().collect();
    let d: HashSet<_> = vec!["Z"].into_iter().collect();

    let mut ds = DisjointSubsets::new();
    ds.add("set1".to_string(), a);
    ds.add("set2".to_string(), b);

    println!("{:?}", ds.disjoint_sets);

    ds.add("set3".to_string(), c);
    println!("{:?}", ds.disjoint_sets);

    ds.add("set4".to_string(), d);
    println!("{:?}", ds.disjoint_sets);

}