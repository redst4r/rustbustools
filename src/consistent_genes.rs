use std::collections::{HashSet, HashMap};
use crate::io::BusRecord;
use std::hash::Hash;


// #[inline(never)]
fn update_intersection_via_retain<T:Hash+Eq>(inter:  &mut HashSet<T>, newset: &HashSet<T>){        
    // delete any elemt in shared_genes not present in current_set
    // i..e set intersection
    // we cant delete while iterating, so remember which elements to delete
    inter.retain(|item| newset.contains(item));
}


pub fn find_consistent(records: &Vec<BusRecord>, ec2gene: &HashMap<u32, HashSet<String>>) ->HashSet<String> {

    /*
    set intersection in Rust is a MESS due to ownership etc!!
    Here's the strategy:
    1. copy the first set s1 entirely (the intersection will only contain elements of this set)
    2. go threought the second set s2: if an element in s2 is NOT present in s1, remove it from s1. 
    */

    // this one does shrink the set, and uses the iterator

    // get the genes per BusRecord
    let mut setlist = records.iter().map(|r| ec2gene.get(&r.EC).unwrap());

    let s1 = setlist.next().unwrap();

    let mut shared_genes;
    if false{
        shared_genes = HashSet::new(); // could spec capacity
    
        // initial element, make a copy of that
        for el in s1{
            shared_genes.insert(el.clone());
        }
    }
    // pretty much just a clone of s1
    else{
        shared_genes = s1.clone();
    }

    // save some time if its only one record
    if records.len() == 1{
        return shared_genes
    }

    for current_set in setlist{
        // delete any elemt in shared_genes not present in current_set
        // i..e set intersection
        update_intersection_via_retain(&mut shared_genes, current_set);

        // to save some time: if the current intersection is already empty, it'll stay empty
        if shared_genes.len() == 0{
            break
        }
    }
    shared_genes
}


#[test]
fn test_consistent(){

    let ec0: HashSet<String> = vec!["A".to_string()].into_iter().collect();
    let ec1: HashSet<String> = vec!["B".to_string()].into_iter().collect();
    let ec2: HashSet<String> = vec!["A".to_string(), "B".to_string()].into_iter().collect();
    let ec3: HashSet<String> = vec!["C".to_string(), "D".to_string()].into_iter().collect();

    let ec_dict: HashMap<u32, HashSet<String>> = HashMap::from([
        (0, ec0), 
        (1, ec1), 
        (2, ec2), 
        (3, ec3), 
        ]);

    // not find_consistent doesnt check the records indeed come from the same CB/UMI

    // single read, consistent with A
    let r1 =BusRecord{CB: 0, UMI: 21, EC: 0, COUNT: 2, FLAG: 0};
    let res1: Vec<String> = find_consistent(&vec![r1], &ec_dict).into_iter().collect();
    assert_eq!(vec!["A"], res1);

    // single read, consistent with A and B, hence multimapped
    let r2 =BusRecord{CB: 0, UMI: 21, EC: 2, COUNT: 2, FLAG: 0};
    let res2: Vec<String> = find_consistent(&vec![r2], &ec_dict).into_iter().collect();
    // assert_eq!(vec!["A", "B"], res2);  WRNING the order matters here, the function might return the vec ordered differently
    assert_eq!(res2.len(), 2);


    // two reads, consistent with A
    let r3 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
    let r4 =BusRecord{CB: 3, UMI: 0, EC: 2, COUNT:  2, FLAG: 0}; 
    let res3: Vec<String> = find_consistent(&vec![r3, r4], &ec_dict).into_iter().collect();
    assert_eq!(vec!["A"], res3);

    // two reads, inconsistent with A, B
    let r5 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
    let r6 =BusRecord{CB: 3, UMI: 0, EC: 1, COUNT:  2, FLAG: 0}; 
    let res4: Vec<String> = find_consistent(&vec![r5, r6], &ec_dict).into_iter().collect();
    assert_eq!(res4.len(), 0);

    // three reads,  A, B, (A,B)
    // inconsintent at the end
    let r7 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
    let r8 =BusRecord{CB: 3, UMI: 0, EC: 1, COUNT:  2, FLAG: 0}; 
    let r9 =BusRecord{CB: 3, UMI: 0, EC: 2, COUNT:  2, FLAG: 0}; 
    let res5: Vec<String> = find_consistent(&vec![r7, r8, r9], &ec_dict).into_iter().collect();
    assert_eq!(res5.len(), 0);

}