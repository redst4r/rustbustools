use std::collections::{HashSet, HashMap};
use crate::io::BusRecord;


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

    // indicator, remembering which elements need to be deleted from the first set
    let mut shared_genes:HashSet<String> = HashSet::new();
    
    // initial element, make a copy of that
    let s1 = setlist.next().unwrap();
    for el in s1{
        shared_genes.insert(el.clone());
    }

    // 
    if records.len() == 1{
        return shared_genes
    }

    for current_set in setlist{
        // delete any elemt in shared_genes not present in current_set
        // i..e set intersection

        // we cant delete while iterating, so remember which elements to delete
        let mut to_delete = Vec::new();
        for el in shared_genes.iter(){
            if !current_set.contains(el){
                to_delete.push(el.clone());
            }
        }
        // delete
        for el in to_delete{
            shared_genes.remove(&el);
        }

        // to save some time: if the current intersection is already empty, it'll stay empty
        if shared_genes.len() == 0{
            break
        }
    }
    shared_genes
}


fn find_consistent2(records: &Vec<BusRecord>, ec2gene: &HashMap<u32, HashSet<String>>) ->HashSet<String> {

    // this one has the enumerate instead of doign the fist iteration separately
    // no indicator, just shrinnk the set as we go

    // seems pretty good to me
    let mut shared_genes:HashSet<String> = HashSet::new();

    for (i, r) in records.iter().enumerate(){
        // all genes this record is consitent with
        let g = ec2gene.get(&r.EC).unwrap();
        if i == 0{
            // copy the g hashset
            for el in g{
                shared_genes.insert(el.clone());
            }
        }
        // otherwise delete any elemt in shared_genes not present in g
        // i..e set intersection
        // we cant delete while iterating, so remember which elements to delete
        let mut to_delete = Vec::new();
        for el in shared_genes.iter(){
            if !g.contains(el){
                // shared_genes.remove(el);
                to_delete.push(el.clone());
                // shared_genes_indicator.insert(el.clone(), false);
            }
        }
        // delete
        for el in to_delete{
            shared_genes.remove(&el);
        }

        // to save some time: if the current intersection is already empty, it'll stay empty
        if shared_genes.len() == 0{
            break
        }
    }
    shared_genes
}



fn find_consistent3(records: &Vec<BusRecord>, ec2gene: &HashMap<u32, HashSet<String>>) ->HashSet<String> {

    /*
    set intersection in Rust is a MESS due to ownership etc!!
    Here's the strategy:
    1. copy the first set s1 entirely (the intersection will only contain elements of this set)
    2. go threought the second set s2: if an element in s2 is NOT present in s1, remove it from s1. However, due to rust restrictions, we cant really modify
       Instead, for each element in s1, we keep an indicator, telling us if that element needs to be deleted
    3. delete any element in s1 thats indicated.
    */

    // pretty much same as above find_consistent2, just using the indicator, which will be slower (no set shrink)


    // get the genes per BusRecord
    let mut setlist = records.iter().map(|r| ec2gene.get(&r.EC).unwrap());

    // indicator, remembering which elements need to be deleted from the first set
    let mut shared_genes_indicator:HashMap<String, bool> = HashMap::new();
    
    // initial element
    let s1 = setlist.next().unwrap();
    for el in s1{
        shared_genes_indicator.insert(el.clone(), true);
    }

    // the remaining sets
    for g in setlist{
        // otherwise delete any elemt in shared_genes not present in g
        // i..e set intersection
        // cant iterate and modify!!

        // this is not allowed (chaning while iterating over it)
        // we could do a retain instead
        /*
        for el in shared_genes_indicator.keys(){
            if !g.contains(el){
                // shared_genes_indicator.insert(el.clone(), false);
                shared_genes_indicator.insert(el.clone(), false);
            }
        }
        */
        // shared_genes_indicator.retain(|k,v| k);
    }

    // 3. remove any element that was marked == False
    let mut shared_genes2:HashSet<String> = HashSet::new();
    for (gene, ind) in shared_genes_indicator{
        if ind{
            shared_genes2.insert(gene);
        }
    };
    shared_genes2

}


fn test_consistent(){

}