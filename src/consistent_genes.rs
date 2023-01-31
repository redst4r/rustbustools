use std::collections::{HashSet, HashMap};
use crate::io::BusRecord;
use std::hash::Hash;

/*
just some wrappers for Strings and ints that we use for genes and equivalence classes
 */
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone, Debug)]
struct EC(u32);
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug)]
pub struct Gene(String);

// #[inline(never)]
fn update_intersection_via_retain<T:Hash+Eq>(inter:  &mut HashSet<T>, newset: &HashSet<T>){        
    // delete any elemt in shared_genes not present in current_set
    // i..e set intersection
    // we cant delete while iterating, so remember which elements to delete
    inter.retain(|item| newset.contains(item));
}

// pub fn find_consistent(records: &Vec<BusRecord>, ec2gene: &HashMap<u32, HashSet<String>>) ->HashSet<String> {
pub fn find_consistent(records: &Vec<BusRecord>, ec2gene: &Ec2GeneMapper) ->HashSet<u32> {

    /*
    set intersection in Rust is a MESS due to ownership etc!!
    Here's the strategy:
    1. copy the first set s1 entirely (the intersection will only contain elements of this set)
    2. go threought the second set s2: if an element in s2 is NOT present in s1, remove it from s1. 
    */

    // this one does shrink the set, and uses the iterator

    // get the genes per BusRecord
    // let mut setlist = records.iter().map(|r| ec2gene.get(&r.EC).unwrap());
    let mut setlist = records.iter().map(|r| ec2gene.get_genes(r.EC));

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
        if shared_genes.is_empty(){
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

    let es_mapper = Ec2GeneMapper::new(ec_dict);

    // not find_consistent doesnt check the records indeed come from the same CB/UMI

    // single read, consistent with A
    let r1 =BusRecord{CB: 0, UMI: 21, EC: 0, COUNT: 2, FLAG: 0};
    let res1: Vec<String> = find_consistent(&vec![r1], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
    assert_eq!(vec!["A"], res1);

    // single read, consistent with A and B, hence multimapped
    let r2 =BusRecord{CB: 0, UMI: 21, EC: 2, COUNT: 2, FLAG: 0};
    let res2: Vec<String> = find_consistent(&vec![r2], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
    // assert_eq!(vec!["A", "B"], res2);  WRNING the order matters here, the function might return the vec ordered differently
    assert_eq!(res2.len(), 2);

    // two reads, consistent with A
    let r3 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
    let r4 =BusRecord{CB: 3, UMI: 0, EC: 2, COUNT:  2, FLAG: 0}; 
    let res3: Vec<String> = find_consistent(&vec![r3, r4], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
    assert_eq!(vec!["A"], res3);

    // two reads, inconsistent with A, B
    let r5 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
    let r6 =BusRecord{CB: 3, UMI: 0, EC: 1, COUNT:  2, FLAG: 0}; 
    let res4: Vec<String> = find_consistent(&vec![r5, r6], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
    assert_eq!(res4.len(), 0);

    // three reads,  A, B, (A,B)
    // inconsintent at the end
    let r7 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
    let r8 =BusRecord{CB: 3, UMI: 0, EC: 1, COUNT:  2, FLAG: 0}; 
    let r9 =BusRecord{CB: 3, UMI: 0, EC: 2, COUNT:  2, FLAG: 0}; 
    let res5: Vec<String> = find_consistent(&vec![r7, r8, r9], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
    assert_eq!(res5.len(), 0);

}



#[derive(Debug)]
pub struct Ec2GeneMapper {
    /*
    A struct to help deal with the complicated relationship between 
    equivalence class (EC) and genes:
    A single EC maps to a SET of genes
    */
    ec2geneid: HashMap<EC, HashSet<u32>>,  // maps each EC to a set of geneids (integers)
    gene_to_int: HashMap<Gene, u32>,  // for a given gene (ENSBEMLE id), get its unique ID (as recorded in ec2geneid)

    // TODO this could be a simple Vec, indexing into it
    pub int_to_gene: HashMap<u32, Gene>,  // for a given gene (ENSBEMLE id), get its unique ID (as recorded in ec2geneid)

    // consistent_ec_sets: HashSet<HashSet<EC>>  // what general EC combinations are allowed
}


fn ec_combinatorics_for_gene(geneid: u32, ec2geneid: &HashMap<EC, HashSet<u32>>){
    todo!();
    // all ECs that contain the gene
    let ecs_sets: HashMap<EC, HashSet<u32>> = ec2geneid.iter()
        .filter(|(_ec, gset)| gset.contains(&geneid))
        .map(|(_ec, gset)| (*_ec, gset.clone())).collect();

    let mut current_consistent_sets :HashMap<EC, HashSet<u32>> = HashMap::new();

    // for now, every candidate is a maybe
    let mut current_maybe_sets :HashMap<EC, HashSet<u32>> = ecs_sets.clone();
    let target_set = HashSet::from_iter(vec![geneid]);

    // "recursion": just keep iterating until there's no more candidates
    while  !current_maybe_sets.is_empty() {
        
        let mut new_maybe_sets :HashMap<EC, HashSet<u32>> = HashMap::new();

        // filter out anything thats not overlaping with target at all
        for (ec, set) in current_maybe_sets{
            let overlap = set.intersection(&target_set).count();
            if overlap > 0{
                if set == target_set{ // the ones that are exacly the target set are added to the solution
                    current_consistent_sets.insert(ec, set);
                }
                else{ // it contains more elements then the target set
                    new_maybe_sets.insert(ec, set);
                }
            }
        }
        // at this point we have a bunch of candidates in new_maybe_sets
        // however, they are all too large and need to be combined with some other set
        // to possibly result in the desired intersection
        // 1. combine them with any of the current solutions: that'll be a sure 

        current_maybe_sets = new_maybe_sets;
    }

}
impl Ec2GeneMapper{

    pub fn new(ec2gene: HashMap<u32, HashSet<String>>) ->  Self{

        // get all genes and sort them
        // could be done with a flatmap?
        // let gene_list2:HashSet<Gene> =  ec2gene
        //     .values()
        //     .flat_map(|set|
        //                 set.iter().map(|s| Gene(s.to_string()))
        //             )
        //     .collect();
        
        let mut gene_list : HashSet<Gene> = HashSet::new();
        for genes in ec2gene.values(){
            for g in genes{
                gene_list.insert(Gene(g.to_string()));
            }
        }
        let mut gene_vector: Vec<Gene> = gene_list.into_iter().collect();
        gene_vector.sort();

        // we're going to encode genes by ints hence we need
        // a mapping from gene-> int
        let gene_to_int: HashMap<Gene, u32> = gene_vector
            .into_iter()
            .enumerate()
            .map(|(i, g)| (g, i as u32))
            .collect();

        // the reverse
        let int_to_gene: HashMap<u32, Gene> = gene_to_int.iter().map(|(g, i)| (*i as u32, g.clone())).collect();

        // mapping from EC to a set of gene_id
        // this is the main struct we'll be using!
        let mut ec2geneid: HashMap<EC, HashSet<u32>>= HashMap::new();
        for (ec, genes) in ec2gene{
            let geneids: HashSet<u32> = genes.iter().map(|gname| *gene_to_int.get(&Gene(gname.to_string())).unwrap()).collect();
            ec2geneid.insert(EC(ec), geneids);
        }

        // iterate through all genes, and see what ECs are consitent with it: what set of ECs (of any size)
        // could uniquely resolve into this gene
        
        Ec2GeneMapper { ec2geneid,  gene_to_int, int_to_gene}
    }

    pub fn get_genes(&self, ec: u32) -> &HashSet<u32>{
        // resolves an EC into a set of gene_ids
        return self.ec2geneid.get(&EC(ec)).unwrap();
    }

    pub fn get_genenames(&self, ec: u32) -> HashSet<String>{
        // resolves an EC into a set of gene_names
        let geneids = self.get_genes(ec);
        let genenames = geneids
            .iter()
            .map(|gid| self.resolve_gene_id(*gid))
            .collect();
        genenames
    }

    pub fn resolve_gene_id(&self, gene_id: u32) -> String{
        // turns a gene_id into a genename (ENSEMBL)
        let r = self.int_to_gene.get(&gene_id).unwrap();
        r.0.to_string()
    }

}
#[test]
fn test_ec_struct(){

    let ec0: HashSet<String> = vec!["A".to_string()].into_iter().collect();
    let ec1: HashSet<String> = vec!["B".to_string()].into_iter().collect();
    let ec2: HashSet<String> = vec!["A".to_string(), "B".to_string()].into_iter().collect();
    let ec3: HashSet<String> = vec!["C".to_string(), "D".to_string()].into_iter().collect();

    let ec_dict: HashMap<u32, HashSet<String>> = HashMap::from([
        (0, ec0.clone()), 
        (1, ec1.clone()), 
        (2, ec2.clone()), 
        (3, ec3.clone()), 
        ]);

    let es = Ec2GeneMapper::new(ec_dict);

    for (i, ec) in vec![ec0, ec1, ec2, ec3].iter().enumerate(){
        let gids = es.get_genes(i as u32);
        let gnames: HashSet<_> = gids.iter().map(|i| es.resolve_gene_id(*i)).collect();
        assert_eq!(&gnames, ec);
    }
}
