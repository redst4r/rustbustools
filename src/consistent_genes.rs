use itertools::izip;
use std::collections::{HashMap, HashSet};

use crate::{disjoint::Intersector, io::BusRecord};
// use crate::disjoint::DisjointSubsets;
use std::hash::Hash;

/*
just some wrappers for Strings and ints that we use for genes and equivalence classes
 */
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone, Debug)]
pub struct EC(pub u32);

#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone, Debug)]
pub struct GeneId(pub u32);

#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone, Debug)]
pub struct CB(pub u64);

#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug)]
pub struct Genename(pub String);

// #[inline(never)]
pub fn update_intersection_via_retain<T: Hash + Eq>(inter: &mut HashSet<T>, newset: &HashSet<T>) {
    // delete any elemt in shared_genes not present in current_set
    // i..e set intersection
    // we cant delete while iterating, so remember which elements to delete
    inter.retain(|item| newset.contains(item));
}

// pub fn find_consistent(records: &Vec<BusRecord>, ec2gene: &HashMap<u32, HashSet<String>>) ->HashSet<String> {
pub fn find_consistent(records: &[BusRecord], ec2gene: &Ec2GeneMapper) -> HashSet<GeneId> {
    /*
    set intersection in Rust is a MESS due to ownership etc!!
    Here's the strategy:
    1. copy the first set s1 entirely (the intersection will only contain elements of this set)
    2. go threought the second set s2: if an element in s2 is NOT present in s1, remove it from s1.
    */

    // this one does shrink the set, and uses the iterator

    // get the genes per BusRecord
    // let mut setlist = records.iter().map(|r| ec2gene.get(&r.EC).unwrap());
    let mut setlist = records.iter().map(|r| ec2gene.get_genes(EC(r.EC)));

    let s1 = setlist.next().unwrap();

    let mut shared_genes;
    if false {
        shared_genes = HashSet::new(); // could spec capacity

        // initial element, make a copy of that
        for el in s1 {
            shared_genes.insert(*el);
        }
    }
    // pretty much just a clone of s1
    else {
        shared_genes = s1.clone();
    }

    // save some time if its only one record
    if records.len() == 1 {
        return shared_genes;
    }

    for current_set in setlist {
        // delete any elemt in shared_genes not present in current_set
        // i..e set intersection
        update_intersection_via_retain(&mut shared_genes, current_set);

        // to save some time: if the current intersection is already empty, it'll stay empty
        if shared_genes.is_empty() {
            break;
        }
    }
    shared_genes
}

#[derive(Debug)]
pub struct Ec2GeneMapper {
    /*
    A struct to help deal with the complicated relationship between
    equivalence class (EC) and genes:
    A single EC maps to a SET of genes
    */
    ec2geneid: HashMap<EC, HashSet<GeneId>>, // maps each EC to a set of geneids (integers)

    // TODO this could be a simple Vec, indexing into it
    int_to_gene: HashMap<GeneId, Genename>, // for a given gene (ENSBEMLE id), get its unique ID (as recorded in ec2geneid)

    // for each gene get a EC that uniquely maps to that gene
    geneid2ec: HashMap<GeneId, EC>, // consistent_ec_sets: HashSet<HashSet<EC>>  // what general EC combinations are allowed
}

// fn ec_combinatorics_for_gene(geneid: u32, ec2geneid: &HashMap<EC, HashSet<u32>>){
//     todo!();
//     // all ECs that contain the gene
//     let ecs_sets: HashMap<EC, HashSet<u32>> = ec2geneid.iter()
//         .filter(|(_ec, gset)| gset.contains(&geneid))
//         .map(|(_ec, gset)| (*_ec, gset.clone())).collect();

//     let mut current_consistent_sets :HashMap<EC, HashSet<u32>> = HashMap::new();

//     // for now, every candidate is a maybe
//     let mut current_maybe_sets :HashMap<EC, HashSet<u32>> = ecs_sets.clone();
//     let target_set = HashSet::from_iter(vec![geneid]);

//     // "recursion": just keep iterating until there's no more candidates
//     while  !current_maybe_sets.is_empty() {

//         let mut new_maybe_sets :HashMap<EC, HashSet<u32>> = HashMap::new();

//         // filter out anything thats not overlaping with target at all
//         for (ec, set) in current_maybe_sets{
//             let overlap = set.intersection(&target_set).count();
//             if overlap > 0{
//                 if set == target_set{ // the ones that are exacly the target set are added to the solution
//                     current_consistent_sets.insert(ec, set);
//                 }
//                 else{ // it contains more elements then the target set
//                     new_maybe_sets.insert(ec, set);
//                 }
//             }
//         }
//         // at this point we have a bunch of candidates in new_maybe_sets
//         // however, they are all too large and need to be combined with some other set
//         // to possibly result in the desired intersection
//         // 1. combine them with any of the current solutions: that'll be a sure

//         current_maybe_sets = new_maybe_sets;
//     }

// }
impl Ec2GeneMapper {
    pub fn new(ec2gene: HashMap<EC, HashSet<Genename>>) -> Self {
        // get all genes and sort them
        // could be done with a flatmap?
        // let gene_list2:HashSet<Gene> =  ec2gene
        //     .values()
        //     .flat_map(|set|
        //                 set.iter().map(|s| Gene(s.to_string()))
        //             )
        //     .collect();

        let mut gene_list: HashSet<Genename> = HashSet::new();
        for genes in ec2gene.values() {
            for g in genes {
                gene_list.insert(g.clone());
            }
        }
        let mut gene_vector: Vec<Genename> = gene_list.into_iter().collect();
        gene_vector.sort();

        // we're going to encode genes by ints hence we need
        // a mapping from gene-> int
        let gene_to_int: HashMap<Genename, GeneId> = gene_vector
            .into_iter()
            .enumerate()
            .map(|(i, g)| (g, GeneId(i as u32)))
            .collect();

        // the reverse
        let int_to_gene: HashMap<GeneId, Genename> =
            gene_to_int.iter().map(|(g, i)| (*i, g.clone())).collect();

        // mapping from EC to a set of gene_id
        // this is the main struct we'll be using!
        let mut ec2geneid: HashMap<EC, HashSet<GeneId>> = HashMap::new();
        for (ec, genes) in ec2gene.iter() {
            let geneids: HashSet<GeneId> = genes
                .iter()
                .map(|gname| *gene_to_int.get(gname).unwrap())
                .collect();
            ec2geneid.insert(*ec, geneids);
        }

        let mut geneid2ec: HashMap<GeneId, EC> = HashMap::new();
        for (ec, genes) in ec2gene {
            if genes.len() == 1 {
                let gname = genes.iter().next().unwrap().clone();
                let gid = *gene_to_int.get(&gname).unwrap();
                geneid2ec.insert(gid, ec);
            }
        }

        Ec2GeneMapper {
            ec2geneid,
            int_to_gene,
            geneid2ec,
        }
    }

    pub fn get_genes(&self, ec: EC) -> &HashSet<GeneId> {
        // resolves an EC into a set of gene_ids
        self.ec2geneid.get(&ec).unwrap()
    }

    pub fn get_genenames(&self, ec: EC) -> HashSet<Genename> {
        // resolves an EC into a set of gene_names
        let geneids = self.get_genes(ec);
        let genenames = geneids
            .iter()
            .map(|gid| self.resolve_gene_id(*gid))
            .collect();
        genenames
    }

    pub fn resolve_gene_id(&self, gene_id: GeneId) -> Genename {
        // turns a gene_id into a genename (ENSEMBL)
        let r = self.int_to_gene.get(&gene_id).unwrap();
        r.clone()
    }

    pub fn get_gene_list(&self) -> Vec<Genename> {
        // returns a list of all genes in the mapping, sorted by int-id
        let ngenes = self.int_to_gene.len();
        let genelist_vector: Vec<Genename> = (0..ngenes)
            .map(|k| self.resolve_gene_id(GeneId(k as u32)))
            .collect();
        genelist_vector
    }

    pub fn resolve_geneid_to_ec_uniquely(&self, geneid: u32) -> Option<EC> {
        // just dereferencing whats inside the Option returned by the dict
        // self.geneid2ec.get(&geneid).map(|ec| *ec)
        self.geneid2ec.get(&GeneId(geneid)).copied()
    }
}

#[derive(Debug)]
#[allow(non_snake_case)]
pub struct CUGset {
    pub CB: u64,                    //8byte
    pub UMI: u64,                   // 8byte
    pub GENESET: HashSet<Genename>, // 4v byte
    pub COUNT: u32,                 // 4v byte
}

pub fn groubygene(records: Vec<BusRecord>, ec2gene: &Ec2GeneMapper) -> Vec<CUGset> {
    // build disjoint set based on samplename and genes
    // let mut disjoint_set = DisjointSubsets::new();
    // for (id, r) in records.iter().enumerate(){
    //     let gset = ec2gene.get_genenames(r.EC);
    //     disjoint_set.add(id.to_string(), gset);
    // }
    let mut emissions: Vec<CUGset> = Vec::with_capacity(records.len()); //with capacity worst case scenario

    // aggregate by overlpping genes
    let mut inter: Intersector<Genename, BusRecord> = Intersector::new();
    for r in records {
        let gset = ec2gene.get_genenames(EC(r.EC));
        inter.add(gset, r);
    }

    for (gset, grouped_records) in izip!(inter.keys, inter.items) {
        let counts: u32 = grouped_records.iter().map(|x| x.COUNT).sum();
        let r1 = grouped_records.get(0).unwrap();

        let new_record = CUGset {
            CB: r1.CB,
            UMI: r1.UMI,
            GENESET: gset,
            COUNT: counts,
        };
        emissions.push(new_record);
    }
    emissions
}

#[cfg(test)]
mod testing {
    use std::collections::{HashMap, HashSet};

    use crate::{
        consistent_genes::{find_consistent, groubygene, Genename},
        io::{BusFolder, BusRecord},
        utils::vec2set,
    };

    use super::{Ec2GeneMapper, EC};

    #[test]
    fn test_consistent() {
        let ec0 = vec2set(vec![Genename("A".to_string())]);
        let ec1 = vec2set(vec![Genename("B".to_string())]);
        let ec2 = vec2set(vec![Genename("A".to_string()), Genename("B".to_string())]);
        let ec3 = vec2set(vec![Genename("C".to_string()), Genename("D".to_string())]);

        let ec_dict: HashMap<EC, HashSet<Genename>> =
            HashMap::from([(EC(0), ec0), (EC(1), ec1), (EC(2), ec2), (EC(3), ec3)]);

        let es_mapper = Ec2GeneMapper::new(ec_dict);

        // not find_consistent doesnt check the records indeed come from the same CB/UMI

        // single read, consistent with A
        let r1 =BusRecord{CB: 0, UMI: 21, EC: 0, COUNT: 2, FLAG: 0};
        let res1: Vec<Genename> = find_consistent(&vec![r1], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
        assert_eq!(vec![Genename("A".to_string())], res1);

        // single read, consistent with A and B, hence multimapped
        let r2 =BusRecord{CB: 0, UMI: 21, EC: 2, COUNT: 2, FLAG: 0};
        let res2: Vec<Genename> = find_consistent(&vec![r2], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
        // assert_eq!(vec!["A", "B"], res2);  WRNING the order matters here, the function might return the vec ordered differently
        assert_eq!(res2.len(), 2);

        // two reads, consistent with A
        let r3 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
        let r4 =BusRecord{CB: 3, UMI: 0, EC: 2, COUNT:  2, FLAG: 0}; 
        let res3: Vec<Genename> = find_consistent(&vec![r3, r4], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
        assert_eq!(vec![Genename("A".to_string())], res3);

        // two reads, inconsistent with A, B
        let r5 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
        let r6 =BusRecord{CB: 3, UMI: 0, EC: 1, COUNT:  2, FLAG: 0}; 
        let res4: Vec<Genename> = find_consistent(&vec![r5, r6], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
        assert_eq!(res4.len(), 0);

        // three reads,  A, B, (A,B)
        // inconsintent at the end
        let r7 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
        let r8 =BusRecord{CB: 3, UMI: 0, EC: 1, COUNT:  2, FLAG: 0}; 
        let r9 =BusRecord{CB: 3, UMI: 0, EC: 2, COUNT:  2, FLAG: 0}; 
        let res5: Vec<Genename> = find_consistent(&vec![r7, r8, r9], &es_mapper).into_iter().map(|gid|es_mapper.resolve_gene_id(gid)).collect();
        assert_eq!(res5.len(), 0);
    }

    #[test]
    fn test_ec_struct() {
        let ec0 = vec2set(vec![Genename("A".to_string())]);
        let ec1 = vec2set(vec![Genename("B".to_string())]);
        let ec2 = vec2set(vec![Genename("A".to_string()), Genename("B".to_string())]);
        let ec3 = vec2set(vec![Genename("C".to_string()), Genename("D".to_string())]);

        let ec_dict: HashMap<EC, HashSet<Genename>> = HashMap::from([
            (EC(0), ec0.clone()),
            (EC(1), ec1.clone()),
            (EC(2), ec2.clone()),
            (EC(3), ec3.clone()),
        ]);

        let es = Ec2GeneMapper::new(ec_dict);

        for (i, ec_set) in vec![ec0, ec1, ec2, ec3].iter().enumerate() {
            let gids = es.get_genes(EC(i as u32));
            let gnames: HashSet<_> = gids.iter().map(|i| es.resolve_gene_id(*i)).collect();
            assert_eq!(&gnames, ec_set);
        }
    }

    #[test]
    fn test_groubygene() {
        let ec0: HashSet<Genename> = vec2set(vec![Genename("A".to_string())]);
        let ec1: HashSet<Genename> = vec2set(vec![Genename("B".to_string())]);
        let ec2: HashSet<Genename> =
            vec2set(vec![Genename("A".to_string()), Genename("B".to_string())]);
        let ec3: HashSet<Genename> =
            vec2set(vec![Genename("C".to_string()), Genename("D".to_string())]);

        let ec_dict: HashMap<EC, HashSet<Genename>> = HashMap::from([
            (EC(0), ec0.clone()),
            (EC(1), ec1.clone()),
            (EC(2), ec2.clone()),
            (EC(3), ec3.clone()),
        ]);

        let es = Ec2GeneMapper::new(ec_dict);

        // first sample: two records, with consistent gene A
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let r2 =BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 2, FLAG: 0};
        
        let r = groubygene(vec![r1,r2], &es);
        assert_eq!(r.len(), 1);
        let r1 = r.first().unwrap();
        assert_eq!(r1.COUNT, 4);
        assert_eq!(
            r1.GENESET,
            vec![Genename("A".to_string())]
                .into_iter()
                .collect::<HashSet<Genename>>()
        );

        // second sample: two records, with consistent gene A the other consistent with gene B
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0}; // A
        let s2 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 4, FLAG: 0}; //B
        let r = groubygene(vec![s1,s2], &es);
        assert_eq!(r.len(), 2);
        let r1 = &r[0];
        let r2 = &r[1];
        assert_eq!(r1.COUNT, 3);
        assert_eq!(r2.COUNT, 4);
    }

    // #[test]
    fn test_ec() {
        let folder = "/home/michi/mounts/TB4drive/ISB_data/201015_NS500720_0063_AHV53GBGXG/kallisto_quant/01_Day2/kallisto/sort_bus/bus_output/";
        let t2g_file = "/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt";
        let b = BusFolder::new(folder, t2g_file);
        let ecmapper = b.ec2gene;

        for g in 0..ecmapper.int_to_gene.len() {
            ecmapper.resolve_geneid_to_ec_uniquely(g as u32).unwrap();
            // let ec = ecmapper.geneid2ec.get(&(g as u32)).unwrap().0;
        }
    }
}
