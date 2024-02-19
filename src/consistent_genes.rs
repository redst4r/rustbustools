//! This module handles the Equivalance class to gene mapping
//!
//! An EC represents a *set* of transcripts which the read is consistent with.
//! This is then mapped to genes:
//! * the EC resolves to a single gene, clear case for a count
//! * the EC resolves to multipel genes (the read aligned to a part of the transcriptome which is ambigous)
//!
//! It gets more complicated if we have multiple ECs for an mRNA,
//! i.e. multiple busrecords with same CB/UMI, but different EC
//! This happens since mRNAs get fragmented after amplification, and different fragments can map to different parts \
//! of the transcriptome, hence yielding different ECs
//!
//! See [Ec2GeneMapper] and [MappingResult]
//!
use crate::io::BusFolder;
use crate::{disjoint::Intersector, io::BusRecord};
use itertools::izip;
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::hash::Hash;
use std::io::{BufReader, BufRead};

/*
just some wrappers for Strings and ints that we use for genes and equivalence classes
 */
/// Thin wrapper aroud u32, representing an Equivalence Class
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone, Debug)]
pub struct EC(pub u32);

/// Gene identifier
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone, Debug)]
pub struct GeneId(pub u32);

/// Thin wrapper around u64, a cell-barcode
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone, Debug)]
pub struct CB(pub u64);

/// Name of a particular gene
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug)]
pub struct Genename(pub String);

/// intersecting both sets, modifying the first one (it'll stay the same or loose elements)
fn update_intersection_via_retain<T: Hash + Eq>(inter: &mut HashSet<T>, newset: &HashSet<T>) {
    // delete any elemt in shared_genes not present in current_set
    // i..e set intersection
    // we cant delete while iterating, so remember which elements to delete
    inter.retain(|item| newset.contains(item));
}

/// MappingResult represents an attempt to unify several BusRecords with matching CB/UMI
/// into an mRNA expressed from a single gene.
///
/// This is not always possible, as the records will have different EC
/// Sometimes the ECs are inconsistent (no overlap in gene space)
/// Sometimes the ECS are ambigous and can resolve to multipel genes
/// TODO: MappingResult isnt the best name (mapping sounds like alignment!)
///
#[derive(Debug, PartialEq)]
pub enum MappingResult {
    /// records mapped to a single Gene
    SingleGene(GeneId),
    /// record multimapped, i.e. records are consistent with multiple genes
    Multimapped(HashSet<GeneId>),
    /// records are inconsistent, pointing to different genes
    Inconsistent,
}

pub enum MappingMode {
    EC(InconsistentResolution), // just the EC as the molecular identiy of a read
    // IgnoreMultipleCbUmi,  // if we come arcoss a CB/UMI with more than one busrecord (,i.e. multiple ECs), just skip that cb/umi
    Gene(Ec2GeneMapper, InconsistentResolution),  // use the gene 
}

/// if we come across a CB/UMI the has inconsistent mapping
/// i.e. mapping to 2 different genes, how to handle
#[derive(Debug)]
pub enum InconsistentResolution {
    IgnoreInconsistent, // simply ignore the entire CB/UMI (as if we'd never seen it)
    AsDistinct, // treat as distinct entities, i.e. a chance collision of two different molecules (mRNAs)
    AsSingle,  // just ignore the fact that the moelcule might be different, aggregate all counts
}


/// For a set of busrecords (coming from the same molecule), this function tries
/// to map those records to genes consistently.
///
///  As all records come from the same mRNA,ideally we get a consistent
/// mapping (e.g. all records map to the same gene)
/// See [MappingResult] for mor explanation
pub fn find_consistent(records: &[BusRecord], ec2gene: &Ec2GeneMapper) -> MappingResult {
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

    let mut shared_genes = s1.clone();

    // save some time if its only one record
    // actually not sure if this saves anything any more (after using MappingResult)
    // if we have one record, shared_Genes will be len==1, we pull out s1 (shared_Genes is empty now)
    // -> check if we can get rid of this block!
    if records.len() == 1 {
        if shared_genes.len() == 1 {
            let elem = *shared_genes.iter().next().unwrap();
            return MappingResult::SingleGene(elem);
        } else {
            return MappingResult::Multimapped(shared_genes);
        }
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
    match shared_genes.len() {
        0 => MappingResult::Inconsistent,
        1 => {
            let elem = *shared_genes.iter().next().unwrap();
            MappingResult::SingleGene(elem)
        }
        _ => MappingResult::Multimapped(shared_genes),
    }
}

/// Deals with the EC (equivalence class) to gene mapping
/// Resolve a given EC into a set of genes consistent with that EC
/// # Example
/// ```rust, no_run
/// # use bustools::io::BusFolder;
/// # use std::collections::HashSet;
/// # use bustools::consistent_genes::{Genename, EC};
/// let bfolder = BusFolder::new("/path/to/busfolder");
/// let ec2g = bfolder.make_mapper("/path/to/transcripts_to_genes.txt");
/// let ec = EC(1234);
/// let genenames: HashSet<Genename> = ec2g.get_genenames(ec);
/// ```
#[derive(Debug, Clone)]
pub struct Ec2GeneMapper {
    /*
    A struct to help deal with the complicated relationship between
    equivalence class (EC) and genes:
    A single EC maps to a SET of genes
    */
    ec2geneid: HashMap<EC, HashSet<GeneId>>, // maps each EC to a set of geneids (integers)

    // ec2geneid_array: Vec<HashSet<GeneId>>, // maps each EC to a set of geneids (integers)
    //
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
    /// make an Ec2GeneMapper from a HashMap of EC-> set(Genenames)
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
        for (ec, genes) in ec2gene.iter() {
            if genes.len() == 1 {
                let gname = genes.iter().next().unwrap().clone();
                let gid = *gene_to_int.get(&gname).unwrap();
                geneid2ec.insert(gid, *ec);
            }
        }

        // //faster
        // // need to mark ec2gene ordered first
        // let mut ordered_ec2gene: BTreeMap<EC, HashSet<Genename>> = BTreeMap::new();
        // for (k,v) in ec2gene {
        //     ordered_ec2gene.insert(k ,v);
        // }
        // let mut _vec = Vec::with_capacity(ordered_ec2gene.len());
        // for (ec, genes) in ordered_ec2gene {
        //     let geneids: HashSet<GeneId> = genes
        //         .iter()
        //         .map(|gname| *gene_to_int.get(gname).unwrap())
        //         .collect();
        //     _vec.push(geneids);
        // }

        Ec2GeneMapper { ec2geneid, int_to_gene, geneid2ec }
    }

    /// resolves an EC into a set of gene_ids
    pub fn get_genes(&self, ec: EC) -> &HashSet<GeneId> {
        self.ec2geneid.get(&ec).unwrap()
        // self.ec2geneid_array.get(ec.0 as usize).unwrap()
    }

    /// returns a set of genenames (usually ENSG) consistent with EC
    pub fn get_genenames(&self, ec: EC) -> HashSet<Genename> {
        // resolves an EC into a set of gene_names
        let geneids = self.get_genes(ec);
        let genenames = geneids
            .iter()
            .map(|gid| self.resolve_gene_id(*gid))
            .collect();
        genenames
    }

    /// resolves a gene_id into its genename (ENSEMBL)
    pub fn resolve_gene_id(&self, gene_id: GeneId) -> Genename {
        let r = self.int_to_gene.get(&gene_id).unwrap();
        r.clone()
    }

    /// returns a list of all genes in the mapping, sorted by int-id
    pub fn get_gene_list(&self) -> Vec<Genename> {
        let ngenes = self.int_to_gene.len();
        let genelist_vector: Vec<Genename> = (0..ngenes)
            .map(|k| self.resolve_gene_id(GeneId(k as u32)))
            .collect();
        genelist_vector
    }

    fn resolve_geneid_to_ec_uniquely(&self, geneid: u32) -> Option<EC> {
        // just dereferencing whats inside the Option returned by the dict
        self.geneid2ec.get(&GeneId(geneid)).copied()
    }
}

pub (crate) fn make_mapper(busfolder: &BusFolder, t2g_file: &str) -> Ec2GeneMapper{
    let t2g_dict = parse_t2g(t2g_file);
    let e2g = build_ec2gene(
        &busfolder.parse_ecmatrix(), 
        &busfolder.parse_transcript(),
        &t2g_dict);
    Ec2GeneMapper::new(e2g)
}


fn build_ec2gene(
    ec_dict: &HashMap<EC, Vec<u32>>,
    transcript_dict: &HashMap<u32, String>,
    t2g_dict: &HashMap<String, Genename>,
) -> HashMap<EC, HashSet<Genename>> {
    let mut ec2gene: HashMap<EC, HashSet<Genename>> = HashMap::new();

    for (ec, transcript_ints) in ec_dict.iter() {
        let mut genes: HashSet<Genename> = HashSet::new();

        for t_int in transcript_ints {
            let t_name = transcript_dict.get(t_int).unwrap();

            // if we can resolve, put the genename, otherwsise use the transcript name instead
            //
            // actually, turns out that kallisto/bustools treats it differently:
            // if the transcript doenst resolve, just drop tha trasncrip from the EC set
            // TODO what happens if non of the EC transcripts resolve
            if let Some(genename) = t2g_dict.get(t_name) {
                genes.insert(genename.clone());
            }
            // else { genes.insert(Genename(t_name.clone())); }
        }
        ec2gene.insert(*ec, genes);
    }
    // // sanity check, make sure no Ec set is empty
    // for (ec, gset) in ec2gene.iter(){
    //     if gset.is_empty(){
    //         println!("{ec:?}'s geneset is empty");
    //     }
    // }

    ec2gene
}

fn parse_t2g(t2g_file: &str) -> HashMap<String, Genename> {
    let mut t2g_dict: HashMap<String, Genename> = HashMap::new();
    let file = File::open(t2g_file).unwrap_or_else(|_| panic!("{} not found", t2g_file));
    let reader = BufReader::new(file);
    for line in reader.lines() {
        if let Ok(l) = line {
            let mut s = l.split_whitespace();
            let transcript_id = s.next().unwrap();
            let ensemble_id = s.next().unwrap();
            let _symbol = s.next().unwrap();

            assert!(!t2g_dict.contains_key(&transcript_id.to_string())); //make sure transcripts dont map to multiple genes
            t2g_dict.insert(transcript_id.to_string(), Genename(ensemble_id.to_string()));
        }
    }
    t2g_dict
}



/// Represents a busrecord (actually an observed CB/UIM combination)
/// with a set of consistent genes (genes that this CB/UMI could map to)
/// and the number of times this combination was seen in the busfile
///
#[derive(Debug)]
#[allow(non_snake_case)]
pub struct CUGset {
    pub CB: u64,                    //8byte
    pub UMI: u64,                   // 8byte
    pub GENESET: HashSet<Genename>, // 4v byte
    pub COUNT: u32,                 // 4v byte
}

/// Group a set of busrecords (same CB/UMI) by genes they are consistent with
pub fn groubygene(records: Vec<BusRecord>, ec2gene: &Ec2GeneMapper) -> Vec<CUGset> {
    let mut emissions: Vec<CUGset> = Vec::with_capacity(records.len()); //with capacity worst case scenario

    // aggregate by overlpping genes
    let mut inter: Intersector<Genename, BusRecord> = Intersector::new();
    for r in records {
        let gset = ec2gene.get_genenames(EC(r.EC));
        inter.add(gset, r);
    }

    // for (gset, grouped_records) in inter.iterate_items(){
    for (gset, grouped_records) in izip!(inter.keys, inter.items) {
        let counts: u32 = grouped_records.iter().map(|x| x.COUNT).sum();
        let r1 = grouped_records.first().unwrap();

        let new_record = CUGset { CB: r1.CB, UMI: r1.UMI, GENESET: gset, COUNT: counts };
        emissions.push(new_record);
    }
    emissions
}

#[cfg(test)]
mod testing {
    use crate::{
        consistent_genes::{find_consistent, groubygene, Genename, MappingResult},
        io::{BusFolder, BusRecord},
        utils::vec2set,
    };
    use std::collections::{HashMap, HashSet};

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
        let r1 = BusRecord { CB: 0, UMI: 21, EC: 0, COUNT: 2, FLAG: 0 };
        let res1: MappingResult = find_consistent(&vec![r1], &es_mapper);
        match res1 {
            MappingResult::SingleGene(g) => {
                assert_eq!(Genename("A".to_string()), es_mapper.resolve_gene_id(g))
            }
            MappingResult::Multimapped(_) | MappingResult::Inconsistent => panic!(),
        }

        // single read, consistent with A and B, hence multimapped
        let r2 = BusRecord { CB: 0, UMI: 21, EC: 2, COUNT: 2, FLAG: 0 };
        let res2: MappingResult = find_consistent(&vec![r2], &es_mapper);
        println!("{:?}", res2);
        match res2 {
            MappingResult::SingleGene(_) | MappingResult::Inconsistent => panic!(),
            MappingResult::Multimapped(gs) => {
                let obs = gs
                    .iter()
                    .map(|g| es_mapper.resolve_gene_id(*g))
                    .collect::<HashSet<_>>();
                assert_eq!(
                    obs,
                    vec2set(vec![Genename("A".to_string()), Genename("B".to_string())])
                )
            }
        }

        // two reads, consistent with A
        let r3 = BusRecord { CB: 1, UMI: 3, EC: 0, COUNT: 2, FLAG: 0 };
        let r4 = BusRecord { CB: 3, UMI: 0, EC: 2, COUNT: 2, FLAG: 0 };
        let res3: MappingResult = find_consistent(&vec![r3, r4], &es_mapper);
        match res3 {
            MappingResult::SingleGene(g) => {
                assert_eq!(Genename("A".to_string()), es_mapper.resolve_gene_id(g))
            }
            MappingResult::Multimapped(_) | MappingResult::Inconsistent => panic!(),
        }

        // // two reads, inconsistent with A, B
        let r5 = BusRecord { CB: 1, UMI: 3, EC: 0, COUNT: 2, FLAG: 0 };
        let r6 = BusRecord { CB: 3, UMI: 0, EC: 1, COUNT: 2, FLAG: 0 };
        let res4 = find_consistent(&vec![r5, r6], &es_mapper);
        assert_eq!(res4, MappingResult::Inconsistent);

        // // three reads,  A, B, (A,B)
        // // inconsintent at the end
        let r7 = BusRecord { CB: 1, UMI: 3, EC: 0, COUNT: 2, FLAG: 0 };
        let r8 = BusRecord { CB: 3, UMI: 0, EC: 1, COUNT: 2, FLAG: 0 };
        let r9 = BusRecord { CB: 3, UMI: 0, EC: 2, COUNT: 2, FLAG: 0 };
        let res5 = find_consistent(&vec![r7, r8, r9], &es_mapper);
        assert_eq!(res5, MappingResult::Inconsistent)
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
    fn test_get_gene_list() {
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
        assert_eq!(
            es.get_gene_list(),
            vec![
                Genename("A".to_string()),
                Genename("B".to_string()),
                Genename("C".to_string()),
                Genename("D".to_string())
            ]
        )
    }

    #[test]
    fn test_resolve_geneid_2_ec() {
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

        // Gene A resolves to EC0
        assert_eq!(es.resolve_geneid_to_ec_uniquely(0), Some(EC(0)));
        // Gene Bresolves to EC1
        assert_eq!(es.resolve_geneid_to_ec_uniquely(1), Some(EC(1)));

        // Gene C cant be resolved uniquely (EC3 mapps to both C,D)
        assert_eq!(es.resolve_geneid_to_ec_uniquely(2), None);
        // Gene D cant be resolved uniquely (EC3 mapps to both C,D)
        assert_eq!(es.resolve_geneid_to_ec_uniquely(3), None);
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
        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 1, EC: 2, COUNT: 2, FLAG: 0 };

        let r = groubygene(vec![r1, r2], &es);
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
        let s1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0 }; // A
        let s2 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 4, FLAG: 0 }; //B
        let r = groubygene(vec![s1, s2], &es);
        assert_eq!(r.len(), 2);
        let r1 = &r[0];
        let r2 = &r[1];
        assert_eq!(r1.COUNT, 3);
        assert_eq!(r2.COUNT, 4);
    }

    // #[test]
    #[allow(dead_code)]
    fn test_ec() {
        let folder = "/home/michi/mounts/TB4drive/ISB_data/201015_NS500720_0063_AHV53GBGXG/kallisto_quant/01_Day2/kallisto/sort_bus/bus_output/";
        let t2g_file = "/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt";
        let b = BusFolder::new(folder);
        let ecmapper = b.make_mapper(t2g_file);

        for g in 0..ecmapper.int_to_gene.len() {
            ecmapper.resolve_geneid_to_ec_uniquely(g as u32).unwrap();
            // let ec = ecmapper.geneid2ec.get(&(g as u32)).unwrap().0;
        }
    }
}
