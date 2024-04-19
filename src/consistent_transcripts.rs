//! 
//! Mapping ECs to transcripts
//! 
use std::collections::{HashMap, HashSet};
use crate::{consistent_genes::{update_intersection_via_retain, EC}, io::{BusFolder, BusRecord}};

#[derive(Debug, PartialEq)]
pub enum MappingResultTranscript {
    /// records mapped to a single Gene
    SingleTranscript(TranscriptId),
    /// record multimapped, i.e. records are consistent with multiple genes
    Multimapped(HashSet<TranscriptId>),
    /// records are inconsistent, pointing to different genes
    Inconsistent,
}

/// transcript identifier
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Copy, Clone, Debug)]
pub struct TranscriptId(pub u32);
/// Name of a particular gene
#[derive(Eq, PartialEq, Hash, Ord, PartialOrd, Clone, Debug)]
pub struct Transcriptname(pub String);


/// Deals with the EC (equivalence class) to transcript mapping
/// Resolve a given EC into a set of transcripts consistent with that EC
#[derive(Debug, Clone)]
pub struct Ec2TranscriptMapper {
    /*
    A struct to help deal with the complicated relationship between
    equivalence class (EC) and genes:
    A single EC maps to a SET of genes
    */
    ec2tid: HashMap<EC, HashSet<TranscriptId>>, // maps each EC to a set of geneids (integers)

    // ec2geneid_array: Vec<HashSet<GeneId>>, // maps each EC to a set of geneids (integers)
    //
    // TODO this could be a simple Vec, indexing into it
    int_to_transcript: HashMap<TranscriptId, Transcriptname>, // for a given gene (ENSBEMLE id), get its unique ID (as recorded in ec2geneid)
}

impl Ec2TranscriptMapper {
    /// make an Ec2GeneMapper from a HashMap of EC-> set(Genenames)
    pub fn new(ec2transcript: HashMap<EC, HashSet<Transcriptname>>) -> Self {

        let mut transcript_list: HashSet<Transcriptname> = HashSet::new();
        for transcripts in ec2transcript.values() {
            for t in transcripts {
                transcript_list.insert(t.clone());
            }
        }
        let mut transcript_vector: Vec<Transcriptname> = transcript_list.into_iter().collect();
        transcript_vector.sort();

        // we're going to encode genes by ints hence we need
        // a mapping from gene-> int
        let transcript_to_int: HashMap<Transcriptname, TranscriptId> = transcript_vector
            .into_iter()
            .enumerate()
            .map(|(i, g)| (g, TranscriptId(i as u32)))
            .collect();

        // the reverse
        let int_to_transcript: HashMap<TranscriptId, Transcriptname> =
            transcript_to_int.iter().map(|(g, i)| (*i, g.clone())).collect();

        // mapping from EC to a set of transcript ids
        // this is the main struct we'll be using!
        let mut ec2tid: HashMap<EC, HashSet<TranscriptId>> = HashMap::new();
        for (ec, genes) in ec2transcript.iter() {
            let geneids: HashSet<TranscriptId> = genes
                .iter()
                .map(|gname| *transcript_to_int.get(gname).unwrap())
                .collect();
            ec2tid.insert(*ec, geneids);
        }
        Ec2TranscriptMapper { ec2tid, int_to_transcript }
    }

    /// resolves an EC into a set of gene_ids
    pub fn get_transcripts(&self, ec: EC) -> &HashSet<TranscriptId> {
        self.ec2tid.get(&ec).unwrap()
        // self.ec2geneid_array.get(ec.0 as usize).unwrap()
    }

    /// returns a set of genenames (usually ENSG) consistent with EC
    pub fn get_genenames(&self, ec: EC) -> HashSet<Transcriptname> {
        // resolves an EC into a set of gene_names
        let tids = self.get_transcripts(ec);
        let genenames = tids
            .iter()
            .map(|gid| self.resolve_tid(*gid))
            .collect();
        genenames
    }

    /// resolves a gene_id into its genename (ENSEMBL)
    pub fn resolve_tid(&self, tid: TranscriptId) -> Transcriptname {
        let r = self.int_to_transcript.get(&tid).unwrap();
        r.clone()
    }

    /// returns a list of all genes in the mapping, sorted by int-id
    pub fn get_transcript_list(&self) -> Vec<Transcriptname> {
        let ntrans = self.int_to_transcript.len();
        let transcriptlist_vector: Vec<Transcriptname> = (0..ntrans)
            .map(|k| self.resolve_tid(TranscriptId(k as u32)))
            .collect();
        transcriptlist_vector
    }

    // fn resolve_tid_to_ec_uniquely(&self, geneid: u32) -> Option<EC> {
    //     // just dereferencing whats inside the Option returned by the dict
    //     self.geneid2ec.get(&GeneId(geneid)).copied()
    // }
}

// For a set of busrecords (coming from the same molecule), this function tries
/// to map those records to transcripts consistently.
///
///  As all records come from the same mRNA,ideally we get a consistent
/// mapping (e.g. all records map to the same transcripts)
/// See [MappingResultTranscript] for mor explanation
pub fn find_consistent_transcripts(records: &[BusRecord], ec2gene: &Ec2TranscriptMapper) -> MappingResultTranscript {
    /*
    set intersection in Rust is a MESS due to ownership etc!!
    Here's the strategy:
    1. copy the first set s1 entirely (the intersection will only contain elements of this set)
    2. go threought the second set s2: if an element in s2 is NOT present in s1, remove it from s1.
    */

    // this one does shrink the set, and uses the iterator

    // get the genes per BusRecord
    // let mut setlist = records.iter().map(|r| ec2gene.get(&r.EC).unwrap());
    let mut setlist = records.iter().map(|r| ec2gene.get_transcripts(EC(r.EC)));

    let s1 = setlist.next().unwrap();

    let mut shared_transcripts = s1.clone();

    // save some time if its only one record
    // actually not sure if this saves anything any more (after using MappingResult)
    // if we have one record, shared_transcripts will be len==1, we pull out s1 (shared_transcripts is empty now)
    // -> check if we can get rid of this block!
    if records.len() == 1 {
        if shared_transcripts.len() == 1 {
            let elem = *shared_transcripts.iter().next().unwrap();
            return MappingResultTranscript::SingleTranscript(elem);
        } else {
            return MappingResultTranscript::Multimapped(shared_transcripts);
        }
    }

    for current_set in setlist {
        // delete any elemt in shared_transcripts not present in current_set
        // i..e set intersection
        update_intersection_via_retain(&mut shared_transcripts, current_set);

        // to save some time: if the current intersection is already empty, it'll stay empty
        if shared_transcripts.is_empty() {
            break;
        }
    }
    match shared_transcripts.len() {
        0 => MappingResultTranscript::Inconsistent,
        1 => {
            let elem = *shared_transcripts.iter().next().unwrap();
            MappingResultTranscript::SingleTranscript(elem)
        }
        _ => MappingResultTranscript::Multimapped(shared_transcripts),
    }
}



fn build_ec2transcript(
    ec_dict: &HashMap<EC, Vec<u32>>,
    transcript_dict: &HashMap<u32, String>,
) -> HashMap<EC, HashSet<Transcriptname>> {
    let mut ec2transcript: HashMap<EC, HashSet<Transcriptname>> = HashMap::new();

    for (ec, transcript_ints) in ec_dict.iter() {
        let mut transcripts: HashSet<Transcriptname> = HashSet::new();

        for t_int in transcript_ints {
            let t_name = transcript_dict.get(t_int).unwrap();
            transcripts.insert(Transcriptname(t_name.clone()));
        }
        ec2transcript.insert(*ec, transcripts);
    }
    // // sanity check, make sure no Ec set is empty
    // for (ec, gset) in ec2gene.iter(){
    //     if gset.is_empty(){
    //         println!("{ec:?}'s geneset is empty");
    //     }
    // }
    ec2transcript
}

pub (crate) fn make_mapper_transcript(busfolder: &BusFolder) -> Ec2TranscriptMapper{
    let e2g = build_ec2transcript(
        &busfolder.parse_ecmatrix(), 
        &busfolder.parse_transcript()
    );
    Ec2TranscriptMapper::new(e2g)
}