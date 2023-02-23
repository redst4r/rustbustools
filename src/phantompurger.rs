use core::panic;
use std::{collections::{HashMap, HashSet}, hash::Hash, fs::File, io::Write};
use crate::{disjoint::{DisjointSubsets, SEPARATOR}, io::BusFolder};
use indicatif::{ProgressBar, ProgressStyle};
use crate::{bus_multi::CellUmiIteratorMulti, io::{BusRecord}, iterators::CbUmiIterator, utils::get_progressbar, consistent_genes::Ec2GeneMapper};

pub fn detect_overlap(busfolders: HashMap<String, String>) -> HashMap<Vec<String>, usize> {
    // deprecated
    // mesures the number of ovelapping CUGs across the experiments

    let mut total = 0;
    // TODO: this doesnt check if the EC overlaps
    for v in busfolders.values(){
        let cbumi_iter_tmp = CbUmiIterator::new(v);
        println!("determine size of iterator");
        // let now = Instant::now();
        let total_records = cbumi_iter_tmp.count();
        if total< total_records{
            total=total_records
        }
    }
    println!("total records {}", total);


    let multi_iter = CellUmiIteratorMulti::new(&busfolders);
    
    let bar = get_progressbar(total as u64);
    
    // let mut counter: HashMap<HashSet<String>, usize> = HashMap::new();
    let mut counter: HashMap<Vec<String>, usize> = HashMap::new();


    for (i,((_cb, _umi), record_dict)) in multi_iter.enumerate(){
        // let the_Set: HashSet<String> = record_dict.keys().cloned().collect();
        let mut the_set: Vec<String> = record_dict.keys().cloned().collect();
        the_set.sort();
        let val = counter.entry(the_set).or_insert(0);
        *val += 1; 

        if i % 1000000== 0{
            bar.inc(1000000);
        }
    }

    counter

}

#[derive(Debug)] 
pub struct FingerprintHistogram{
    order: Vec<String>,
    histogram: HashMap<Vec<u32>, usize>,
    drop_multimapped: bool
}

impl FingerprintHistogram{

    pub fn new(sample_order: &[String]) -> Self{
        let hist = HashMap::new();
        FingerprintHistogram{ 
            order: sample_order.to_owned(), // not sure why, linter is suggesting it
            histogram : hist,
            drop_multimapped: true
        }
    }

    pub fn add(&mut self, record_dict: HashMap<String, Vec<BusRecord>>, ecmapper_dict:&HashMap<String, &Ec2GeneMapper>){

        if self.drop_multimapped{
            // filter out the cases where a CB/UMI has more than one EC
            let filtered_dict: HashMap<String, BusRecord> = (record_dict).into_iter()
                .filter(|(_k,v)| v.len()==1)
                .map(|(k, v)| (k, v[0].clone())) // clone -> pop since we dont use vec after; not critical though
                .collect();

            for rdict in groupby_gene_simple(filtered_dict, ecmapper_dict){
                let fp_hash = make_fingerprint_simple(&rdict);
                
                // turn the hashmap into a vector, sorted acc to order
                let fp:  Vec<_>  = self.order.iter()
                    .map(|s| fp_hash.get(s).unwrap_or(&0)).cloned().collect();

                // update frequency
                let v = self.histogram.entry(fp).or_insert(0);
                *v += 1;
     
            }
        }
        else{
            panic!("not supported")
        }

        // let fp_hash = make_fingerprint(&record_dict);
        // // turn the hashmap into a vector, sorted acc to order
        // let fp:  Vec<_>  = self.order.iter()
        //     .map(|s| fp_hash.get(s).unwrap_or(&0)).cloned().collect();

        // // update frequency
        // let v = self.histogram.entry(fp).or_insert(0);
        // *v+=1;
    }

    pub fn to_csv(&self, outfile: &str){
        let mut fh = File::create(outfile).unwrap();

        let mut header = self.order.join(",");
        header.push_str(", frequency");

        writeln!(fh, "{}", header).unwrap();

        for (fingerprint, freq) in self.histogram.iter(){
            // concat with commas
            let mut s = fingerprint.iter().map(|i| i.to_string()).collect::<Vec<String>>().join(",");
            s.push_str(&format!(", {}", freq));
            writeln!(fh, "{}", s).unwrap();
        }
    }

}

pub fn make_fingerprint_histogram(busfolders: HashMap<String, BusFolder>) -> FingerprintHistogram{
    // main function here: takes a dict of busfolders, creates fingerprints for each molecule 
    // returns the fingerprints histogram (how often each fingerprint was observed)
    // and the ordering of the fingerprint (i.e. which element corresponds to which experiment)

    // create the EC2gene mappers
    let ecmapper_dict = busfolders.iter()
        .map(|(samplename, bfolder)|
            (samplename.clone(), &bfolder.ec2gene)   //#todo remove clone
        ).collect();

    // a list of busfile-names for the iterator
    let busnames =busfolders.iter()
        .map(|(s, bfolder)| (s.clone(), bfolder.get_busfile()))
        .collect();

    _make_fingerprint_histogram(&busnames, &ecmapper_dict)

}

pub fn _make_fingerprint_histogram(busnames: &HashMap<String, String>, ecmapper_dict: &HashMap<String, &Ec2GeneMapper>) -> FingerprintHistogram{

    // the actual workhorse, make_fingerprint_histogram is just a convenient wrapper
    let multi_iter = CellUmiIteratorMulti::new(busnames);
    
    let bar = ProgressBar::new_spinner();
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise}] {bar:40.cyan/blue} {pos} {per_sec}")
        .progress_chars("##-"));
    // let bar = get_progressbar(total as u64);
        
    let mut order: Vec<_> = busnames.keys().cloned().collect();
    order.sort();

    let mut fp_histo = FingerprintHistogram::new(&order);

    for (i,((_cb, _umi), record_dict)) in multi_iter.enumerate(){
        // println!("adding {:?}", record_dict);
        fp_histo.add(record_dict, ecmapper_dict);
        // println!("fphist {:?}", fp_histo);

        if i % 1000000== 0{
            bar.inc(1000000);
        }
    }
    fp_histo
}

type Fingerprint = HashMap<String, u32>;
fn make_fingerprint(record_dict: &HashMap<String,Vec<BusRecord>>) -> Fingerprint{
    // for a molecule detected over several experiments, 
    // get is fingerprint, i.e. freq across the experiments
    // if there's multiple ECs supporting the read, jsut sum their counts
    let fingerprint: Fingerprint = record_dict
        .iter()
        .map(|(k, v)| {
            (k.clone(), v.iter().map(|r|r.COUNT).sum::<u32>())
        }).collect();
    fingerprint
}

fn make_fingerprint_simple(record_dict: &HashMap<String,BusRecord>) -> Fingerprint{
    // for a molecule detected over several experiments, 
    // get is fingerprint, i.e. freq across the experiments
    let fingerprint: Fingerprint = record_dict
        .iter()
        .map(|(k, v)| {
            (k.clone(), v.COUNT)
        }).collect();
    fingerprint
}

#[derive(Eq, PartialEq, Debug, Hash)]
struct BusTmp {cb: u64, umi: u64, ec: u32, count:u32, flag:u32, samplename:String }
impl BusTmp {
    fn parse(item: &str) -> Option<Self> {

        let s: Vec<_> = item.split("@@").collect();
        if s.len() == 6{
            let cb:u64 = s[0].parse().unwrap();
            let umi:u64 = s[1].parse().unwrap();
            let ec: u32 = s[2].parse().unwrap();
            let count: u32 = s[3].parse().unwrap();
            let flag:u32 = s[4].parse().unwrap();
            let samplename = s[5].to_string();
            Some(BusTmp {cb , umi , ec, count, flag, samplename})
        }
        else{
            None
        }
        // BusTmp { value: item }
    }
    fn to_string(&self) -> String{
        let mut s = String::new();
        s.push_str(&self.cb.to_string());
        s.push_str("@@");

        s.push_str(&self.umi.to_string());
        s.push_str("@@");

        s.push_str(&self.ec.to_string());
        s.push_str("@@");

        s.push_str(&self.count.to_string());
        s.push_str("@@");

        s.push_str(&self.flag.to_string());
        s.push_str("@@");

        s.push_str(&self.samplename.to_string());
        s
    }
}

// #[inline(never)]
pub fn groupby_gene_simple(record_dict: HashMap<String,BusRecord>, ecmapper_dict: &HashMap<String, &Ec2GeneMapper>) -> Vec<HashMap<String, BusRecord>>
{
    // assuming no multimapped reads, hence record dict contains only single entries (not a list)

    // check if the genes are consistent across samples
    // if so yield the record as is
    // otherwise split into multiple records
    if record_dict.len() == 1{
        return vec![record_dict];
        // return
    };
    
    let mut records:Vec<BusTmp> = Vec::new();
    for (sname, r) in record_dict{
        records.push(
                BusTmp { cb: r.CB, umi: r.UMI, ec: r.EC, count: r.COUNT, flag: r.FLAG, samplename: sname.clone() }
            )
    }

    let mut genes: HashMap<BusTmp, &HashSet<u32>> = HashMap::new();
    for b in records{
        let ecmapper = ecmapper_dict.get(&b.samplename).unwrap();
        let g = ecmapper.get_genes(b.ec);
        genes.insert(b, g);
    }

    let mut disjoint_set = DisjointSubsets::new();
    for (b, gset) in genes.drain(){
        disjoint_set.add(b.to_string(), gset.clone());
    }

    let mut emit_vector: Vec<HashMap<String, BusRecord>> = Vec::new();
    for bus_string_concat in disjoint_set.disjoint_sets.keys(){
        // these are now all records across the samples that map to the same gene
        let bus_tuple: Vec<BusTmp> = bus_string_concat.split(SEPARATOR).map(|bstring|BusTmp::parse(bstring).unwrap()).collect();

        let mut emited_dict: HashMap<String, BusRecord> = HashMap::new();
        for b in bus_tuple{

            if emited_dict.contains_key(&b.samplename){
                panic!("cant happen, each sample only has one record")
            } 
            else{
                let brecord = BusRecord{ CB: b.cb, UMI: b.umi, EC:0, COUNT: b.count, FLAG: b.flag};
                emited_dict.insert(b.samplename, brecord);
            }
        }
        emit_vector.push(emited_dict);
    }
    emit_vector
}


fn groupby_gene(record_dict: HashMap<String,Vec<BusRecord>>, ecmapper_dict: &HashMap<String, Ec2GeneMapper>) 
    -> Vec<HashMap<String, BusRecord>>
    {
    // multiple experiments yielded the same CB/UMI
    // turn this into a list of record_dicts, where the gene is consistent

    let mut records:Vec<BusTmp> = Vec::new();
    for (sname, recordlist) in record_dict{
        for r in recordlist{
            records.push(
                BusTmp { cb: r.CB, umi: r.UMI, ec: r.EC, count: r.COUNT, flag: r.FLAG, samplename: sname.clone() }
            )
        }
    }
    // if records.len() == 1{
    //     return record_dict
    //     // return
    // };

    // get genes for each record.
    let mut genes: HashMap<BusTmp, &HashSet<u32>> = HashMap::new();
    for b in records{
        let ecmapper = ecmapper_dict.get(&b.samplename).unwrap();
        let g = ecmapper.get_genes(b.ec);
        genes.insert(b, g);
    }

    // also build the disjoint set based on the genes
    let mut disjoint_set = DisjointSubsets::new();
    for (b, gset) in genes.drain(){
        disjoint_set.add(b.to_string(), gset.clone());
    }

    let mut emit_vector: Vec<HashMap<String, BusRecord>> = Vec::new();
    for bus_string_concat in disjoint_set.disjoint_sets.keys(){
        // these are now all records across the samples that map to the same gene
        let bus_tuple: Vec<BusTmp> = bus_string_concat.split('_').map(|bstring|BusTmp::parse(bstring).unwrap()).collect();

        let mut emited_dict: HashMap<String, BusRecord> = HashMap::new();
        for b in bus_tuple{

            // if we already have a record in that sample, just update the count
            if emited_dict.contains_key(&b.samplename){
                let brecord = emited_dict.get_mut(&b.samplename).unwrap();
                brecord.COUNT += b.count;
            } 
            else{
                let brecord = BusRecord{ CB: b.cb, UMI: b.umi, EC:0, COUNT: b.count, FLAG: b.flag};
                emited_dict.insert(b.samplename, brecord);
            }
        }
        emit_vector.push(emited_dict);
    }
    emit_vector
}

#[cfg(test)]
mod tests{
    use std::collections::{HashSet, HashMap};
    use crate::{consistent_genes::Ec2GeneMapper, io::{BusRecord, BusFolder}, phantompurger::{groupby_gene_simple, _make_fingerprint_histogram, make_fingerprint_simple}};
    use super::make_fingerprint_histogram;

    fn create_dummy_ec() ->Ec2GeneMapper{
        let ec0: HashSet<String> = vec!["A".to_string()].into_iter().collect();
        let ec1: HashSet<String> = vec!["B".to_string()].into_iter().collect();
        let ec2: HashSet<String> = vec!["A".to_string(), "B".to_string()].into_iter().collect();
        let ec3: HashSet<String> = vec!["C".to_string(), "D".to_string()].into_iter().collect();

        let ec_dict: HashMap<u32, HashSet<String>> = HashMap::from([
            (0, ec0), // A
            (1, ec1), // B
            (2, ec2), // A,B
            (3, ec3), // C,D
            ]);

        Ec2GeneMapper::new(ec_dict)
    }


    #[test]
    fn test_groupby_gene_simple(){

        // --------------------------------------------
        // two records that share the same EC, should be grouped into a single emission
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0};

        let es1 = create_dummy_ec();
        let es2 = create_dummy_ec();
        let es3 = create_dummy_ec();
        let es_dict: HashMap<String, &Ec2GeneMapper> = vec![
            ("s1".to_string(), &es1),
            ("s2".to_string(), &es2),
            ("s3".to_string(), &es3)
            ]
        .into_iter().collect();


        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
        ].into_iter().collect();
        let res = groupby_gene_simple(record_dict, &es_dict);
        println!("{:?}", res);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].len(), 2);

        // --------------------------------------------
        // two records, different EC, but same gene, should be grouped into a single emission
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 1, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 1, FLAG: 0};
        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
        ].into_iter().collect();
        let res = groupby_gene_simple(record_dict, &es_dict);
        println!("{:?}", res);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].len(), 2);

        // --------------------------------------------
        // two records, different EC, and inconsistnet gene, should be grouped into a TWO emissions
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 1, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 1, FLAG: 0};
        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
        ].into_iter().collect();
        let res = groupby_gene_simple(record_dict, &es_dict);
        println!("{:?}", res);
        assert_eq!(res.len(), 2);
        assert_eq!(res[0].len(), 1);
        assert_eq!(res[1].len(), 1);

        // --------------------------------------------
        // three records, A, B, (A,B). should be yielded as a single emssion
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 1, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 1, FLAG: 0};
        let t1 = BusRecord{CB: 0, UMI: 1, EC: 2, COUNT: 1, FLAG: 0};
        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
            ("s3".to_string(), t1),
        ].into_iter().collect();
        let res = groupby_gene_simple(record_dict, &es_dict);
        println!("{:?}", res);
        assert_eq!(res.len(), 1);
        assert_eq!(res[0].len(), 3);

    }

    #[test]
    fn test_make_fingerprint_histogram(){
        use crate::io::setup_busfile;

        // a pair, same EC
        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0};
        
        // singleton in sample1
        let r2 =BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0}; 
        
        // singleton in sample2
        let s2 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 10, FLAG: 0};


        // pair but only gene overlap
        let r3 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
        let s3 = BusRecord{CB: 1, UMI: 3, EC: 2, COUNT:  3, FLAG: 0}; 

        //overall we should get 
        // [12, 0] = 1
        // [0, 10] = 1
        // [2, 3] = 2

        let v1 = vec![r1.clone(),r2.clone(),r3.clone()];
        let v2 = vec![s1.clone(),s2.clone(),s3.clone()];

        // write the records to file
        let (busname1, _dir1) =setup_busfile(&v1);
        let (busname2, _dir2) =setup_busfile(&v2);

        let hashmap = HashMap::from([
            ("s1".to_string(), busname1.to_string()),
            ("s2".to_string(), busname2.to_string())
        ]);

        let es1 = create_dummy_ec();
        let es2 = create_dummy_ec();
        let es_dict: HashMap<String, &Ec2GeneMapper> = vec![
            ("s1".to_string(), &es1),
            ("s2".to_string(), &es2)]
            .into_iter().collect();

        let s = _make_fingerprint_histogram(&hashmap, &es_dict);
        println!("{:?}", s);

        let e = *s.histogram.get(&vec![12, 0]).unwrap();
        assert_eq!(e, 1);

        let e = *s.histogram.get(&vec![0, 10]).unwrap();
        assert_eq!(e, 1);

        let e = *s.histogram.get(&vec![2, 3]).unwrap();
        assert_eq!(e, 2);

        s.to_csv("/tmp/finger.csv")

        // let s = detect_overlap(hashmap);
    }

    #[test]

    fn test_make_fingerprint_simple(){

        let r1 =BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 2, FLAG: 0};
        let s1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 3, FLAG: 0};

        let record_dict = vec![
            ("s1".to_string(), r1),
            ("s2".to_string(), s1),
        ].into_iter().collect();
        let res = make_fingerprint_simple(&record_dict);
        println!("{:?}", res);

        let exp: HashMap<_,_> = vec![("s1".to_string(), 2), ("s2".to_string(), 3)].into_iter().collect();
        assert_eq!(res, exp);
    }


    // #[test]
    pub fn testing2(){
        
        let t2g = "/home/michi/bus_testing/transcripts_to_genes.txt";
        // let b1 = BusFolder::new("/home/michi/bus_testing/bus_output/", t2g);
        let b1 = BusFolder::new("/home/michi/bus_testing/bus_output_short/", t2g);
        let b2 = BusFolder::new("/home/michi/bus_testing/bus_output_short/", t2g);

        // let hashmap = HashMap::from([
        //     ("full".to_string(), "/home/michi/bus_testing/bus_output/output.corrected.sort.bus".to_string()),
        //     ("short".to_string(), "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus".to_string())
        // ]);
        let hashmap = HashMap::from([
            ("full".to_string(), b1),
            ("short".to_string(), b2)
        ]);

        let s = make_fingerprint_histogram(hashmap);

        s.to_csv("/tmp/testing2.csv");
        println!("{:?}", s)
    }
}