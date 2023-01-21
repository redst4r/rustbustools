
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::Write;
use crate::io::{BusFolder,CellIterator, BusRecord, group_record_by_cb_umi, BusIteratorBuffered};
use sprs;
use indicatif::{ProgressBar, ProgressStyle};
// use std::iter::FromIterator;
use statrs::distribution::{Multinomial, Binomial};  // warning the statrs::Binomial has very slow sampling (sum of Bernullis)
use probability;  // use probability::distribution::Binomial instead, which does inverse cdf sampling
use probability::distribution::Sample;
use probability::prelude::*;


use rand;
use rand::distributions::Distribution;


pub fn count_bayesian(bfolder: BusFolder) {

    let bfile = format!("{}/{}", bfolder.foldername, bfolder.busfile);
    println!("{}",bfile);

    // let ec2gene = bfolder.ec2gene;
    
    // // the total amount of reads in the busfile
    // println!("Getting total reads");
    // let busiter = BusIteratorBuffered::new(&bfile);
    // let total_reads: u32 = busiter.map(|record| record.COUNT).sum();

    // // get the fraction needs another iterator
    // // TODO add the 0.5 pseudo
    // println!("Getting thetas");
    // let busiter = BusIteratorBuffered::new(&bfile);
    // let thetas: Vec<f64> = busiter.map(|record| record.COUNT as f64 / total_reads as f64 ).collect();

    // let norm_constant: f64 =  thetas.iter().sum();

    // println!("{}, {}", total_reads,norm_constant);

    let busiter = BusIteratorBuffered::new(&bfile);
    println!("Getting counts");
    let counts: Vec<f64> = busiter.map(|record| record.COUNT as f64).collect();
    let total_counts:f64 = counts.iter().sum();

    println!("counts total {}, Entries {}",total_counts as u64 ,counts.len());
    println!("summing");


    // let total: f64 = counts.iter().sum();
    // let mut r = rand::thread_rng();
    // println!("RV");
    // let n = Multinomial::new(&counts, total as u64).unwrap().sample(&mut r);
    // println!("{:?}", n[0]);

}


pub fn multinomial_sample_statrs(n: u64, pvec: Vec<f64>) -> Vec<f64>{
    let mut r = rand::thread_rng();
    let mut x :Vec<f64> = Vec::new();

    // normalize the pvec
    let _sum: f64 = pvec.iter().sum();
    let pvec_norm: Vec<f64> = pvec.iter().map(|x| x/_sum).collect();
    let dim = pvec_norm.len();

    let mut remaining_p = 1.0;
    let mut remaining_n = n;
    for (counter, p) in pvec_norm.iter().enumerate(){
        // binomial with BIN(p/remaining_p, remaining_n)
        if remaining_n == 0{
            x.push(0.0);
        }
        else if counter == dim - 1{
            //last element, p will be 1 (or due to errors, a litte >1)
            x.push(remaining_n as f64);
        }
        else{
            let _ptmp = p / remaining_p;
            if !(0.0 < _ptmp && _ptmp < 1.0){
                println!("{:?}", &pvec_norm[counter..]);
                panic!("0<{}<1, counter {}", _ptmp, counter);
            }            
            let b = Binomial::new(_ptmp, remaining_n).expect(&format!("p={} remaining_p={} ptmp={} n={}", p, remaining_p ,_ptmp, remaining_n)).sample(&mut r);

            x.push(b);
            remaining_n -= b as u64;
            remaining_p -= p;
        }
    }
    x
}

pub fn multinomial_sample(n: u64, pvec: Vec<f64>) -> Vec<f64>{
    /*
    my own multinomial sampling, using the fact that all marginals are binomial

    statrs version does the same algorithm, but relies internally on a statrs::distribution::Binomial
    which is extremely slow.
    */
    let mut source = source::default(42);

    let mut x :Vec<f64> = Vec::new();

    // normalize the pvec
    let _sum: f64 = pvec.iter().sum();
    let pvec_norm: Vec<f64> = pvec.iter().map(|x| x/_sum).collect();
    let dim = pvec_norm.len();

    let mut remaining_p = 1.0;
    let mut remaining_n = n;
    for (counter, p) in pvec_norm.iter().enumerate(){
        // binomial with BIN(p/remaining_p, remaining_n)
        if remaining_n == 0{
            x.push(0.0);
        }
        else if counter == dim - 1{
            //lastt else if element, p will be 1 (or due to errors, a litte >1)
            x.push(remaining_n as f64);
        }
        else{
            let _ptmp = p / remaining_p;

            if !(0.0 < _ptmp && _ptmp < 1.0){
                println!("{:?}", &pvec_norm[counter..]);
                panic!("0<{}<1, counter {}", _ptmp, counter);
            }
            let btmp = probability::distribution::Binomial::new(remaining_n as usize, _ptmp ).sample(&mut source); //.expect(&format!("p={} remaining_p={} ptmp={} n={}", p, remaining_p ,_ptmp, remaining_n)).sample(&mut r);
            let b = btmp as f64;


            x.push(b);
            remaining_n -= b as u64;
            remaining_p -= p;
        }
    }
    x
}

pub fn count(bfolder: BusFolder) -> sprs::CsMat<usize>{
    /*
    busfile to count matrix, analogous to "bustools count"
    */
    let bfile = format!("{}/{}", bfolder.foldername, bfolder.busfile);
    println!("{}",bfile);

    let ec2gene = bfolder.ec2gene;
    let cb_iter = CellIterator::new(&bfile);

    let mut all_expression_vector: HashMap<u64, HashMap<String, u32>> = HashMap::new();

    // let bar = ProgressBar::new_spinner();
    let bar = ProgressBar::new(1_000_000);
    bar.set_style(ProgressStyle::default_bar()
        .template("[{elapsed_precise} ETA {eta}] {bar:40.cyan/blue} {pos}/{len} {per_sec}")
        .progress_chars("##-"));

    for (cb, record_list) in cb_iter{
        // println!("{} {}",cb, record_list.len());
        let s = records_to_expression_vector(record_list, &ec2gene);

        // this will also insert emtpy cells (i.e. their records are all multimapped)
        all_expression_vector.insert(cb, s);
        bar.inc(1)
    }

    //collect all genes
    let mut genelist = HashSet::new(); 
    for glist in ec2gene.values(){
        genelist.extend(glist)
    }
    let genelist_vector :Vec<&String>= genelist.into_iter().collect::<Vec<&String>>();

    let countmatrix = expression_vectors_to_matrix(all_expression_vector, genelist_vector );
    println!("{:?} nnz {}", countmatrix.shape(), countmatrix.nnz());
    countmatrix

}


fn find_consistent_genes(record_list: Vec<BusRecord>, ec2gene: &HashMap<u32, Vec<String>>) -> HashSet<String>{
    // is there a single gene consistent with all those busrecords?

    let mut genes :Vec<HashSet<String>>= Vec::new();
    // let mut genes :HashSet<String>= HashSet::new();
    for r in record_list {
        let ec_genes = &ec2gene[&r.EC];
        let mut ec_genes_set = HashSet::new();
        for g in ec_genes.into_iter(){  // todo: not sure how to do this nicely, just transform vec into set
            ec_genes_set.insert(g.clone());
        }
        genes.push(ec_genes_set);
    
        // if  genes.len() == 0{
        //     // if nothing is in the set, add
        //     genes = ec_genes_set;
        // }
        // else{
        //     genes = genes.intersection(&ec_genes_set).cloned().collect();
        // }
        // now, genes should contain only those genes that are consisten with all records of that CB/UMI
    }
    let consistent_genes = intersect_sets(&mut genes);
    consistent_genes
}


fn records_to_expression_vector_groupby(record_list: Vec<BusRecord>, ec2gene: &HashMap<u32, Vec<String>>) -> HashMap<String, u32>{
    let mut expression_vector: HashMap<String, u32> = HashMap::new(); // gene -> count
    let mut _multimapped = 0;
    let mut _inconsistant= 0;

    // first, group the records by UMI
    // TODO: EXPENSIVE!! 25k/s

    use itertools::Itertools;
    let cb_umi_grouped = record_list.into_iter().group_by(|a|(a.CB, a.UMI));

    for (_key, group) in cb_umi_grouped.into_iter(){
        let records: Vec<BusRecord> = group.collect();
        // all records coresponding to the same UMI
        // TODO: EXPENSIVE!! 25k/s
        let consistent_genes= find_consistent_genes(records, ec2gene);

        // let consistent_genes = ec2gene.get(&records[0].EC).unwrap();
        // let consistent_genes = vec![records[0].EC];


        if consistent_genes.len() > 1{
            //multimapped
            _multimapped += 1;
        }
        else if consistent_genes.len() == 1 {
            //single gene
            // let g = consistent_genes.drain().next().unwrap();

            // let v: Vec<String> = consistent_genes.into_iter().collect();  // Set to Vec
            // let g = &v[0];

            let g = consistent_genes.iter().next().unwrap();  // Set to first element


            let val = expression_vector.entry(g.to_string()).or_insert(0);
            *val += 1;        
        }
        else{
            // inconsistant
            _inconsistant += 1
        }
    }
    // println!("{}, {}",_multimapped, _inconsistant);

    expression_vector

}


fn records_to_expression_vector(record_list: Vec<BusRecord>, ec2gene: &HashMap<u32, Vec<String>>) -> HashMap<String, u32>{
    /*
    turn the list of records of a single CB into a expression vector: per gene, how many umis are observed
    TODO this doesnt consider multiple records with same umi/cb, but EC mapping to different genes
    */
    let mut expression_vector: HashMap<String, u32> = HashMap::new(); // gene -> count
    let mut _multimapped = 0;
    let mut _inconsistant= 0;

    // first, group the records by UMI
    // TODO: EXPENSIVE!! 25k/s
    let cb_umi_grouped = group_record_by_cb_umi(record_list);

    // let records = record_list;
    for ((_cb, _umi), records) in cb_umi_grouped{
        // all records coresponding to the same UMI
        // TODO: EXPENSIVE!! 25k/s
        let consistent_genes= find_consistent_genes(records, ec2gene);

        // let consistent_genes = ec2gene.get(&records[0].EC).unwrap();
        // let consistent_genes = vec![records[0].EC];


        if consistent_genes.len() > 1{
            //multimapped
            _multimapped += 1;
        }
        else if consistent_genes.len() == 1 {
            //single gene
            // let g = consistent_genes.drain().next().unwrap();

            // let v: Vec<String> = consistent_genes.into_iter().collect();  // Set to Vec
            // let g = &v[0];

            let g = consistent_genes.iter().next().unwrap();  // Set to first element


            let val = expression_vector.entry(g.to_string()).or_insert(0);
            *val += 1;        
        }
        else{
            // inconsistant
            _inconsistant += 1
        }
    }
    // println!("{}, {}",_multimapped, _inconsistant);

    expression_vector
}

fn expression_vectors_to_matrix(all_expression_vector: HashMap<u64, HashMap<String, u32>>, genelist: Vec<&String>) -> sprs::CsMat<usize>{
    /*
    turn an collection of expression vector (from many cells)
    into a sparse count matrix
    */

    // sparse matrix indices
    let mut ii: Vec<usize> = Vec::new();
    let mut jj: Vec<usize> = Vec::new();
    let mut vv: Vec<usize> = Vec::new();

    // the cell barcodes, same order as in the matrix
    let mut cbs: Vec<u64> = Vec::new(); 

    // mapping for gene-> index/order
    let mut gene2index: HashMap<&String, usize> = HashMap::new();
    for (i, g) in genelist.iter().enumerate(){
        gene2index.insert(g, i);
    }

    for (i, (cb, expr_vec)) in all_expression_vector.iter().enumerate(){
        for (gene, count) in expr_vec{
            ii.push(i);
            jj.push(gene2index[gene]);
            vv.push(*count as usize)
        }
        cbs.push(*cb)
    }
   
    let c: sprs::TriMat<usize> = sprs::TriMat::from_triplets(
        (cbs.len(), genelist.len()),
        ii,
        jj,
        vv
    );
    let b: sprs::CsMat<_> = c.to_csr();
    b

}

fn intersect_sets(list_of_sets: &mut Vec<HashSet<String>>) -> HashSet<String>{
    // warning this mutates the list by poping
    // todo we probably can do this with iterator reduce
    let mut intersection_set = list_of_sets.pop().unwrap();

    for s in list_of_sets{
        intersection_set = intersection_set.intersection(&s).map(|x|x.to_string()).collect();  //todo ugly conversion to string
    }
    intersection_set


    // let t= list_of_sets.iter()
                                                //  .reduce(|acc, e| &acc.intersection(e).map(|x|x.to_string()).collect());

    // let a:<Vec<String>> = list_of_sets[0].collect();

    // return t.unwrap();
}

// #[test]
pub fn test_multinomial(dim: i32){
    let n = 1000;
    let p = (1..dim).map(|x| x as f64).collect();
    // println!("{:?}", p);

    let x = multinomial_sample(n, p);
    let n2: f64 = x.iter().sum();
    // println!("{:?} {}", x, n2);
}
// #[test]
pub fn test_multinomial_stats(dim: i32){
    let n = 1000;
    let p: Vec<f64> = (1..dim).map(|x| x as f64).collect();
    // println!("{:?}", *p);

    let mut r = rand::thread_rng();
    let n = Multinomial::new(&p, n).unwrap();
    let x = n.sample(&mut r);
    let n2: f64 = x.iter().sum();
    // println!("{:?} {}", x, n2);

}


#[test]    
fn test_itt(){  
    // let t2g_file = String::from("/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt");
    // let foldername = String::from("/home/michi/mounts/TB4drive/ISB_data/MNGZ01/MS_processed/S1/kallisto/sort_bus/bus_output");
    let t2g_file = String::from("/home/michi/bus_testing/transcripts_to_genes.txt");
    let foldername = String::from("/home/michi/bus_testing/bus_output");


    let b = BusFolder::new(foldername, t2g_file);
    let count_matrix = count(b);

    write_sprs_to_file(count_matrix, "/tmp/test.mtx");
    // count_bayesian(b)
}

#[test]
// use crate::count::records_to_expression_vector;
fn test_records_to_expression_vector(){

    let ec_dict = HashMap::from([
        (0, vec!["G1".to_string(), "G2".to_string()]),
        (1, vec!["G1".to_string()]),
        (2, vec!["G2".to_string()]),
        (3, vec!["G1".to_string(), "G2".to_string()]),
    ]);

    // those three records are consistent with G1
    let r1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0};
    let r2 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

    // those records are consistent with G1 and G2, hence multimapped
    let r4 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 2, FLAG: 0};
    let r5 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 2, FLAG: 0};
    let r6 = BusRecord{CB: 0, UMI: 2, EC: 3, COUNT: 2, FLAG: 0};

    // those records are inconsistent with G1 vs G2
    let r7= BusRecord{CB: 0, UMI: 3, EC: 0, COUNT: 2, FLAG: 0};
    let r8= BusRecord{CB: 0, UMI: 3, EC: 1, COUNT: 2, FLAG: 0};
    let r9= BusRecord{CB: 0, UMI: 3, EC: 2, COUNT: 2, FLAG: 0};

    // those records are consistent with G1
    let r10= BusRecord{CB: 0, UMI: 5, EC: 1, COUNT: 2, FLAG: 0};
    let r11= BusRecord{CB: 0, UMI: 5, EC: 0, COUNT: 2, FLAG: 0};

    // those records are consistent with G2
    let r12 = BusRecord{CB: 0, UMI: 4, EC: 2, COUNT: 2, FLAG: 0};
    let r13= BusRecord{CB: 0, UMI: 4, EC: 0, COUNT: 2, FLAG: 0};

    let records0 = vec![r1,r2];
    let c0 = records_to_expression_vector(records0, &ec_dict);
    assert_eq!(c0, HashMap::from([("G1".to_string(), 1)]));

    let records1 = vec![r1,r2, r10, r11];
    let c1= records_to_expression_vector(records1, &ec_dict);
    assert_eq!(c1, HashMap::from([("G1".to_string(), 2)]));


    let records2 = vec![r4,r5,r6];
    let c2= records_to_expression_vector(records2, &ec_dict);
    assert_eq!(c2, HashMap::from([]));


    let records3 = vec![r1,r2, r4,r5,r6,r7,r8,r9,r10, r11,r12,r13];
    let c3 = records_to_expression_vector(records3, &ec_dict);
    assert_eq!(c3, HashMap::from([("G1".to_string(), 2), ("G2".to_string(), 1)]));



}


