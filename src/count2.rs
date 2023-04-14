use std::collections::{HashMap, HashSet, BTreeSet};
use std::time::Instant;
use sprs::{DenseVector};
use crate::consistent_genes::{find_consistent, Ec2GeneMapper, GeneId, CB};
use crate::countmatrix::CountMatrix;
use crate::iterators::CbUmiGroupIterator;
use crate::io::{BusFolder, BusRecord};
use crate::multinomial::multinomial_sample;
use crate::utils::{get_progressbar, int_to_seq};

fn countmap_to_matrix(countmap: &HashMap<(CB, GeneId), usize>, gene_vector: Vec<String>) -> CountMatrix{

    // get all CBs, a BTreeSet gives us order for free
    // let cb_set: BTreeSet<u64> = BTreeSet::new();
    println!("getting all CBs");
    let all_cbs = countmap.keys().map(|(cb, _gene)| cb).collect::<BTreeSet<_>>();
    // println!("getting all genes");
    // let all_genes = countmap.keys().map(|(cb, gene)| *gene.clone()).collect::<BTreeSet<_>>();

    println!("building index");
    // some issues with the cb.clone: clippy complains!
    // let cb_ix = all_cbs.iter().enumerate().map(|(ix, cb)|(cb.clone(), ix)).collect::<HashMap<_,_>>();
    let cb_ix = all_cbs.iter().enumerate().map(|(ix, cb)|(**cb, ix)).collect::<HashMap<_,_>>();

    // sparse matrix indices
    let mut ii: Vec<usize> = Vec::new();
    let mut jj: Vec<usize> = Vec::new();
    let mut vv: Vec<usize> = Vec::new();

    for ((cb, geneid), counter) in countmap{

        let cbi = cb_ix.get(cb).unwrap();
        let genei = (*geneid).0 as usize;
        ii.push(*cbi);
        jj.push(genei);
        vv.push(*counter);
    }

    let c: sprs::TriMat<usize> = sprs::TriMat::from_triplets(
        (cb_ix.len(), gene_vector.len()),
        ii,
        jj,
        vv
    );

    let b: sprs::CsMat<_> = c.to_csr();

    let cbs_seq: Vec<String> = all_cbs.into_iter().map(|x|int_to_seq((*x).0, 16)).collect();
    // let gene_seq: Vec<String> = gene_vector.into_iter().map(|x|x.clone()).collect();
    let gene_seq: Vec<String> = gene_vector.into_iter().collect(); //not sure if this does anything
    
    CountMatrix::new(b, cbs_seq, gene_seq)

}

pub fn baysian_count(bfolder: BusFolder, ignore_multimapped:bool, n_samples: usize){
    let bfile = bfolder.get_busfile();
    println!("{}",bfile);

    println!("determine size of iterator");
    let cbumi_iter_tmp = bfolder.get_iterator().groupby_cbumi();
    let now = Instant::now();
    let total_records = cbumi_iter_tmp.count();

    let cbumi_iter_tmp = bfolder.get_iterator().groupby_cbumi();

    let max_length_records = cbumi_iter_tmp.map(|(_cbumi, rlist)| rlist.len()).max().unwrap();

    let elapsed_time = now.elapsed();
    println!("determined size of iterator {} in {:?}. Longest element: {} in a single CB/UMI", total_records, elapsed_time, max_length_records);

    // handles the mapping between EC and gene
    let egm = &bfolder.ec2gene;

    // prep for the multinomial sample
    println!("Preparing the probability vector for mutlinomial");
    let cbumi_iter_tmp = bfolder.get_iterator().groupby_cbumi();

    let count_vec:Vec<_> = cbumi_iter_tmp
        .flat_map(|(_cbumi, rlist)| rlist.into_iter().map(|r|r.COUNT as f64))
        .collect();

    let total_counts: f64 = count_vec.iter().sum();
    let p_vec: Vec<f64> = count_vec.into_iter().map(|c| c/total_counts).collect();
    println!("Done: {} rercods, {} counts", p_vec.len(), total_counts);

    use probability::prelude::*;
    let mut random_source = source::default(42);

    let mut counter = 0;
    for i in 0..n_samples{

        // CB,gene_id -> count
        let mut all_expression_vector: HashMap<(CB, GeneId), usize> = HashMap::new();
        let mut n_mapped = 0;
        let mut n_multi_inconsistent = 0;

        // subsample the count vector
        println!("Iteration {}: Mutlinomial sample", i);
        let new_count_sample = multinomial_sample(total_counts as u64, &p_vec, &mut random_source);
        println!("Done");

        let cbumi_iter = bfolder.get_iterator().groupby_cbumi();

        let now = Instant::now();
        let bar = get_progressbar(total_records as u64);
        let mut current_record_counter: usize = 0;

        for ((cb, _umi), rlist) in cbumi_iter{

            // inject the sampled numbers into the records

            let indices = current_record_counter..current_record_counter+rlist.len();
            let injected_counts: Vec<u32> = indices.map(|idx| *new_count_sample.index(idx) as u32 ).collect(); // wrning f64->u32
            // let mut injected_records: Vec<BusRecord> = Vec::with_capacity(rlist.len());
            let mut injected_records: Vec<BusRecord> = rlist.clone();

            for i in 0..injected_records.len(){
                // let mut r = injected_records.get_mut(i).expect(&format!("injected_records {}", i));
                let mut r = injected_records.get_mut(i).unwrap_or_else(|| panic!("injected_records {}", i));
                let c = injected_counts.get(i).unwrap_or_else(|| panic!("injected_counts {}", i));
                r.COUNT = *c;
            }

            injected_records.retain(|r|r.COUNT>0);

            // for (r, new_count) in injected_records.iter_mut().zip(injected_counts.into_iter()){
            //     r.COUNT = new_count;
            //     injected_records.push(r);
            // }
            current_record_counter += rlist.len();

            if injected_records.is_empty(){
                continue;
            }

            if let Some(g) = count_from_record_list(&injected_records, egm, ignore_multimapped){
                // the records could be made into a single count for gene g
                let key = (CB(cb), g);
                let current_count = all_expression_vector.entry(key).or_insert(0);
                *current_count+=1;
    
                n_mapped += 1;
            }
            else{
                // multimapped, or not consistently mapped
                n_multi_inconsistent += 1;
            }
    
            if counter % 1_000_000 == 0{
                bar.inc(1_000_000);
            }
            counter += 1;
        }

        let elapsed_time = now.elapsed();
        let fraction_mapped = n_multi_inconsistent as f64 /(n_mapped as f64+n_multi_inconsistent as f64);
        println!("Iteration {}: Mapped {}, multi-discard {} ({}%) in {:?}",i, n_mapped, n_multi_inconsistent, 100.0*fraction_mapped, elapsed_time);

        let genelist_vector: Vec<String> = egm.get_gene_list();
        // this is how genes are ordered as by EGM
        // i.e. countmap[cb, i] corresponds to the number of count of genelist_vector[i]
    
        let countmatrix = countmap_to_matrix(&all_expression_vector, genelist_vector );
        println!("{:?} nnz {}", countmatrix.matrix.shape(), countmatrix.matrix.nnz()); 
        println!("finished iteration {}", i)         
    }
}

fn count_from_record_list(records: &Vec<BusRecord>, egmapper: &Ec2GeneMapper, ignore_multimapped:bool) -> Option<GeneId>{

    // given a set of Records from the same CB/UMI, are they consistent with a particular gene
    // which would yield a signel count
    let consistent_genes: HashSet<GeneId>;

    if !ignore_multimapped{
            consistent_genes= find_consistent(records, egmapper);
    }else if records.len() > 1{
        return None
    }
    else{
        // consistent_genes = ec2gene.get(&records[0].EC).unwrap().clone();
        consistent_genes = egmapper.get_genes(records[0].EC).clone();
    }
       
    // uniquely mapped to a single gene
    if consistent_genes.len() == 1{
        let g = consistent_genes.into_iter().next().unwrap();
        Some(g)
    }
    else{
        // multimapped or inconsistent
        None
    }

}

pub fn count(bfolder: &BusFolder, ignore_multimapped:bool) -> CountMatrix {
    /*
    busfile to count matrix, analogous to "bustools count"
    */
    let bfile = bfolder.get_busfile();
    println!("{}",bfile);

    let cbumi_iter = bfolder.get_iterator().groupby_cbumi();
    let cbumi_iter_tmp = bfolder.get_iterator().groupby_cbumi();

    println!("determine size of iterator");
    let now = Instant::now();
    let total_records = cbumi_iter_tmp.count();

    let cbumi_iter_tmp = bfolder.get_iterator().groupby_cbumi();

    let max_length_records = cbumi_iter_tmp.map(|(_cbumi, rlist)| rlist.len()).max().unwrap();

    let elapsed_time = now.elapsed();

    println!("determined size of iterator {} in {:?}. Longest element: {} in a single CB/UMI", total_records, elapsed_time, max_length_records);

    // CB,gene_id -> count
    let mut all_expression_vector: HashMap<(CB, GeneId), usize> = HashMap::new();
    let bar = get_progressbar(total_records as u64);

    let mut n_mapped = 0;
    let mut n_multi_inconsistent = 0;

    let now = Instant::now();

    for (counter, ((cb, _umi), record_list)) in cbumi_iter.enumerate() {

        // try to map the records of this CB/UMI into a single gene
        if let Some(g) = count_from_record_list(&record_list, &bfolder.ec2gene, ignore_multimapped){
            // the records could be made into a single count for gene g
            let key = (CB(cb), g);
            let current_count = all_expression_vector.entry(key).or_insert(0);
            *current_count+=1;

            n_mapped += 1;
        }
        else{
            // multimapped, or not consistently mapped
            n_multi_inconsistent += 1;
            let cbstr = int_to_seq(cb, 16);
            let umistr = int_to_seq(_umi, 12);
            println!("not countable {cbstr}/{umistr} {:?}", record_list);
            let cgeneids= find_consistent(&record_list, &bfolder.ec2gene);
            let cgenes: Vec<_> = cgeneids.iter().map(|gid| bfolder.ec2gene.resolve_gene_id(*gid)).collect();
            println!("{cgenes:?}")

        }

        if counter % 1_000_000 == 0{
            bar.inc(1_000_000);
        }
    }

    let elapsed_time = now.elapsed();
    println!("Mapped {}, multi-discard {} in {:?}", n_mapped, n_multi_inconsistent, elapsed_time); 

    let genelist_vector: Vec<String> = bfolder.ec2gene.get_gene_list();

    // this is how genes are ordered as by EGM
    // i.e. countmap[cb, i] corresponds to the number of count of genelist_vector[i]

    let countmatrix = countmap_to_matrix(&all_expression_vector, genelist_vector );

    let (shape1,shape2) = countmatrix.matrix.shape();
    println!("{}x{} nnz {}\n\n\n",shape1, shape2 , countmatrix.matrix.nnz());
    countmatrix

}


#[cfg(test)]
mod test{
    use std::collections::HashMap;
    use ndarray::{arr2};
    use crate::consistent_genes::{GeneId, CB};

    use super::countmap_to_matrix;

    #[test]
    fn test_countmatrix(){

        let mut countmap: HashMap<(CB, GeneId), usize> = HashMap::new();
        countmap.insert((CB(0), GeneId(0)), 10);
        countmap.insert((CB(0), GeneId(1)), 1);
        countmap.insert((CB(1), GeneId(0)), 0);  // lets see what happens with empty counts
        countmap.insert((CB(1), GeneId(1)), 5);

        let gene_vector = vec!["geneA".to_string(), "geneB".to_string()];

        let cmat = countmap_to_matrix(&countmap, gene_vector);

        let dense_mat = cmat.matrix.to_dense();
        let expected = arr2(&[[ 10, 1],
                                                                  [ 0,  5]]);
        assert_eq!(dense_mat, expected);

        assert_eq!(cmat.cbs, vec!["AAAAAAAAAAAAAAAA".to_string(),"AAAAAAAAAAAAAAAC".to_string()]);
    }

    #[test]
    fn test_countmatrix_equal(){

        // two cells two genes
        // first cell has both genes
        // second cell has only the second gene
        let mut countmap1: HashMap<(CB, GeneId), usize> = HashMap::new();
        countmap1.insert((CB(0), GeneId(0)), 10);
        countmap1.insert((CB(0), GeneId(1)), 1);
        countmap1.insert((CB(1), GeneId(0)), 0);  // lets see what happens with empty counts
        countmap1.insert((CB(1), GeneId(1)), 5);

        let gene_vector = vec!["geneA".to_string(), "geneB".to_string()];

        let cmat1 = countmap_to_matrix(&countmap1, gene_vector);

        // a version with permuated genes
        let mut countmap2: HashMap<(CB, GeneId), usize> = HashMap::new();
        countmap2.insert((CB(0), GeneId(1)), 10);
        countmap2.insert((CB(0), GeneId(0)), 1);
        countmap2.insert((CB(1), GeneId(1)), 0);  // lets see what happens with empty counts
        countmap2.insert((CB(1), GeneId(0)), 5);

        let gene_vector = vec!["geneB".to_string(), "geneA".to_string()];
        let cmat2 = countmap_to_matrix(&countmap2, gene_vector);

        println!("{:?}", cmat1.to_map());
        println!("{:?}", cmat2.to_map());

        assert!(cmat1.is_equal(&cmat2))

    }

    fn test_against_bustools_count(){

    }
}