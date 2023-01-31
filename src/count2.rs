use std::collections::{HashMap, HashSet, BTreeSet};
use std::fs::File;
use std::time::Instant;
use sprs::DenseVector;
use sprs::io::write_matrix_market;

use crate::consistent_genes::{find_consistent, Ec2GeneMapper};
use crate::iterators::CbUmiIterator;
use crate::io::{BusFolder, BusRecord};
use crate::multinomial::multinomial_sample;
use crate::utils::{get_progressbar, int_to_seq};

use std::io::Write;

pub struct CountMatrix{
    /*
    Represents a count matrix of Cells vs Genes
    Cells are encoded a Strings already!
    */
    pub matrix: sprs::CsMat<usize>,
    pub cbs: Vec<String>,
    pub genes: Vec<String>,
}
impl CountMatrix {
    pub fn new(matrix: sprs::CsMat<usize>, cbs:Vec<String>, genes: Vec<String>) ->CountMatrix{
        CountMatrix{
            matrix,
            cbs,
            genes,
        }
    }

    pub fn write(self, foldername: &str){

        let mfile = format!("{}/gene.mtx", foldername);
        let cbfile = format!("{}/gene.barcodes.txt", foldername);
        let genefile = format!("{}/gene.genes.txt", foldername);


        write_matrix_market(mfile, &self.matrix).unwrap();

        let mut fh_cb = File::create(cbfile).unwrap();
        let mut fh_gene = File::create(genefile).unwrap();

        for cb in self.cbs{
            fh_cb.write_all(format!("{}\n", cb).as_bytes()).unwrap();
        }

        for g in self.genes{
            fh_gene.write_all(format!("{}\n", g).as_bytes()).unwrap();
        }
    }
}

fn countmap_to_matrix(countmap: &HashMap<(u64, u32), usize>, gene_vector: Vec<String>) -> CountMatrix{

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
        let genei = *geneid as usize;
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

    let cbs_seq: Vec<String> = all_cbs.into_iter().map(|x|int_to_seq(*x, 16)).collect();
    // let gene_seq: Vec<String> = gene_vector.into_iter().map(|x|x.clone()).collect();
    let gene_seq: Vec<String> = gene_vector.into_iter().collect(); //not sure if this does anything
    
    CountMatrix::new(b, cbs_seq, gene_seq)

}


pub fn baysian_count(bfolder: BusFolder, ignore_multimapped:bool, n_samples: usize){
    let bfile = format!("{}/{}", bfolder.foldername, bfolder.busfile);
    println!("{}",bfile);

    let ec2gene = bfolder.ec2gene;

    let cbumi_iter_tmp = CbUmiIterator::new(&bfile);
    println!("determine size of iterator");
    let now = Instant::now();
    let total_records = cbumi_iter_tmp.count();

    let cbumi_iter_tmp = CbUmiIterator::new(&bfile);
    let max_length_records = cbumi_iter_tmp.map(|(_cbumi, rlist)| rlist.len()).max().unwrap();

    let elapsed_time = now.elapsed();
    println!("determined size of iterator {} in {:?}. Longest element: {} in a single CB/UMI", total_records, elapsed_time, max_length_records);

    // handles the mapping between EC and gene
    let egm = Ec2GeneMapper::new(ec2gene);


    // prep for the multinomial sample
    println!("Preparing the probability vector for mutlinomial");
    let cbumi_iter_tmp = CbUmiIterator::new(&bfile);
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
        let mut all_expression_vector: HashMap<(u64, u32), usize> = HashMap::new();
        let mut n_mapped = 0;
        let mut n_multi_inconsistent = 0;

        // subsample the count vector
        println!("Iteration {}: Mutlinomial sample", i);
        let new_count_sample = multinomial_sample(total_counts as u64, &p_vec, &mut random_source);
        println!("Done");

        let cbumi_iter = CbUmiIterator::new(&bfile);
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


            if let Some(g) = count_from_record_list(injected_records, &egm, ignore_multimapped){
                // the records could be made into a single count for gene g
                let key = (cb, g);
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

        let ngenes = egm.int_to_gene.len();
        let genelist_vector: Vec<String> = (0..ngenes).map(|k| egm.resolve_gene_id(k as u32)).collect();
        // this is how genes are ordered as by EGM
        // i.e. countmap[cb, i] corresponds to the number of count of genelist_vector[i]
    
        let countmatrix = countmap_to_matrix(&all_expression_vector, genelist_vector );
        println!("{:?} nnz {}", countmatrix.matrix.shape(), countmatrix.matrix.nnz()); 
        println!("finished iteration {}", i)         
    }
}



fn count_from_record_list(records: Vec<BusRecord>, egmapper: &Ec2GeneMapper, ignore_multimapped:bool) -> Option<u32>{

    // given a set of Records from the same CB/UMI, are they consistent with a particular gene
    // which would yield a signel count
    let consistent_genes: HashSet<u32>;

    if !ignore_multimapped{
            consistent_genes= find_consistent(&records, egmapper);
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

pub fn count(bfolder: BusFolder, ignore_multimapped:bool) -> CountMatrix {
    /*
    busfile to count matrix, analogous to "bustools count"
    */
    let bfile = format!("{}/{}", bfolder.foldername, bfolder.busfile);
    println!("{}",bfile);

    let ec2gene = bfolder.ec2gene;
    let cbumi_iter = CbUmiIterator::new(&bfile);

    let cbumi_iter_tmp = CbUmiIterator::new(&bfile);
    println!("determine size of iterator");
    let now = Instant::now();
    let total_records = cbumi_iter_tmp.count();

    let cbumi_iter_tmp = CbUmiIterator::new(&bfile);
    let max_length_records = cbumi_iter_tmp.map(|(_cbumi, rlist)| rlist.len()).max().unwrap();

    let elapsed_time = now.elapsed();

    println!("determined size of iterator {} in {:?}. Longest element: {} in a single CB/UMI", total_records, elapsed_time, max_length_records);

    // handles the mapping between EC and gene
    let eg_mapper = Ec2GeneMapper::new(ec2gene);


    // CB,gene_id -> count
    let mut all_expression_vector: HashMap<(u64, u32), usize> = HashMap::new();
    let bar = get_progressbar(total_records as u64);

    let mut n_mapped = 0;
    let mut n_multi_inconsistent = 0;

    let now = Instant::now();

    for (counter, ((cb, _umi), record_list)) in cbumi_iter.enumerate() {

        // try to map the records of this CB/UMI into a single gene
        if let Some(g) = count_from_record_list(record_list, &eg_mapper, ignore_multimapped){
            // the records could be made into a single count for gene g
            let key = (cb, g);
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
    }

    let elapsed_time = now.elapsed();
    println!("Mapped {}, multi-discard {} in {:?}", n_mapped, n_multi_inconsistent, elapsed_time); 

    let ngenes = eg_mapper.int_to_gene.len();
    let genelist_vector: Vec<String> = (0..ngenes).map(|k| eg_mapper.resolve_gene_id(k as u32)).collect();
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
    use super::countmap_to_matrix;

    #[test]
    fn test_countmatrix(){

        let mut countmap: HashMap<(u64, u32), usize> = HashMap::new();
        countmap.insert((0, 0), 10);
        countmap.insert((0, 1), 1);
        countmap.insert((1, 0), 0);  // lets see what happens with empty counts
        countmap.insert((1, 1), 5);

        let gene_vector = vec!["geneA".to_string(), "geneB".to_string()];

        let cmat = countmap_to_matrix(&countmap, gene_vector);

        let dense_mat = cmat.matrix.to_dense();
        let expected = arr2(&[[ 10, 1],
                                                                  [ 0,  5]]);
        assert_eq!(dense_mat, expected);

        assert_eq!(cmat.cbs, vec!["AAAAAAAAAAAAAAAA".to_string(),"AAAAAAAAAAAAAAAC".to_string()]);
    }
}