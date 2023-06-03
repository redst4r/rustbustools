use crate::consistent_genes::{find_consistent, Ec2GeneMapper, GeneId, Genename, CB, MappingResult};
use crate::countmatrix::CountMatrix;
use crate::io::{BusFolder, BusRecord};
use crate::iterators::CbUmiGroupIterator;
use crate::multinomial::multinomial_sample;
use crate::utils::{get_progressbar, int_to_seq};
use sprs::DenseVector;
use std::collections::{BTreeSet, HashMap};
use std::time::Instant;

///
/// ## emulating bustools count
/// This turns a busfolder into a count matrix.
/// 
/// Slightly different strategy as count.rs:
/// 1. iterate over CB/UMI, turn into (possibly) count for a gene (if not multimapped) via count_from_record_list()
/// 2. this creates  HashMap<(CB, GeneId), usize> directly, to be turned into a sparse CountMatrix
pub fn countmap_to_matrix(
    countmap: &HashMap<(CB, GeneId), usize>,
    gene_vector: Vec<Genename>,
) -> CountMatrix {
    // get all CBs, a BTreeSet gives us order for free
    // let cb_set: BTreeSet<u64> = BTreeSet::new();
    println!("getting all CBs");
    let all_cbs = countmap
        .keys()
        .map(|(cb, _gene)| cb)
        .collect::<BTreeSet<_>>();
    // println!("getting all genes");
    // let all_genes = countmap.keys().map(|(cb, gene)| *gene.clone()).collect::<BTreeSet<_>>();

    println!("building index");
    // some issues with the cb.clone: clippy complains!
    // let cb_ix = all_cbs.iter().enumerate().map(|(ix, cb)|(cb.clone(), ix)).collect::<HashMap<_,_>>();
    let cb_ix = all_cbs
        .iter()
        .enumerate()
        .map(|(ix, cb)| (**cb, ix))
        .collect::<HashMap<_, _>>();

    // sparse matrix indices
    let mut ii: Vec<usize> = Vec::new();
    let mut jj: Vec<usize> = Vec::new();
    let mut vv: Vec<usize> = Vec::new();

    for ((cb, geneid), counter) in countmap {
        let cbi = cb_ix.get(cb).unwrap();
        let genei = geneid.0 as usize;
        ii.push(*cbi);
        jj.push(genei);
        vv.push(*counter);
    }

    let c: sprs::TriMat<usize> =
        sprs::TriMat::from_triplets((cb_ix.len(), gene_vector.len()), ii, jj, vv);

    let b: sprs::CsMat<_> = c.to_csr();

    let cbs_seq: Vec<String> = all_cbs
        .into_iter()
        .map(|x| int_to_seq(x.0, 16))
        .collect();
    // let gene_seq: Vec<String> = gene_vector.into_iter().map(|x|x.clone()).collect();
    let gene_seq: Vec<String> = gene_vector.into_iter().map(|x| x.0).collect(); //not sure if this does anything

    CountMatrix::new(b, cbs_seq, gene_seq)
}

pub fn baysian_count(bfolder: BusFolder, ignore_multimapped: bool, n_samples: usize) {
    let bfile = bfolder.get_busfile();
    println!("{}", bfile);

    println!("determine size of iterator");
    let now = Instant::now();
    let total_records = bfolder.get_cbumi_size();
    let elapsed_time: std::time::Duration = now.elapsed();
    println!("determined size of iterator {} in {:?}", total_records, elapsed_time);

    let elapsed_time = now.elapsed();
    println!(
        "determined size of iterator {} in {:?}.",
        total_records, elapsed_time
    );

    // handles the mapping between EC and gene
    let egm = &bfolder.ec2gene;

    // prep for the multinomial sample
    println!("Preparing the probability vector for mutlinomial");
    let cbumi_iter_tmp = bfolder.get_iterator().groupby_cbumi();

    let count_vec: Vec<_> = cbumi_iter_tmp
        .flat_map(|(_cbumi, rlist)| rlist.into_iter().map(|r| r.COUNT as f64))
        .collect();

    let total_counts: f64 = count_vec.iter().sum();
    let p_vec: Vec<f64> = count_vec.into_iter().map(|c| c / total_counts).collect();
    println!("Done: {} rercods, {} counts", p_vec.len(), total_counts);

    use probability::prelude::*;
    let mut random_source = source::default(42);

    let mut counter = 0;
    for i in 0..n_samples {
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

        for ((cb, _umi), rlist) in cbumi_iter {
            // inject the sampled numbers into the records

            let indices = current_record_counter..current_record_counter+rlist.len();
            let injected_counts: Vec<u32> = indices.map(|idx| *new_count_sample.index(idx) as u32 ).collect(); // wrning f64->u32
            // let mut injected_records: Vec<BusRecord> = Vec::with_capacity(rlist.len());
            let mut injected_records: Vec<BusRecord> = rlist.clone();

            for i in 0..injected_records.len() {
                // let mut r = injected_records.get_mut(i).expect(&format!("injected_records {}", i));
                let mut r = injected_records.get_mut(i).unwrap_or_else(|| panic!("injected_records {}", i));
                let c = injected_counts.get(i).unwrap_or_else(|| panic!("injected_counts {}", i));
                r.COUNT = *c;
            }

            injected_records.retain(|r| r.COUNT > 0);

            // for (r, new_count) in injected_records.iter_mut().zip(injected_counts.into_iter()){
            //     r.COUNT = new_count;
            //     injected_records.push(r);
            // }
            current_record_counter += rlist.len();

            if injected_records.is_empty() {
                continue;
            }

            match count_from_record_list(&injected_records, egm, ignore_multimapped) {
                MappingResult::SingleGene(g) => {
                    let key = (CB(cb), g);
                    let current_count = all_expression_vector.entry(key).or_insert(0);
                    *current_count += 1;
                    n_mapped += 1;
                },
                MappingResult::Multimapped(_) | MappingResult::Inconsistent => n_multi_inconsistent += 1,
            }

            if counter % 1_000_000 == 0 {
                bar.inc(1_000_000);
            }
            counter += 1;
        }

        let elapsed_time = now.elapsed();
        let fraction_mapped =
            n_multi_inconsistent as f64 / (n_mapped as f64 + n_multi_inconsistent as f64);
        println!(
            "Iteration {}: Mapped {}, multi-discard {} ({}%) in {:?}",
            i,
            n_mapped,
            n_multi_inconsistent,
            100.0 * fraction_mapped,
            elapsed_time
        );

        let genelist_vector: Vec<Genename> = egm.get_gene_list();
        // this is how genes are ordered as by EGM
        // i.e. countmap[cb, i] corresponds to the number of count of genelist_vector[i]

        let countmatrix = countmap_to_matrix(&all_expression_vector, genelist_vector);
        println!(
            "{:?} nnz {}",
            countmatrix.matrix.shape(),
            countmatrix.matrix.nnz()
        );
        println!("finished iteration {}", i)
    }
}

fn count_from_record_list(
    records: &Vec<BusRecord>,
    egmapper: &Ec2GeneMapper,
    ignore_multimapped: bool,
) -> MappingResult {
    // given a set of Records from the same CB/UMI, are they consistent with a particular gene
    // which would yield a signel count

    // watch out for the convoluted logic with ignore_multimapped!!
    
    if ignore_multimapped {
        // means: If the records map to more than one gene, just treat as unmappable
        match records.len() {
            1 => find_consistent(records, egmapper),  // single record, still has to resolve to a single gene!
            0 => panic!(),
            _ => MappingResult::Inconsistent  // if theres more than one record just skip (we dont even try to resolve)
        }
    } else {
        find_consistent(records, egmapper) 
    }
}

pub fn count(bfolder: &BusFolder, ignore_multimapped: bool) -> CountMatrix {
    /*
    busfile to count matrix, analogous to "bustools count"
    */
    let bfile = bfolder.get_busfile();
    println!("{}", bfile);

    let cbumi_iter = bfolder.get_iterator().groupby_cbumi();

    println!("determine size of iterator");
    let now = Instant::now();
    let total_records = bfolder.get_cbumi_size();
    let elapsed_time: std::time::Duration = now.elapsed();
    println!("determined size of iterator {} in {:?}", total_records, elapsed_time);

    // CB,gene_id -> count
    let mut all_expression_vector: HashMap<(CB, GeneId), usize> = HashMap::new();
    let bar = get_progressbar(total_records as u64);

    let mut n_mapped = 0;
    let mut n_multi_inconsistent = 0;

    let now = Instant::now();

    for (counter, ((cb, _umi), record_list)) in cbumi_iter.enumerate() {
        // try to map the records of this CB/UMI into a single gene
        // if let Some(g) = count_from_record_list(&record_list, &bfolder.ec2gene, ignore_multimapped)
        match count_from_record_list(&record_list, &bfolder.ec2gene, ignore_multimapped)
        {
            MappingResult::SingleGene(g) => {
                let key = (CB(cb), g);
                let current_count = all_expression_vector.entry(key).or_insert(0);
                *current_count += 1;
                n_mapped += 1;
            },
            MappingResult::Multimapped(_) | MappingResult::Inconsistent => {
                // multimapped, or not consistently mapped
                n_multi_inconsistent += 1;
                // let cbstr = int_to_seq(cb, 16);
                // let umistr = int_to_seq(_umi, 12);
                //println!("not countable {cbstr}/{umistr} {:?}", record_list);
                // let cgeneids= find_consistent(&record_list, &bfolder.ec2gene);
                // let cgenes: Vec<_> = cgeneids.iter().map(|gid| bfolder.ec2gene.resolve_gene_id(*gid)).collect();
                //println!("{cgenes:?}")
            },
        }

        if counter % 1_000_000 == 0 {
            bar.inc(1_000_000);
        }
    }

    let elapsed_time = now.elapsed();
    println!(
        "Mapped {}, multi-discard {} in {:?}",
        n_mapped, n_multi_inconsistent, elapsed_time
    );

    let genelist_vector: Vec<Genename> = bfolder.ec2gene.get_gene_list();

    // this is how genes are ordered as by EGM
    // i.e. countmap[cb, i] corresponds to the number of count of genelist_vector[i]

    let countmatrix = countmap_to_matrix(&all_expression_vector, genelist_vector);

    let (shape1, shape2) = countmatrix.matrix.shape();
    println!(
        "{}x{} nnz {}\n\n\n",
        shape1,
        shape2,
        countmatrix.matrix.nnz()
    );
    countmatrix
}
