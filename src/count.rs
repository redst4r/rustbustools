use crate::consistent_genes::{find_consistent, Ec2GeneMapper, Genename, MappingResult, CB};
use crate::countmatrix::CountMatrix;
use crate::io::{group_record_by_cb_umi, BusFolder, BusReader, BusRecord};
use crate::iterators::CellGroupIterator;
use crate::utils::{get_progressbar, int_to_seq};
use sprs;
use std::collections::HashMap;
use std::time::Instant;

///
/// ## emulating bustools count
/// This turns a busfolder into a count matrix.
///
/// The strategy:
/// 1. iterate over CBs, turn each cell (Busrecords from same cell) into an ExpressionVector (HashMap<Genename, u32>).
///    Those all have slightly different key sets
/// 2. Determine ALL genes: from the EC2Gene file
/// 3. turn into a big sparse matrix via expression_vectors_to_matrix()

type ExpressionVector = HashMap<Genename, u32>;

pub fn count_bayesian(bfolder: BusFolder) {
    let bfile = bfolder.get_busfile();
    println!("{}", bfile);

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

    let busiter = BusReader::new(&bfile);
    println!("Getting counts");
    let counts: Vec<f64> = busiter.map(|record| record.COUNT as f64).collect();
    let total_counts: f64 = counts.iter().sum();

    println!(
        "counts total {}, Entries {}",
        total_counts as u64,
        counts.len()
    );
    println!("summing");

    // let total: f64 = counts.iter().sum();
    // let mut r = rand::thread_rng();
    // println!("RV");
    // let n = Multinomial::new(&counts, total as u64).unwrap().sample(&mut r);
    // println!("{:?}", n[0]);
}

pub fn count(bfolder: &BusFolder, ignore_multimapped: bool) -> CountMatrix {
    /*
    busfile to count matrix, analogous to "bustools count"
    */

    let cb_iter = bfolder.get_iterator().groupby_cb();

    println!("determine size of iterator");
    let now = Instant::now();
    let total_records = bfolder.get_cb_size();
    let elapsed_time: std::time::Duration = now.elapsed();
    println!(
        "determined size of iterator {} in {:?}",
        total_records, elapsed_time
    );

    let mut all_expression_vector: HashMap<CB, ExpressionVector> = HashMap::new();
    let now = Instant::now();

    let bar = get_progressbar(total_records as u64);

    for (counter, (cb, record_list)) in cb_iter.enumerate() {
        //}.take(1_000_000){
        let s = records_to_expression_vector(record_list, &bfolder.ec2gene, ignore_multimapped);

        // this will also insert emtpy cells (i.e. their records are all multimapped)
        all_expression_vector.insert(CB(cb), s);

        if counter % 10_000 == 0 {
            bar.inc(10_000)
        }
    }

    let elapsed_time = now.elapsed();
    println!("done in {:?}", elapsed_time);

    //collect all genes
    let genelist_vector: Vec<Genename> = bfolder.ec2gene.get_gene_list();
    println!(" genes {}", genelist_vector.len());

    // todo: whats the point of this conversion from Vec<Genename> -> Vec<&Genename>
    let mut genelist_vector2 = genelist_vector.iter().collect::<Vec<&Genename>>();

    genelist_vector2.sort();

    assert!(genelist_vector2.contains(&&Genename("ENSG00000000003.14".to_string())));

    let countmatrix = expression_vectors_to_matrix(all_expression_vector, genelist_vector2);
    println!(
        "{:?} nnz {}",
        countmatrix.matrix.shape(),
        countmatrix.matrix.nnz()
    );
    countmatrix
}

fn records_to_expression_vector(
    record_list: Vec<BusRecord>,
    // ec2gene: &HashMap<u32, HashSet<String>>,
    eg_mapper: &Ec2GeneMapper,
    ignore_multimapped: bool,
) -> ExpressionVector {
    /*
    turn the list of records of a single CB into a expression vector: per gene, how many umis are observed
    TODO this doesnt consider multiple records with same umi/cb, but EC mapping to different genes
    */
    let mut expression_vector: ExpressionVector = HashMap::new(); // gene -> count
    let mut _multimapped = 0_u32;
    let mut _inconsistant = 0_u32;

    // first, group the records by UMI
    // TODO: EXPENSIVE!! 25k/s
    let cb_umi_grouped = group_record_by_cb_umi(record_list);

    for ((_cb, _umi), records) in cb_umi_grouped {
        // all records coresponding to the same UMI

        let m: MappingResult = if ignore_multimapped {
            // means: If the records map to more than one gene, just treat as unmappable
            match records.len() {
                1 => find_consistent(&records, eg_mapper), // single record, still has to resolve to a single gene!
                0 => panic!(),
                _ => MappingResult::Inconsistent, // if theres more than one record just skip (we dont even try to resolve)
            }
        } else {
            find_consistent(&records, eg_mapper)
        };

        match m {
            // mapped to a single gene: update count!
            MappingResult::SingleGene(g) => {
                let gname = eg_mapper.resolve_gene_id(g);
                let val = expression_vector.entry(gname).or_insert(0);
                *val += 1;
            }
            MappingResult::Multimapped(_) => _multimapped += 1,
            MappingResult::Inconsistent => _inconsistant += 1,
        }

        // let consistent_genes: HashSet<GeneId>;
        // // let consistent_genes: HashSet<String>;
        // if !ignore_multimapped {
        //     consistent_genes = find_consistent(&records, eg_mapper);
        //     // consistent_genes= find_consistent(&records, ec2gene);
        // } else if records.len() > 1 {
        //     continue;
        // } else {
        //     // consistent_genes = ec2gene.get(&records[0].EC).unwrap().clone();
        //     consistent_genes = eg_mapper.get_genes(EC(records[0].EC)).clone();
        // }

        // if consistent_genes.len() > 1 {
        //     //multimapped
        //     _multimapped += 1;
        // } else if consistent_genes.len() == 1 {
        //     //single gene
        //     let g = consistent_genes.iter().next().unwrap(); // Set to first element
        //     let gname = eg_mapper.resolve_gene_id(*g);
        //     let val = expression_vector.entry(gname).or_insert(0);
        //     *val += 1;
        // } else {
        //     // inconsistant
        //     _inconsistant += 1
        // }
    }
    expression_vector
}

fn expression_vectors_to_matrix(
    all_expression_vector: HashMap<CB, ExpressionVector>,
    genelist: Vec<&Genename>,
) -> CountMatrix {
    /*
    turn an collection of expression vector (from many cells)
    into a sparse count matrix
    */

    // sparse matrix indices
    let mut ii: Vec<usize> = Vec::new();
    let mut jj: Vec<usize> = Vec::new();
    let mut vv: Vec<i32> = Vec::new();

    // the cell barcodes, same order as in the matrix
    let mut cbs: Vec<CB> = Vec::new();

    // we need matrix to be sorted
    // for each CB, map it to an index in sorted order
    // let cbs: Vec<u64> = all_expression_vector.keys().;

    // mapping for gene-> index/order
    // this tells us which column in the matrix we need to insert
    let mut gene2index: HashMap<&Genename, usize> = HashMap::new();
    for (i, g) in genelist.iter().enumerate() {
        gene2index.insert(g, i);
    }

    for (i, (cb, expr_vec)) in all_expression_vector.iter().enumerate() {
        for (gene, count) in expr_vec {
            ii.push(i);

            let gindex = gene2index
                .get(gene)
                .unwrap_or_else(|| panic!("{:?} not found", gene));
            jj.push(*gindex);
            vv.push(*count as i32)
        }
        cbs.push(*cb)
    }

    let c: sprs::TriMat<i32> =
        sprs::TriMat::from_triplets((cbs.len(), genelist.len()), ii, jj, vv);
    let b: sprs::CsMat<_> = c.to_csr();

    let cbs_seq: Vec<String> = cbs.into_iter().map(|x| int_to_seq(x.0, 16)).collect();
    // let gene_seq: Vec<String> = genelist.into_iter().map(|x|x.clone()).collect();
    let gene_seq: Vec<String> = genelist.into_iter().map(|x| x.0.to_string()).collect();
    CountMatrix::new(b, cbs_seq, gene_seq)
}

#[cfg(test)]
mod test {
    use super::count;
    use crate::{
        consistent_genes::{Ec2GeneMapper, Genename, EC},
        count::records_to_expression_vector,
        io::{BusFolder, BusRecord},
        utils::vec2set,
    };
    use std::collections::{HashMap, HashSet};

    // #[test]
    fn test_itt() {
        use sprs::io::write_matrix_market;
        // let t2g_file = String::from("/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt");
        // let foldername = String::from("/home/michi/mounts/TB4drive/ISB_data/MNGZ01/MS_processed/S1/kallisto/sort_bus/bus_output");
        let t2g_file = "/home/michi/bus_testing/transcripts_to_genes.txt";
        let foldername = "/home/michi/bus_testing/bus_output";

        let b = BusFolder::new(foldername, t2g_file);
        let count_matrix = count(&b, false);

        // write_sprs_to_file(count_matrix.matrix, "/tmp/test.mtx");
        write_matrix_market("/tmp/test.mtx", &count_matrix.matrix).unwrap();

        // count_bayesian(b)
    }

    #[test]
    fn test_records_to_expression_vector() {
        let ec0: HashSet<Genename> =
            vec2set(vec![Genename("G1".to_string()), Genename("G2".to_string())]);
        let ec1: HashSet<Genename> = vec2set(vec![Genename("G1".to_string())]);
        let ec2: HashSet<Genename> = vec2set(vec![Genename("G2".to_string())]);
        let ec3: HashSet<Genename> =
            vec2set(vec![Genename("G1".to_string()), Genename("G2".to_string())]);

        let ec_dict: HashMap<EC, HashSet<Genename>> = HashMap::from([
            (EC(0), ec0.clone()),
            (EC(1), ec1.clone()),
            (EC(2), ec2.clone()),
            (EC(3), ec3.clone()),
        ]);

        let es = Ec2GeneMapper::new(ec_dict);

        // those three records are consistent with G1
        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };

        // // those records are consistent with G1 and G2, hence multimapped
        let r4 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 2, FLAG: 0 };
        let r5 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 2, FLAG: 0 };
        let r6 = BusRecord { CB: 0, UMI: 2, EC: 3, COUNT: 2, FLAG: 0 };

        // // those records are inconsistent with G1 vs G2
        let r7 = BusRecord { CB: 0, UMI: 3, EC: 0, COUNT: 2, FLAG: 0 };
        let r8 = BusRecord { CB: 0, UMI: 3, EC: 1, COUNT: 2, FLAG: 0 };
        let r9 = BusRecord { CB: 0, UMI: 3, EC: 2, COUNT: 2, FLAG: 0 };

        // // those records are consistent with G1
        let r10 = BusRecord { CB: 0, UMI: 5, EC: 1, COUNT: 2, FLAG: 0 };
        let r11 = BusRecord { CB: 0, UMI: 5, EC: 0, COUNT: 2, FLAG: 0 };

        // // those records are consistent with G2
        let r12 = BusRecord { CB: 0, UMI: 4, EC: 2, COUNT: 2, FLAG: 0 };
        let r13 = BusRecord { CB: 0, UMI: 4, EC: 0, COUNT: 2, FLAG: 0 };

        let records0 = vec![r1.clone(), r2.clone()];
        let c0 = records_to_expression_vector(records0, &es, false);
        assert_eq!(c0, HashMap::from([(Genename("G1".to_string()), 1)]));

        let records1 = vec![r1.clone(), r2.clone(), r10.clone(), r11.clone()];
        let c1 = records_to_expression_vector(records1, &es, false);
        assert_eq!(c1, HashMap::from([(Genename("G1".to_string()), 2)]));

        let records2 = vec![r4.clone(), r5.clone(), r6.clone()];
        let c2 = records_to_expression_vector(records2, &es, false);
        assert_eq!(c2, HashMap::from([]));

        let records3 = vec![r1, r2, r4, r5, r6, r7, r8, r9, r10, r11, r12, r13];
        let c3 = records_to_expression_vector(records3, &es, false);
        assert_eq!(
            c3,
            HashMap::from([
                (Genename("G1".to_string()), 2),
                (Genename("G2".to_string()), 1)
            ])
        );
    }
}
