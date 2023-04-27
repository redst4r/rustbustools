use std::collections::{HashMap, HashSet};
use std::time::Instant;
use crate::countmatrix::CountMatrix;
use crate::io::{BusFolder, BusRecord, group_record_by_cb_umi, BusReader};
use crate::iterators::CellGroupIterator;
use crate::utils::{get_progressbar, int_to_seq};
use sprs;
use crate::consistent_genes::{find_consistent, Ec2GeneMapper, GeneId, Genename, EC, CB};


pub fn count_bayesian(bfolder: BusFolder) {

    let bfile = bfolder.get_busfile();
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

    let busiter = BusReader::new(&bfile);
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



pub fn count(bfolder: &BusFolder, ignore_multimapped: bool) -> CountMatrix {
    /*
    busfile to count matrix, analogous to "bustools count"
    */

    let cb_iter = bfolder.get_iterator().groupby_cb();
    let cb_iter_tmp = bfolder.get_iterator().groupby_cb();

    println!("determine size of iterator");
    let now = Instant::now();
    let total_records = cb_iter_tmp.count();
    let elapsed_time = now.elapsed();
    println!("determined size of iterator {} in {:?}", total_records, elapsed_time);

    let mut all_expression_vector: HashMap<CB, HashMap<Genename, u32>> = HashMap::new();
    let now = Instant::now();

    // let bar = ProgressBar::new_spinner();
    let bar = get_progressbar(total_records as u64);

    for (counter, (cb, record_list)) in cb_iter.enumerate() {//}.take(1_000_000){
        let s = records_to_expression_vector(record_list, &bfolder.ec2gene, ignore_multimapped);

        // this will also insert emtpy cells (i.e. their records are all multimapped)
        all_expression_vector.insert(CB(cb), s);

        if counter % 10_000 == 0{
            bar.inc(10_000)
        }
    }

    let elapsed_time = now.elapsed();
    println!("done in {:?}", elapsed_time); 
    
    //collect all genes
    let genelist_vector: Vec<Genename> = bfolder.ec2gene.get_gene_list();
    println!(" genes {}", genelist_vector.len());

    let mut genelist_vector2 = genelist_vector.iter().collect::<Vec<&Genename>>();

    genelist_vector2.sort();

    assert!(genelist_vector2.contains(&&Genename("ENSG00000000003.14".to_string())));

    
    let countmatrix = expression_vectors_to_matrix(all_expression_vector, genelist_vector2 );
    println!("{:?} nnz {}", countmatrix.matrix.shape(), countmatrix.matrix.nnz());
    countmatrix

}



// fn records_to_expression_vector_groupby(record_list: Vec<BusRecord>, ec2gene: &HashMap<u32, HashSet<String>>) -> HashMap<String, u32>{
//     let mut expression_vector: HashMap<String, u32> = HashMap::new(); // gene -> count
//     let mut _multimapped = 0;
//     let mut _inconsistant= 0;

//     // first, group the records by UMI
//     // TODO: EXPENSIVE!! 25k/s

//     use itertools::Itertools;
//     let cb_umi_grouped = record_list.into_iter().group_by(|a|(a.CB, a.UMI));

//     for (_key, group) in cb_umi_grouped.into_iter(){
//         let records: Vec<BusRecord> = group.collect();
//         // all records coresponding to the same UMI
//         // TODO: EXPENSIVE!! 25k/s
//         let consistent_genes= find_consistent(&records, ec2gene);

//         // let consistent_genes = ec2gene.get(&records[0].EC).unwrap();
//         // let consistent_genes = vec![records[0].EC];


//         if consistent_genes.len() > 1{
//             //multimapped
//             _multimapped += 1;
//         }
//         else if consistent_genes.len() == 1 {
//             //single gene
//             // let g = consistent_genes.drain().next().unwrap();

//             // let v: Vec<String> = consistent_genes.into_iter().collect();  // Set to Vec
//             // let g = &v[0];

//             let g = consistent_genes.iter().next().unwrap();  // Set to first element


//             let val = expression_vector.entry(g.to_string()).or_insert(0);
//             *val += 1;        
//         }
//         else{
//             // inconsistant
//             _inconsistant += 1
//         }
//     }
//     // println!("{}, {}",_multimapped, _inconsistant);

//     expression_vector

// }


fn records_to_expression_vector(
    record_list: Vec<BusRecord>, 
    // ec2gene: &HashMap<u32, HashSet<String>>, 
    eg_mapper: &Ec2GeneMapper, 
    ignore_multimapped: bool
) -> HashMap<Genename, u32>{
    /*
    turn the list of records of a single CB into a expression vector: per gene, how many umis are observed
    TODO this doesnt consider multiple records with same umi/cb, but EC mapping to different genes
    */
    let mut expression_vector: HashMap<Genename, u32> = HashMap::new(); // gene -> count
    let mut _multimapped = 0_u32;
    let mut _inconsistant= 0_u32;

    // first, group the records by UMI
    // TODO: EXPENSIVE!! 25k/s
    let cb_umi_grouped = group_record_by_cb_umi(record_list);

    for ((_cb, _umi), records) in cb_umi_grouped{
        // all records coresponding to the same UMI

        let consistent_genes: HashSet<GeneId>;
        // let consistent_genes: HashSet<String>;
        if ! ignore_multimapped{
            consistent_genes= find_consistent(&records, eg_mapper);
            // consistent_genes= find_consistent(&records, ec2gene);
        }else if records.len() > 1{
                continue;
        }
        else{
            // consistent_genes = ec2gene.get(&records[0].EC).unwrap().clone();
            consistent_genes = eg_mapper.get_genes(EC(records[0].EC)).clone();
        }

        if consistent_genes.len() > 1{
            //multimapped
            _multimapped += 1;
        }
        else if consistent_genes.len() == 1 {
            //single gene
            let g = consistent_genes.iter().next().unwrap();  // Set to first element
            let gname = eg_mapper.resolve_gene_id(*g);
            let val = expression_vector.entry(gname).or_insert(0);
            *val += 1;        
        }
        else{
            // inconsistant
            _inconsistant += 1
        }
    }
    expression_vector
}

fn expression_vectors_to_matrix(
    all_expression_vector: HashMap<CB, HashMap<Genename, u32>>, 
    genelist: Vec<&Genename>
) -> CountMatrix {
    /*
    turn an collection of expression vector (from many cells)
    into a sparse count matrix
    */

    // sparse matrix indices
    let mut ii: Vec<usize> = Vec::new();
    let mut jj: Vec<usize> = Vec::new();
    let mut vv: Vec<usize> = Vec::new();

    // the cell barcodes, same order as in the matrix
    let mut cbs: Vec<CB> = Vec::new(); 

    // we need matrix to be sorted
    // for each CB, map it to an index in sorted order
    // let cbs: Vec<u64> = all_expression_vector.keys().; 


    // mapping for gene-> index/order
    let mut gene2index: HashMap<&Genename, usize> = HashMap::new();
    for (i, g) in genelist.iter().enumerate(){
        gene2index.insert(g, i);
    }

    for (i, (cb, expr_vec)) in all_expression_vector.iter().enumerate(){
        for (gene, count) in expr_vec{
            ii.push(i);

            let gindex = gene2index.get(gene).unwrap_or_else(|| panic!("{:?} not found", gene));
            jj.push(*gindex);
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


    let cbs_seq: Vec<String> = cbs.into_iter().map(|x|int_to_seq(x.0, 16)).collect();
    // let gene_seq: Vec<String> = genelist.into_iter().map(|x|x.clone()).collect();
    let gene_seq: Vec<String> = genelist.into_iter().map(|x| x.0.to_string()).collect();
    CountMatrix::new(b, cbs_seq, gene_seq)
}



#[cfg(test)]
mod test{
    use std::{collections::{HashMap, HashSet}, fs, time::Instant};
    use crate::{io::{BusRecord, BusFolder}, count::records_to_expression_vector, countmatrix::{CountMatrix}, count2, consistent_genes::{Ec2GeneMapper, Genename, EC}, utils::vec2set};
    use super::count;
    
    // #[test]    
    fn test_itt(){
        use sprs::io::write_matrix_market;
        // let t2g_file = String::from("/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt");
        // let foldername = String::from("/home/michi/mounts/TB4drive/ISB_data/MNGZ01/MS_processed/S1/kallisto/sort_bus/bus_output");
        let t2g_file = "/home/michi/bus_testing/transcripts_to_genes.txt";
        let foldername ="/home/michi/bus_testing/bus_output";
    
        let b = BusFolder::new(foldername, t2g_file);
        let count_matrix = count(&b, false);
    
        // write_sprs_to_file(count_matrix.matrix, "/tmp/test.mtx");
        write_matrix_market("/tmp/test.mtx", &count_matrix.matrix).unwrap();

        // count_bayesian(b)
    }
    
    #[test]
    fn test_records_to_expression_vector(){

        let ec0: HashSet<Genename> = vec2set(vec![Genename("G1".to_string()), Genename("G2".to_string())]);
        let ec1: HashSet<Genename> = vec2set(vec![Genename("G1".to_string())]);
        let ec2: HashSet<Genename> = vec2set(vec![Genename("G2".to_string())]);
        let ec3: HashSet<Genename> = vec2set(vec![Genename("G1".to_string()), Genename("G2".to_string())]);

        let ec_dict: HashMap<EC, HashSet<Genename>> = HashMap::from([
            (EC(0), ec0.clone()), 
            (EC(1), ec1.clone()), 
            (EC(2), ec2.clone()), 
            (EC(3), ec3.clone()), 
            ]);
    
        let es = Ec2GeneMapper::new(ec_dict);

        // those three records are consistent with G1
        let r1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

        // // those records are consistent with G1 and G2, hence multimapped
        let r4 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 2, FLAG: 0};
        let r5 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 2, FLAG: 0};
        let r6 = BusRecord{CB: 0, UMI: 2, EC: 3, COUNT: 2, FLAG: 0};

        // // those records are inconsistent with G1 vs G2
        let r7= BusRecord{CB: 0, UMI: 3, EC: 0, COUNT: 2, FLAG: 0};
        let r8= BusRecord{CB: 0, UMI: 3, EC: 1, COUNT: 2, FLAG: 0};
        let r9= BusRecord{CB: 0, UMI: 3, EC: 2, COUNT: 2, FLAG: 0};

        // // those records are consistent with G1
        let r10= BusRecord{CB: 0, UMI: 5, EC: 1, COUNT: 2, FLAG: 0};
        let r11= BusRecord{CB: 0, UMI: 5, EC: 0, COUNT: 2, FLAG: 0};

        // // those records are consistent with G2
        let r12 = BusRecord{CB: 0, UMI: 4, EC: 2, COUNT: 2, FLAG: 0};
        let r13= BusRecord{CB: 0, UMI: 4, EC: 0, COUNT: 2, FLAG: 0};

        let records0 = vec![r1.clone(),r2.clone()];
        let c0 = records_to_expression_vector(records0, &es, false);
        assert_eq!(c0, HashMap::from([(Genename("G1".to_string()), 1)]));

        let records1 = vec![r1.clone(),r2.clone(), r10.clone(), r11.clone()];
        let c1= records_to_expression_vector(records1, &es, false);
        assert_eq!(c1, HashMap::from([(Genename("G1".to_string()), 2)]));


        let records2 = vec![r4.clone(),r5.clone(),r6.clone()];
        let c2= records_to_expression_vector(records2, &es, false);
        assert_eq!(c2, HashMap::from([]));


        let records3 = vec![r1,r2, r4,r5,r6,r7,r8,r9,r10, r11,r12,r13];
        let c3 = records_to_expression_vector(records3, &es, false);
        assert_eq!(c3, HashMap::from([(Genename("G1".to_string()), 2), (Genename("G2".to_string()), 1)]));



    }

    #[test]
    fn test_vs_bustools(){
        // comparing our implementation vs kallisto-bustools on a real bus file
        // run: 
        // RUST_BACKTRACE=1 cargo test --release --package rustbustools --lib -- --nocapture count::test::test_vs_bustools
        use std::process::Command;

        let outfolder = "/tmp/bustest_rust";
        let outfolder_kallisto = "/tmp/bustest_kallisto";

        let busfolder = "/home/michi/bus_testing/bus_output_shorter";
        // let t2g = "/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt";
        let t2g = "/home/michi/bus_testing/transcripts_to_genes.txt";
        
        let bfolder = BusFolder::new(busfolder, t2g);
        let tfile = bfolder.get_transcript_file();
        let ecfile = bfolder.get_ecmatrix_file();
        let bfile = bfolder.get_busfile();

        fs::create_dir_all(outfolder).expect("Failed to create outfolder");
        fs::create_dir_all(outfolder_kallisto).expect("Failed to create outfolder_kallisto");
        
        // -------------------
        // Doing our own count
        // -------------------
        let now = Instant::now();
        let c = count(&bfolder, false);
        let elapsed_time = now.elapsed();
        println!("count::count in in {:?}", elapsed_time);
        c.write(outfolder);

        let now = Instant::now();
        let c2 = count2::count(&bfolder, false);
        let elapsed_time = now.elapsed();
        println!("count2::count in in {:?}", elapsed_time);
        assert!(c2.is_equal(&c));

        // -------------------
        // Bustools count
        // -------------------
        // Command
        // println!("source ~/miniconda3_newtry/bin/activate June2021");
        println!("bustools count -o {outfolder_kallisto}/gene -e {ecfile} -g {t2g} -t {tfile} --genecounts  {bfile}");

        let cmat_kallisto = CountMatrix::from_disk(
            &format!("{outfolder_kallisto}/gene.mtx"),
            &format!("{outfolder_kallisto}/gene.barcodes.txt"),
            &format!("{outfolder_kallisto}/gene.genes.txt"),
        );

        let cmat_rust = CountMatrix::from_disk(
            &format!("{outfolder}/gene.mtx"),
            &format!("{outfolder}/gene.barcodes.txt"),
            &format!("{outfolder}/gene.genes.txt"),
        );
        // println!("kallisto \n{:?}", cmat_kallisto.to_map());
        // println!("rust\n{:?}", cmat_rust.to_map());

        let sum1: usize = cmat_kallisto.matrix.iter().map(|(v,_s)| *v).sum();
        let sum2: usize = cmat_rust.matrix.iter().map(|(v,_s)| *v).sum();
        println!("{sum1} {sum2}");

        let h1 =cmat_kallisto.to_map();
        let h2 =cmat_rust.to_map();

        for (k,v) in h1.iter(){
            match h2.get(k){
                Some(v2) => {
                    if *v != *v2{
                        println!("h2 wrong value {k:?}, {v}  <-> {k:?}, {v2}")
                    }
                },
                None => {
                    println!("h1 unique: {k:?}, {v}")
                }
            }
        }
        for (k,v) in h2.iter(){
            match h1.get(k){
                Some(v2) => {
                    if *v != *v2{
                        println!("h1 wrong value {k:?}, {v}  <-> {k:?}, {v2}")
                    }
                },
                None => {
                    println!("h2 unique: {k:?}, {v}")
                }
            }
        }        

        assert!(cmat_kallisto.is_equal(&cmat_rust))

    }
}