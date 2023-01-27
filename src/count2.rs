use std::collections::{HashMap, HashSet, BTreeSet};
use std::fs::File;
use std::time::Instant;
use sprs::io::write_matrix_market;

use crate::consistent_genes::{find_consistent, Ec2GeneMapper};
use crate::iterators::CbUmiIterator;
use crate::io::BusFolder;
use crate::utils::{get_progressbar, int_to_seq};

use std::io::Write;

pub struct CountMatrix{
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
            fh_cb.write(format!("{}\n", cb).as_bytes()).unwrap();
        }

        for g in self.genes{
            fh_gene.write(format!("{}\n", g).as_bytes()).unwrap();
        }
    }
}

fn countmap_to_matrix(countmap: &HashMap<(u64, String), usize>, gene_vector: Vec<&String>) -> CountMatrix{

    // get all CBs, a BTreeSet gives us order for free
    // let cb_set: BTreeSet<u64> = BTreeSet::new();
    println!("getting all CBs");
    let all_cbs = countmap.keys().map(|(cb, _gene)| cb).collect::<BTreeSet<_>>();
    // println!("getting all genes");
    // let all_genes = countmap.keys().map(|(cb, gene)| *gene.clone()).collect::<BTreeSet<_>>();

    println!("building index");
    let cb_ix = all_cbs.iter().enumerate().map(|(ix, cb)|(cb.clone(), ix)).collect::<HashMap<_,_>>();
    let gene_ix = gene_vector.iter().enumerate().map(|(ix, cb)|(cb.clone(), ix)).collect::<HashMap<_,_>>();

    // sparse matrix indices
    let mut ii: Vec<usize> = Vec::new();
    let mut jj: Vec<usize> = Vec::new();
    let mut vv: Vec<usize> = Vec::new();

    for ((cb, gene), counter) in countmap{

        let cbi = cb_ix.get(&cb).unwrap();
        let genei = gene_ix.get(&gene).expect(&format!("{} missing", gene));
        ii.push(*cbi);
        jj.push(*genei);
        vv.push(*counter);
    }

    let c: sprs::TriMat<usize> = sprs::TriMat::from_triplets(
        (cb_ix.len(), gene_ix.len()),
        ii,
        jj,
        vv
    );

    let b: sprs::CsMat<_> = c.to_csr();

    let cbs_seq: Vec<String> = all_cbs.into_iter().map(|x|int_to_seq(*x, 16)).collect();
    let gene_seq: Vec<String> = gene_vector.into_iter().map(|x|x.clone()).collect();
    let cmat = CountMatrix::new(b, cbs_seq, gene_seq);

    cmat

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
    let EGM = Ec2GeneMapper::new(ec2gene);


    // CB,gene -> count
    let mut all_expression_vector: HashMap<(u64, u32), usize> = HashMap::new();
    let bar = get_progressbar(total_records as u64);

    let mut n_mapped = 0;
    let mut n_multi_inconsistent = 0;

    let mut counter = 0;

    let now = Instant::now();

    for ((cb, _umi), record_list) in cbumi_iter {
        let consistent_genes: HashSet<u32>;

        if !ignore_multimapped{
                consistent_genes= find_consistent(&record_list, &EGM);
        }else{
            if record_list.len() > 1{
                continue
            }
            else{
                // consistent_genes = ec2gene.get(&record_list[0].EC).unwrap().clone();
                consistent_genes = EGM.get_genes(record_list[0].EC).clone();
            }
        }
       
        // uniquely mapped to a single gene
        if consistent_genes.len() == 1{
            let g = consistent_genes.into_iter().next().unwrap().clone();
            let key = (cb, g);
            
            let current_count = all_expression_vector.entry(key).or_insert(0);
            *current_count+=1;

            n_mapped += 1;
        }
        else{
            // multimapped or inconsistent
            n_multi_inconsistent += 1;
        }

        if counter % 1_000_000 == 0{
            bar.inc(1_000_000);
        }
        counter += 1;
    }

    let elapsed_time = now.elapsed();
    println!("Mapped {}, multi-discard {} in {:?}", n_mapped, n_multi_inconsistent, elapsed_time); 

    // turn from geneid to String
    let EEE: HashMap<(u64, String), usize> = all_expression_vector
        .into_iter()
        .map(|((cb, geneid), v)| ((cb, EGM.resolve_gene_id(geneid)), v)).collect();

    //collect all genes
    // let mut genelist = HashSet::new(); 
    // for glist in ec2gene.values(){
    //     genelist.extend(glist)
    // }
    // let mut genelist_vector :Vec<&String>= genelist.into_iter().collect::<Vec<&String>>();
    let ngenes = EGM.int_to_gene.len();
    let genelist_vector: Vec<String> = (0..ngenes).map(|k| EGM.resolve_gene_id(k as u32)).collect();
    let mut genelist_vector2 = genelist_vector.iter().collect::<Vec<&String>>();

    genelist_vector2.sort();

    let countmatrix = countmap_to_matrix(&EEE, genelist_vector2 );
    println!("{:?} nnz {}", countmatrix.matrix.shape(), countmatrix.matrix.nnz());
    countmatrix

}
