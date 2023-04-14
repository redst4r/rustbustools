use std::{collections::HashMap, io::{BufReader, BufRead, Write}, fs::File};

use sprs::{io::{read_matrix_market, write_matrix_market}, TriMat};

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

    pub fn to_map(&self) -> HashMap<(String, String), usize> {
        let mut h1: HashMap<(String, String), usize> = HashMap::new();

        for (value, (i,j)) in self.matrix.iter(){
            h1.insert((self.cbs[i].clone(), self.genes[j].clone()), *value);
        }
        h1
    }

    pub fn is_equal(&self, other: &Self) -> bool{
        let h1= self.to_map();
        let h2= other.to_map();
        h1 == h2
    }

    pub fn from_disk(mtx_file: &str, cbfile: &str, genefile: &str) -> CountMatrix{
        let mat: TriMat<usize> = read_matrix_market(mtx_file).unwrap_or_else(|_| panic!("{} not found", mtx_file));
        let matrix: sprs::CsMat<usize> = mat.to_csr();

        let fh = File::open(cbfile).unwrap_or_else(|_| panic!("{} not found", cbfile)); 
        // Read the file line by line, and return an iterator of the lines of the file.
        let cbs: Vec<String> = BufReader::new(fh).lines().collect::<Result<_, _>>().unwrap(); 

        let fh = File::open(genefile).unwrap_or_else(|_| panic!("{} not found", genefile)); 
        let genes: Vec<String> = BufReader::new(fh).lines().collect::<Result<_, _>>().unwrap(); 

        CountMatrix{
            matrix,
            cbs,
            genes,
        }
    }

    pub fn write(&self, foldername: &str){

        let mfile = format!("{}/gene.mtx", foldername);
        let cbfile = format!("{}/gene.barcodes.txt", foldername);
        let genefile = format!("{}/gene.genes.txt", foldername);

        write_matrix_market(mfile, &self.matrix).unwrap();

        let mut fh_cb = File::create(cbfile).unwrap();
        let mut fh_gene = File::create(genefile).unwrap();

        for cb in self.cbs.iter(){
            fh_cb.write_all(format!("{}\n", cb).as_bytes()).unwrap();
        }

        for g in self.genes.iter(){
            fh_gene.write_all(format!("{}\n", g).as_bytes()).unwrap();
        }
    }
}