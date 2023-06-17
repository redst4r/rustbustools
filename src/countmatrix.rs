//! A sparse countmatrix: cells vs genes
//!
//! Usually the result of a bustools count call or [crate::count::count()]
//!
//! ## Example
//! ```rust, no_run
//! // load from disk
//! # use rustbustools::countmatrix::CountMatrix;
//! let path = "/path/to/bus/count/folder";
//! let cmat = CountMatrix::from_disk(
//!     &format!("{}/gene.mtx", path),
//!     &format!("{}/gene.barcodes.txt", path),
//!     &format!("{}/gene.genes.txt", path),
//!  );
//!  // shorter, assuming standard file names
//!  let cmat = CountMatrix::from_folder(path);
//!
//! // write to disk again
//! // note that the folder must exist already
//! let outpath = "/tmp/rustbustools_test_read_write";
//! if !std::path::Path::new(&outpath).exists() {
//!     std::fs::create_dir(outpath).unwrap();
//! }
//! cmat.write(outpath);
//! ```
use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Write},
};

use sprs::{
    io::{read_matrix_market, write_matrix_market},
    TriMat,
};

/// Countmatrix, cells-by-genes
///
/// Cells and genes are indexed via their string reprensentation
///
pub struct CountMatrix {
    /// sparse count matrix
    pub matrix: sprs::CsMat<i32>,
    cbs: Vec<String>,
    genes: Vec<String>,
}
impl CountMatrix {
    /// create a CountMatrix from a sparse matrix type ([sprs::CsMat]) and name the rows (cells) and columns (genes)
    pub fn new(matrix: sprs::CsMat<i32>, cbs: Vec<String>, genes: Vec<String>) -> CountMatrix {
        CountMatrix { matrix, cbs, genes }
    }

    /// turns the count-matrix into a HashMap for easier comparision to other countmatrices
    fn to_map(&self) -> HashMap<(String, String), i32> {
        // transforms the sparse count matrix into a Hashmap (CB,Gene)-> count
        let mut h1: HashMap<(String, String), i32> = HashMap::new();

        for (value, (i, j)) in self.matrix.iter() {
            h1.insert((self.cbs[i].clone(), self.genes[j].clone()), *value);
        }
        h1
    }

    /// comparing countmatrices. True if they represnet the same cb/gene counts irrespective of ordering
    pub fn is_equal(&self, other: &Self) -> bool {
        // the tricky to be invariant to ordering: transform to hashmap
        let h1 = self.to_map();
        let h2 = other.to_map();
        h1 == h2
    }

    /// get the matrix's shape (nrows, ncols)
    pub fn get_shape(&self) -> (usize, usize) {
        self.matrix.shape()
    }

    /// String repr
    pub fn to_string(&self) -> String {
        format!("Shape: {:?};  nnz {}", self.get_shape(), self.matrix.nnz())
    }

    /// load a countmatrix from disk (kallisto format: mtx + barcodes.txt + genes)
    pub fn from_disk(mtx_file: &str, cbfile: &str, genefile: &str) -> Self {
        // load countmatrix from disk, from matrix-market format
        let mat: TriMat<i32> =
            read_matrix_market(mtx_file).unwrap_or_else(|_| panic!("{} not found", mtx_file));
        let matrix: sprs::CsMat<i32> = mat.to_csr();

        let fh = File::open(cbfile).unwrap_or_else(|_| panic!("{} not found", cbfile));
        // Read the file line by line, and return an iterator of the lines of the file.
        let cbs: Vec<String> = BufReader::new(fh)
            .lines()
            .collect::<Result<_, _>>()
            .unwrap();

        let fh = File::open(genefile).unwrap_or_else(|_| panic!("{} not found", genefile));
        let genes: Vec<String> = BufReader::new(fh)
            .lines()
            .collect::<Result<_, _>>()
            .unwrap();

        CountMatrix { matrix, cbs, genes }
    }

    /// load the countmatrix from a folder, assuming standatd file naming
    pub fn from_folder(foldername: &str) -> Self {
        let mfile = &format!("{}/gene.mtx", foldername);
        let cbfile = &format!("{}/gene.barcodes.txt", foldername);
        let genefile = &format!("{}/gene.genes.txt", foldername);
        CountMatrix::from_disk(mfile, cbfile, genefile)
    }

    /// write the matrix to disk in MatrixMarket format + cell and gene metadata (just like kallisto)
    ///
    /// creates 3 files:
    /// * gene.mtx: the sparse matrix
    /// * gene.barcodes.txt: String representation fo the cell barcodes
    /// * gene.genes.txt: Gene names
    pub fn write(&self, foldername: &str) {
        let mfile = format!("{}/gene.mtx", foldername);
        let cbfile = format!("{}/gene.barcodes.txt", foldername);
        let genefile = format!("{}/gene.genes.txt", foldername);

        write_matrix_market(mfile, &self.matrix).unwrap();

        let mut fh_cb = File::create(cbfile).unwrap();
        let mut fh_gene = File::create(genefile).unwrap();

        for cb in self.cbs.iter() {
            fh_cb.write_all(format!("{}\n", cb).as_bytes()).unwrap();
        }

        for g in self.genes.iter() {
            fh_gene.write_all(format!("{}\n", g).as_bytes()).unwrap();
        }
    }
}

#[cfg(test)]
mod test {
    use crate::{
        consistent_genes::{GeneId, Genename, CB},
        count2::countmap_to_matrix,
    };
    use ndarray::arr2;
    use std::collections::HashMap;

    use super::CountMatrix;

    #[test]
    fn test_countmatrix() {
        let mut countmap: HashMap<(CB, GeneId), usize> = HashMap::new();
        countmap.insert((CB(0), GeneId(0)), 10);
        countmap.insert((CB(0), GeneId(1)), 1);
        countmap.insert((CB(1), GeneId(0)), 0); // lets see what happens with empty counts
        countmap.insert((CB(1), GeneId(1)), 5);

        let gene_vector = vec![Genename("geneA".to_string()), Genename("geneB".to_string())];

        let cmat = countmap_to_matrix(&countmap, gene_vector);

        let dense_mat = cmat.matrix.to_dense();
        let expected = arr2(&[[ 10, 1],
                              [ 0,  5]]);
        assert_eq!(dense_mat, expected);

        assert_eq!(
            cmat.cbs,
            vec![
                "AAAAAAAAAAAAAAAA".to_string(),
                "AAAAAAAAAAAAAAAC".to_string()
            ]
        );
    }

    #[test]
    fn test_read_write() {
        let mut countmap: HashMap<(CB, GeneId), usize> = HashMap::new();
        countmap.insert((CB(0), GeneId(0)), 10);
        countmap.insert((CB(0), GeneId(1)), 1);
        countmap.insert((CB(1), GeneId(0)), 0); // lets see what happens with empty counts
        countmap.insert((CB(1), GeneId(1)), 5);

        let gene_vector = vec![Genename("geneA".to_string()), Genename("geneB".to_string())];

        let cmat = countmap_to_matrix(&countmap, gene_vector);

        let path = "/tmp/rustbustools_test_read_write";
        if !std::path::Path::new(&path).exists() {
            std::fs::create_dir(path).unwrap();
        }

        cmat.write(path);

        let cmat2 = CountMatrix::from_disk(
            &format!("{}/gene.mtx", path),
            &format!("{}/gene.barcodes.txt", path),
            &format!("{}/gene.genes.txt", path),
        );

        assert!(cmat.is_equal(&cmat2));
    }

    #[test]
    fn test_countmatrix_equal() {
        //testing the .is_equal function, which should be order invariant (doesnt matter how genes are ordered)

        // two cells two genes
        // first cell has both genes
        // second cell has only the second gene
        let mut countmap1: HashMap<(CB, GeneId), usize> = HashMap::new();
        countmap1.insert((CB(0), GeneId(0)), 10);
        countmap1.insert((CB(0), GeneId(1)), 1);
        countmap1.insert((CB(1), GeneId(0)), 0); // lets see what happens with empty counts
        countmap1.insert((CB(1), GeneId(1)), 5);

        let gene_vector = vec![Genename("geneA".to_string()), Genename("geneB".to_string())];

        let cmat1 = countmap_to_matrix(&countmap1, gene_vector);

        // a version with permuated genes
        let mut countmap2: HashMap<(CB, GeneId), usize> = HashMap::new();
        countmap2.insert((CB(0), GeneId(1)), 10);
        countmap2.insert((CB(0), GeneId(0)), 1);
        countmap2.insert((CB(1), GeneId(1)), 0); // lets see what happens with empty counts
        countmap2.insert((CB(1), GeneId(0)), 5);

        let gene_vector = vec![Genename("geneB".to_string()), Genename("geneA".to_string())];
        let cmat2 = countmap_to_matrix(&countmap2, gene_vector);

        println!("{:?}", cmat1.to_map());
        println!("{:?}", cmat2.to_map());

        assert!(cmat1.is_equal(&cmat2))
    }
}
