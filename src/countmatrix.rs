use std::{
    collections::HashMap,
    fs::File,
    io::{BufRead, BufReader, Write},
};

use sprs::{
    io::{read_matrix_market, write_matrix_market},
    TriMat,
};

pub struct CountMatrix {
    /*
    Represents a count matrix of Cells vs Genes
    Cells are encoded a Strings already!
    */
    pub matrix: sprs::CsMat<i32>,
    pub cbs: Vec<String>,
    pub genes: Vec<String>,
}
impl CountMatrix {
    pub fn new(matrix: sprs::CsMat<i32>, cbs: Vec<String>, genes: Vec<String>) -> CountMatrix {
        CountMatrix { matrix, cbs, genes }
    }

    pub fn to_map(&self) -> HashMap<(String, String), i32> {
        // transforms the sparse count matrix into a Hashmap (CB,Gene)-> count
        let mut h1: HashMap<(String, String), i32> = HashMap::new();

        for (value, (i, j)) in self.matrix.iter() {
            h1.insert((self.cbs[i].clone(), self.genes[j].clone()), *value);
        }
        h1
    }

    pub fn is_equal(&self, other: &Self) -> bool {
        let h1 = self.to_map();
        let h2 = other.to_map();
        h1 == h2
    }

    pub fn from_disk(mtx_file: &str, cbfile: &str, genefile: &str) -> CountMatrix {
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
