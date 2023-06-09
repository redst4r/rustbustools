use crate::{io::BusFolder, iterators::CbUmiGroupIterator, consistent_genes::find_consistent};
use std::{collections::HashMap, fs::File, io::Write};

#[derive(Debug)]
pub struct CUHistogram {
    // amplification (nReads for a single molecule) vs frequency
    histogram: HashMap<usize, usize>,
}
impl CUHistogram {
    pub fn new(h: HashMap<usize, usize>) -> Self {
        CUHistogram { histogram: h }
    }
    /// return the number of reads (#molecules * number of copies) in the busfile
    pub fn get_nreads(&self) -> usize {
        self.histogram
            .iter()
            .map(|(nreads, freq)| nreads * freq)
            .sum()
    }

    /// return the number of molecules (distince CB/UMI pairs) in the busfile
    pub fn get_numis(&self) -> usize {
        self.histogram.iter().map(|(_nreads, freq)| freq).sum()
    }

    /// calcualte the fraction of single-copy molecules (FSCM) in the busfile
    pub fn get_fscm(&self) -> f64 {
        let n1 = *self.histogram.get(&1).unwrap_or(&0);
        (n1 as f64) / (self.get_numis() as f64)
    }

    /// write the CU histogram into a csv on disk
    pub fn to_disk(&self, fname: &str) {
        let mut fh = File::create(fname).unwrap();
        for (n_reads, n_umis) in self.histogram.iter() {
            fh.write_all(format!("{},{}\n", n_reads, n_umis).as_bytes())
                .unwrap();
        }
    }
}

pub fn make_ecs(busfolder: &BusFolder, collapse_ec: bool) -> CUHistogram {

    let mut h: HashMap<usize, usize> = HashMap::new();

    let mut multimapped = 0;
    let mut inconsistent = 0;
    let mut total = 0;

    for ((_cb, _umi), recordlist) in busfolder.get_iterator().groupby_cbumi() {
        total+=1;
        if collapse_ec {
            // check if we can uniquely match those read to the same gene
            // if not its either multimapped or inconsistent (could be a CB/UMI collision)
            match find_consistent(&recordlist, &busfolder.ec2gene) {
                crate::consistent_genes::MappingResult::SingleGene(_) => {
                    // increment our histogram
                    let nreads: usize = recordlist.iter().map(|x| x.COUNT as usize).sum();
                    let freq = h.entry(nreads).or_insert(0);
                    *freq += 1;
                },
                crate::consistent_genes::MappingResult::Multimapped(_) => {multimapped+=1},
                crate::consistent_genes::MappingResult::Inconsistent => {inconsistent+=1},
            }
        } else {
            let nreads: usize = recordlist.iter().map(|x| x.COUNT as usize).sum();
            let freq = h.entry(nreads).or_insert(0);
            *freq += 1;
        }
    }

    println!("Total CB-UMI {}, Multimapped {} ({}%), Discarded/Inconsistent {} ({}%)", total, multimapped, (multimapped as f32)/ (total as f32), inconsistent, (inconsistent as f32)/ (total as f32));
    CUHistogram { histogram: h }
}

#[cfg(test)]
mod testing {
    use crate::{
        butterfly::{make_ecs, CUHistogram},
        io::BusFolder,
    };
    use statrs::assert_almost_eq;
    use std::collections::HashMap;

    #[test]
    pub fn testing() {
        let h: HashMap<usize, usize> = vec![(1, 2), (3, 3)].into_iter().collect();
        let c = CUHistogram::new(h);

        assert_eq!(c.get_nreads(), 11);
        assert_eq!(c.get_numis(), 5);
        assert_almost_eq!(c.get_fscm(), 2.0 / 5.0, 0.00000000000000001);
    }

    #[test]
    pub fn testing2() {
        // let s = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus".to_string();
        let b = BusFolder::new(
            "/home/michi/bus_testing/bus_output_short",
            "/home/michi/bus_testing/transcripts_to_genes.txt",
        );
        let h = make_ecs(&b, true);
        println!("{:?}", h);
    }
}
