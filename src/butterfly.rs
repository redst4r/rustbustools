//! `butterfly` provides quanification of amplification for a kallisto-bus scRNAseq experiment
//! 
//! # Introduction
//! In scRNAseq, each mRNA is tagged uniquely (up to random collision) with CB+UMI.
//! Those are then amplified and sequenced. 
//! If we see the same CB+UMI in multiple reads, we conclude that they are copies of the same original mRNA
//! For each unique mRNA we quantify its amplification factor and record the absolute 
//! frequency of the ampification (i.e. how often do we see a 5x amplification, 
//! 5reads for the same CB+UMI)
//! 
//! This module quantifies the amplifcation (very fast!). 
//! Further processing (where speed is not essential) is typically done in python,
//! e.g. saturation curves, unseen species estimators.
//! 
//! 
//! # Unseen species
//! Considering the CB+UMI as `species` and the reads as `observations`, this relates to the `unseen species` problem
//! How many unobserved `species` (CB+UMI) are there in the library given the amplification profile we've seen so far
//! While the module doesn't provide an unseen species estimator, it can easily be build on the [CUHistogram]
//! 
//! # References
//! The whole concept is described (amongst other things) in this 
//! [paper](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-021-02386-z)
//! 
//! # Examples
//! ```rust, no_run
//! # use rustbustools::io::BusFolder;
//! # use rustbustools::butterfly::make_ecs;
//! let b = BusFolder::new(
//!     "/path/to/bus/output",
//!     "/path/to/transcripts_to_genes.txt",
//! );
//! let h = make_ecs(&b, true);
//! // save the resulting frequency of frequency histogram to disk
//! // can be read in python for further processing (e.g. plot the saturation curves)
//! h.to_disk("/tmp/CU.csv")
//! ``` 

#![deny(missing_docs)]
use crate::{io::BusFolder, iterators::CbUmiGroupIterator, consistent_genes::find_consistent};
use std::{collections::HashMap, fs::File, io::Write};

/// The basic unit of this module, a frequency of frequency histogram
/// 
/// Records how many copies (reads) per mRNA (CB-UMI) we see in a busfile.
/// Should be constructed with [make_ecs]
#[derive(Debug)]
pub struct CUHistogram {
    // amplification (nReads for a single molecule) vs frequency
    histogram: HashMap<usize, usize>,
}
impl CUHistogram {
    // todo: useless!
    // fn new(h: HashMap<usize, usize>) -> Self {
    //     CUHistogram { histogram: h }
    // }

    /// return the number of reads (#molecules * number of copies) in the busfile
    pub fn get_nreads(&self) -> usize {
        self.histogram
            .iter()
            .map(|(nreads, freq)| nreads * freq)
            .sum()
    }

    /// return the number of molecules (distince CB/UMI pairs) in the busfile
    pub fn get_numis(&self) -> usize {
        self.histogram.values().sum()
    }

    /// calcualte the fraction of single-copy molecules (FSCM) in the busfile
    pub fn get_fscm(&self) -> f64 {
        let n1 = *self.histogram.get(&1).unwrap_or(&0);
        (n1 as f64) / (self.get_numis() as f64)
    }

    /// write the CU histogram into a csv on disk
    pub fn to_disk(&self, fname: &str) {
        let mut fh = File::create(fname).unwrap();

        fh.write("Amplification,Frequency\n".as_bytes()).unwrap();

        for (n_reads, n_umis) in self.histogram.iter() {
            fh.write_all(format!("{},{}\n", n_reads, n_umis).as_bytes())
                .unwrap();
        }
    }
}

/// Main function of this module: Quantities the amplification in the given busfolder
/// # Arguments
/// * `busfolder`: The folder containing the busfile, matric.ec etc...
/// * `collapse_ec`: How to handle identical CB-UMI but different EC:
///     - false: just ignore and lump the reads together irresepctive of EC
///     - true: check if they ECs are consistent (if yes, count as aggregate), if no, discard
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
        let c = CUHistogram { histogram: h };

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
