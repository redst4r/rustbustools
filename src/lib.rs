//! # Rustbustools
//! 
//! This is a reimplementation of [bustools](https://github.com/BUStools/bustools) 
//! in rust for scRNAseq data processing. 
//! 
//! At this point, it's **far from complete and correct**, but rather a project to learn rust.
//! 
//! While this is intended as a CLI, there's some useful functionality in the library itself, 
//! iterating over bus files etc. See some [examples](#basics-of-the-library) below.
//! 
//! # CLI
//! `rustbustools <command>`
//! * `correct`: Correct busfiles via a whitelist
//! * `sort`: Sort the busfile by CB/UMI/EC
//! * `count`: Create a count-matrix (CB vs gene)
//! * `inspect`: Basic stats about a busfile (#records, #CBs etc..)
//! 
//! Check the CLI help for arguments.
//! 
//! # Basics of the library
//! The basic unit is the [io::BusRecord], which represents a single entry in a busfile,
//! consisting of CB, UMI, EC, COUNT and Flag.
//! 
//! ## Iterate over a busfile
//! [io] contains the code to read and write from busfiles.
//! In particular it defines a simpe iterator over [io::BusRecord]s via [io::BusReader].
//! BusReader implements the trait [io::CUGIterator], a marker trait for anything that 
//! iterates/produced streams of [io::BusRecord]s in our library.
//! ```rust, no_run
//! # use rustbustools::io::BusReader;
//! let breader = BusReader::new("/path/to/some.bus");
//! for record in breader {
//!     // record.CB == ...
//! }
//! ```
//! 
//! ## Advanced Iterators over busfiles
//! While [io::BusReader] lets you iterate over single [io::BusRecord]s, 
//! it is often convenient to group the records by CB (all records from the same cell)
//! or by CB+UMI (all records from the same mRNA).
//! [iterators] contains the code to enable `chaining` iterators over BusRecords. 
//! 
//! Note that the bus file must be **sorted** (by CB/UMI) to enable these iterators (they will panic if used on an unsorted busfile).
//! 
//! ### Iterate over cells
//! To iterate over a *sorted* busfile, grouping all records by CB:
//! ```rust, no_run
//! # use rustbustools::io::BusReader;
//! use rustbustools::iterators::CellGroupIterator; //need to bring that trait into scope
//! 
//! let breader = BusReader::new("/path/to/some.bus");
//! for (cb, vector_of_records) in breader.groupby_cb() {
//!     // Example: the number of records in that cell
//!     let n_molecules: usize = vector_of_records.len();
//!     
//! }
//! ```
//! 
//! ### Iterate over molecules
//! To iterate over a `sorted` busfile, grouping all records by CB+UMI:
//! ```rust, no_run
//! # use rustbustools::io::BusReader; 
//! use rustbustools::iterators::CbUmiGroupIterator; //need to bring that trait into scope
//! 
//! let breader = BusReader::new("/path/to/some.bus");
//! for ((cb, umi), vector_of_records) in breader.groupby_cbumi() {
//!     // Example: the number of reads of that molecule (CB/UMI)
//!     let n_reads: u32 = vector_of_records.iter().map(|r| r.COUNT).sum();
//! }
//! ```
//! ## EC to gene mapping
//! More convenient features are provided by [io::BusFolder], 
//! which wraps around the busfile, the matric.ec and transcripts.txt created by `kallisto bus`.
//! Those files tell us what a particular bus record `(CB,UMI,EC,Count,flag)` 
//! actually maps to as specified by its EC (equivalence class, a set of transcripts).
//! This automatically constructs a mapper from equivalence class to gene via [consistent_genes::Ec2GeneMapper]
//! which allows to resolve ECs to genes
//! 
//! ```rust, no_run
//! # use rustbustools::io::BusFolder;
//! # use rustbustools::consistent_genes::EC;
//! 
//! let bfolder = BusFolder::new("/path/to/busfolder", "/path/to/transcripts_to_genes.txt");
//! let gene_names = bfolder.ec2gene.get_genenames(EC(1));
//! ```
// #![deny(missing_docs)]
pub mod io;
pub mod iterators;
// pub mod iterators_clone;
// pub mod play;
pub mod busmerger;
pub mod bus_multi;
pub mod count;
pub mod count2;
pub mod utils;
pub mod multinomial;
pub mod consistent_genes;
pub mod disjoint;
pub mod butterfly;
pub mod inspect;
pub mod countmatrix;
// pub mod channel;
// pub mod buffered_channels;
// pub mod new_channel;
pub mod merger;
pub mod sort;
pub mod correct;