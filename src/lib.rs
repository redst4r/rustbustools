//! # bustools
//! 
//! This library allows interaction with the Bus format (see [bustools](https://github.com/BUStools/bustools)) 
//! for scRNAseq data processing. 
//! 
//! At this point, the package is pretty mature, but there might be some minor features missing compared to the original bustools.
//! 
//! # Basics of the library
//! The basic unit is the [`io::BusRecord`], which represents a single entry in a busfile,
//! consisting of CB, UMI, EC, COUNT and Flag.
//! 
//! [`io::BusReader`] and [`io::BusWriter`] are the primary means to actually read and write `.bus` files.
//! For compressed busfiles (`.busz`), use [`busz::BuszReader`] and [`busz::BuszWriter`] instead.
//! 
//! ## Iterate over a busfile
//! [`io`] contains the code to read and write from busfiles.
//! In particular it defines a simpe iterator over [`io::BusRecord`]s via [`io::BusReader`].
//! BusReader implements the trait [`io::CUGIterator`], a marker trait for anything that 
//! iterates/produced streams of [`io::BusRecord`]s in our library.
//! ```rust, no_run
//! # use bustools::io::BusReader;
//! let breader = BusReader::new("/path/to/some.bus");
//! for record in breader {
//!     // record.CB == ...
//! }
//! ```
//! 
//! ## Advanced Iterators over busfiles
//! While [`io::BusReader`] lets you iterate over single [`io::BusRecord`]s, 
//! it is often convenient to group the records by CB (all records from the same cell)
//! or by CB+UMI (all records from the same mRNA).
//! [`iterators`] contains the code to enable `chaining` iterators over BusRecords. 
//! 
//! Note that the bus file must be **sorted** (by CB/UMI) to enable these iterators (they will panic if used on an unsorted busfile).
//! 
//! ### Iterate over cells
//! To iterate over a *sorted* busfile, grouping all records by CB:
//! ```rust, no_run
//! # use bustools::io::BusReader;
//! use bustools::iterators::CellGroupIterator; //need to bring that trait into scope
//! 
//! let breader = BusReader::new("/path/to/some.bus");
//! for (cb, vector_of_records) in breader.groupby_cb() {
//!     // Example: the number of records in that cell
//!     let n_molecules: usize = vector_of_records.len();
//! }
//! ```
//! 
//! ### Iterate over molecules
//! To iterate over a **sorted** busfile, grouping all records by CB+UMI:
//! ```rust, no_run
//! # use bustools::io::BusReader; 
//! use bustools::iterators::CbUmiGroupIterator; //need to bring that trait into scope
//! 
//! let breader = BusReader::new("/path/to/some.bus");
//! for ((cb, umi), vector_of_records) in breader.groupby_cbumi() {
//!     // Example: the number of reads of that molecule (CB/UMI)
//!     let n_reads: u32 = vector_of_records.iter().map(|r| r.COUNT).sum();
//! }
//! ```
//! ## EC to gene mapping
//! More convenient features are provided by [`io::BusFolder`], 
//! which wraps around the `.bus` file, the `matric.ec` and `transcripts.txt` created by the `kallisto bus` command.
//! Those files tell us what a particular [`io::BusRecord`]  
//! actually maps to as specified by its EC (equivalence class, a set of transcripts).
//! This automatically constructs a mapper from equivalence class to gene via [`consistent_genes::Ec2GeneMapper`]
//! which allows to resolve ECs to genes.
//! 
//! ```rust, no_run
//! # use bustools::io::BusFolder;
//! # use bustools::consistent_genes::EC;
//! let bfolder = BusFolder::new("/path/to/busfolder");
//! let ec_mapper = bfolder.make_mapper("/path/to/transcripts_to_genes.txt");
//! let gene_names = ec_mapper.get_genenames(EC(1));
//! ```

// #![deny(missing_docs)]
pub mod io;
// pub mod io_dyn;
// pub mod io_generic;
pub mod iterators;

#[deprecated(note="please use `merger` instead", since="0.11.1")]
pub mod bus_multi;
pub mod utils;
pub mod consistent_genes;
pub mod disjoint;
pub mod merger;
pub mod busz;
// mod runlength_codec;
// pub mod channel;
// pub mod buffered_channels;
// pub mod new_channel;