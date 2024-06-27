//! The io module of bustools deals with reading and writing busfiles.
//!
//! The most important components are:
//! - [`BusRecord`]: a single record of a Busfile, containing CB/UMI/EC/COUNT
//! - [`BusFolder`]: representing the directory structure of kallisto quantification
//! - [`BusReader`]: to iterate over Records of a busfile
//! - [`BusWriter`]: writing records to busfile
//!
//! # Example
//! ```rust, no_run
//! # use bustools::io::{BusReader, BusParams, BusWriter, BusRecord};
//! let bus = BusReader::new("/tmp/some.bus");
//! let params = BusParams {cb_len: 16, umi_len: 12};
//! let mut writer = BusWriter::new("/tmp/out.bus", params);
//! let records: Vec<BusRecord> = bus.collect();
//! writer.write_iterator(records.into_iter());
//!
#![allow(non_snake_case)]

use crate::busz::{BuszReader, BuszWriter};
use crate::consistent_genes::{Ec2GeneMapper, EC, make_mapper};
use crate::consistent_transcripts::{make_mapper_transcript, Ec2TranscriptMapper, TranscriptId, Transcriptname};
use crate::iterators::{CbUmiGroupIterator, CellGroupIterator};
use bincode;
use serde;
use std::collections::HashMap;
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Write};
use tempfile::TempDir;
use rkyv::{self, Deserialize};

pub (crate) const BUS_ENTRY_SIZE: usize = 32;
pub const BUS_HEADER_SIZE: usize = 20;

/// Basic unit of a busfile as created by kallisto.
/// Represents a single scRNAseq molecule, CB/UMI/EC and COUNT
///
/// from python implemenation/specification
/// BUS_ENTRY_SIZE = 32  # each entry is 32 bytes!!
///  unpack_str = 'QQiIIxxxx'
/// Q: 8byte unsigned long,long int
/// i: 4byte int
/// I: unsigned int, 4byte
#[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq, Eq, Clone, rkyv::Archive, rkyv::Deserialize, rkyv::Serialize)] // TODO take away the copy and clone
#[archive(check_bytes)]
pub struct BusRecord {
    pub CB: u64,    // 8byte
    pub UMI: u64,   // 8byte
    pub EC: u32,    // 4v byte
    pub COUNT: u32, // 4v byte
    pub FLAG: u32,  // 4v byte    including this, we have 28 bytes, missing 4 to fill up
    // PAD: u32     // just padding to fill 32bytes
}
impl BusRecord {
    /// converts a busrecord to byte representation, e.g. to write the bytes to a busfile
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut binrecord = bincode::serialize(self).expect("FAILED to serialze record");
        // let mut binrecord = bincode::serde::encode_to_vec(self, bincode::config::legacy()).unwrap(); //.expect("FAILED to deserialze header");

        // the struct is only 28bytes, so we need 4 padding bytes
        for _i in 0..4 {
            //TODO feel like theres a better way to do this
            binrecord.push(0);
        }
        binrecord
    }

    /// coverts bytes into burecord, e.g. when reading a busfile
    pub fn from_bytes(bytes: &[u8]) -> Self {
        // bincode::serde::decode_from_slice(bytes, bincode::config::legacy()).expect("FAILED to deserialze record").0 //.expect("FAILED to deserialze header");

        // let archived = rkyv::check_archived_root::<BusRecord>(&bytes[..]).unwrap();
        let archived = unsafe { rkyv::archived_root::<BusRecord>(bytes) };
        let s:BusRecord = archived.deserialize(&mut rkyv::Infallible).unwrap();
        s
        // bincode::deserialize(bytes).expect("deserial error")
    }
}

#[derive(Debug, PartialEq, Eq, Clone)]
pub struct BusParams {
    pub cb_len: u32,
    pub umi_len: u32,    
}

/// Header of a busfile, as specified by kallisto
/// the only things that are variable are:
/// - cb_len: The number of bases in the cellbarcode (usually 16BP)
/// - umi_len: The number of bases in the cellbarcode (usually 12BP)
/// # Example
/// ```
/// // use bustools::io::BusHeader;
/// // let header = BusHeader::new(16, 12, 1);
/// // can also be obtained from an existing busfile
/// // let header = BusHeader::from_file("somefile.bus");
/// ```
#[derive(serde::Serialize, serde::Deserialize, Debug, PartialEq, Eq, Clone)]
pub struct BusHeader {
    //4sIIII: 20bytes
    pub(crate) magic: [u8; 4],
    pub(crate) version: u32,
    pub(crate) cb_len: u32,
    pub(crate) umi_len: u32,
    pub(crate) tlen: u32,
}

impl BusHeader {
    /// construct a busheader from CB/UMI length
    pub fn new(cb_len: u32, umi_len: u32, tlen: u32, compressed: bool) -> BusHeader {
        let magic: [u8; 4] = if compressed {
            *b"BUS\x01"
        } else {
            *b"BUS\x00"
        };
        BusHeader { cb_len, umi_len, tlen, magic, version: 1 }
    }

    /// creates the header struct, extracting it from an existing bus file
    pub fn from_file(fname: &str) -> BusHeader {
        // getting 20 bytes out of the file, which is the header
        let file = std::fs::File::open(fname).unwrap_or_else(|_| panic!("file not found: {fname}"));
        let mut header = Vec::with_capacity(BUS_HEADER_SIZE);
        let _n = file
            .take(BUS_HEADER_SIZE as u64)
            .read_to_end(&mut header)
            .unwrap();
        BusHeader::from_bytes(&header)
    }

    /// desearializes a BusHeader from Bytes; when reading busfiles
    pub fn from_bytes(bytes: &[u8]) -> BusHeader {
        let header_struct: BusHeader =
            bincode::deserialize(bytes).expect("FAILED to deserialze header");
        header_struct
    }

    /// seialize the header to bytes
    pub fn to_bytes(&self) -> Vec<u8> {
        bincode::serialize(self).expect("FAILED to serialze header")
    }

    /// return the length of the variable part of the busheader
    pub fn get_tlen(&self) -> u32 {
        self.tlen
    }

    // /// return the length of the Cell Barcode in the busfile
    // pub fn get_cb_len(&self) -> u32 {
    //     self.cb_len
    // }
}

/// A marker trait for iterators over CB/UMI/gene_EC iterators
/// to create an iterator compatible with this framework,
/// just implement the usual iterator trait (giving the next() function)
/// and also tag the iterator with CUGIterator: `impl CUGIterator for Dummy {}`
/// ## Example
/// ```
/// struct Dummy{
///     f: u32
/// }
/// use bustools::io::{BusRecord, CUGIterator};
/// impl Iterator for Dummy {
///     type Item=BusRecord;
///
///     fn next(&mut self) -> Option<Self::Item> {
///         Some(BusRecord{CB: 0,UMI: 0, EC: 0, COUNT: self.f, FLAG:0})
///     }
/// }
/// impl CUGIterator for Dummy {}
/// let mut x =  Dummy{f:1};
/// let s = x.next();
/// println!("{s:?}")
/// ```
///
/// pretty much a type alias for Iterator<Item=BusRecord>
pub trait CUGIterator: Iterator<Item = BusRecord> {}
pub const DEFAULT_BUF_SIZE: usize = 800 * 1024; // 800  KB


















/// The main struct to read (un)compressed bus files.
/// Should always be used in favour of the concrete [`BusReaderPlain`] (plain .bus files) and [`BuszReader`] (compressed .busz files),
/// as this allows code to be transparent/agnositc of wether the underlying file is compressed or not. 
/// 
/// # Example
/// Uncompressed busfile
/// ```rust, no_run
/// # use bustools::io::BusReader;
/// let breader = BusReader::new("somefile.bus");
/// for record in breader{
///     let cb= record.CB;
/// };
/// ```
/// Compressed busfile (works just the same way, note the different file extension)
/// ```rust, no_run
/// # use bustools::io::BusReader;
/// let breader = BusReader::new("somefile.busz");
/// for record in breader{
///     let cb= record.CB;
/// };
/// ```
pub enum BusReader <'a> {
    Plain(BusReaderPlain<'a>),
    Compressed(BuszReader<'a>),
}

impl <'a>BusReader<'a> {
    /// Creates a new reader of a bus/busz file. Depending on the file extension,
    /// it infers wether the file is plain-binary (.bus) or compressed (.busz).
    /// 
    /// # TODO: make the plain/compressed choice explicit?
    pub fn new(fname: &str) -> Self {
        if fname.ends_with(".busz") {
            Self::new_compressed(fname)
        } else {
            Self::new_plain(fname)
        }
    } 

    /// construct the BusReader on any input stream which yields bytes (implements read)
    /// if this is based on a File, highly recommended to wrap it in a BufferedReader for performance
    /// otherwise every iteration will cause a read/system call
    pub fn from_read_plain(reader: impl Read + 'a) -> Self {
        BusReader::Plain(BusReaderPlain::from_read(reader))
    }

    pub fn from_read_compressed(reader: impl Read + 'a) -> Self {
        BusReader::Compressed(BuszReader::from_read(reader))
    }

    // THIS ONE IS DEPRECATED: IF YOU WANT SEPCIFIC BUFFER SIZE, FEED a BufferedReader into :from_read_*
    // /// Creates a new reader with specified buffer size (800KB = 800 * 1024) is
    // /// a good choice. 
    // pub fn new_with_capacity(fname: &str, bufsize: usize) -> Self {
    //     if fname.ends_with(".busz") {
    //         BusReader::Compressed(BuszReader::new_with_capacity(fname, bufsize))
    //     } else {
    //         BusReader::Plain(BusReaderPlain::new_with_capacity(fname, bufsize))
    //     }
    // } 

    /// specifically return a Reader for plain bus files
    pub fn new_plain(fname: &str) -> Self{
        BusReader::Plain(BusReaderPlain::new(fname))
    }

    /// specifically return a Reader for compressed busfiles
    pub fn new_compressed(fname: &str) -> Self{
        BusReader::Compressed(BuszReader::new(fname))
    }

    // /// get the [`BusParams`] of the file.
    pub fn get_params(&self) -> &BusParams { 
        match self {
            BusReader::Plain(reader) => reader.get_params(),
            BusReader::Compressed(reader) => reader.get_params(),
        }
    }
}

impl <'a> Iterator for BusReader<'a> {
    type Item = BusRecord;
    fn next(&mut self) -> Option<Self::Item> {
        match self {
            BusReader::Plain(reader) => reader.next(),
            BusReader::Compressed(reader) => reader.next(),
        }
    }
}

impl <'a> CUGIterator for BusReader <'a> {}

//=================================================================================================




/// Main reader for plain Busfiles
/// Allows to iterate over the BusRecords in the file
///
/// # Example
/// ```rust, no_run
/// # use bustools::io::BusReaderPlain;
/// let breader = BusReaderPlain::new("somefile.bus");
/// for record in breader{
///     let cb= record.CB;
/// };
/// ```
/// # From `Read`
/// The [`BusReaderPlain`] can operate on anything implementing the `Read`-trait.
/// For example, one can create a [`BusReaderPlain`] from a [`File`]:
/// ```rust, no_run
/// # use bustools::io::BusReaderPlain;
/// # use std::io::BufReader;
/// # use std::fs::File;
/// // note: Buffering is highly recommended for performance reasons
/// let bufReader = BufReader::with_capacity(10000, File::open("somefile.bus").unwrap());
/// let breader = BusReaderPlain::from_read(bufReader);
/// ```
/// 
/// Similarly one can also construct a [`BusReaderPlain`] for in-memory data:
/// ```rust
/// # use bustools::io::{BusReaderPlain, BusRecord, BusHeader};
/// // create an in memory representation of a busfile (in bytes)
/// let header = BusHeader::new(16, 12, 1, false);
/// let r1 = BusRecord{CB: 1, UMI:1, EC:1, COUNT: 10, FLAG: 0};
/// let mut v = Vec::new();
/// v.extend_from_slice(&header.to_bytes());
/// v.extend_from_slice(&r1.to_bytes());
/// 
/// // v contains a busfile as a Vec<u8> which we can read with BusReader
/// let breader = BusReaderPlain::from_read(&v[..]);
/// let records:Vec<_> = breader.collect();
/// ```
pub struct BusReaderPlain<'a> {
    /// containing info about CB-length, UMI length for decoding
    pub params: BusParams,
    reader:  Box<dyn Read+ 'a>, //ugly way to store any Read-like object in here, BufferedReader, File, Cursor, or just a vec<u8>!
}
impl <'a> BusReaderPlain <'a> {

    /// Creates a BusReader for a file on disk.
    /// Turns the file into a bufferedReader.
    pub fn new(fname: &str) -> Self {
        let file_handle = File::open(fname).expect("Busfile should exist on disk!");       
        let buf = BufReader::with_capacity(DEFAULT_BUF_SIZE, file_handle);
        Self::from_read(buf)
    }

    /// construct the BusReader on any input stream which yields bytes (implements read)
    /// if this is based on a File, highly recommended to wrap it in a BufferedReader for performance
    /// otherwise every iteration will cause a read/system call
    pub fn from_read(mut reader: impl Read + 'a) -> Self {

        // parse header
        let mut header_bytes = [0_u8; BUS_HEADER_SIZE];
        reader.read_exact(&mut header_bytes).expect("failed to read header");
        let header = BusHeader::from_bytes(&header_bytes);
        let params = BusParams { cb_len: header.cb_len, umi_len: header.umi_len };
        
        assert_eq!(
            &header.magic, b"BUS\x00",
            "Header struct not matching; MAGIC is wrong"
        );

        // the variable header
        let mut buffer = Vec::with_capacity(header.tlen as usize);
        for _i in 0..header.tlen {
            buffer.push(0_u8);
        }
        reader.read_exact(&mut buffer).expect("failed to read variable header");

        // reader is now positioned at the first record
        BusReaderPlain { params, reader: Box::new(reader) }
    }
    
    pub fn get_params(&self) -> &BusParams {
        &self.params
    }
}

impl <'a> Iterator for BusReaderPlain <'a> {
    type Item = BusRecord;

    fn next(&mut self) -> Option<Self::Item> {
        /* 
        The previous version using only `read()` had a bit of a bug:
        from unknown reasons (but consitent with the `read() documentation)
        `read` will not fill the entire buffer, which would case the previous implementation
        to throw an error.
        We can instead use `read_exact`. Drawback: 
        This is a bit of a dangerous version: We dont check what happens at the EOF. 
        If the file contains a truncated last record, we would just skip over that 

        However: It's much faster than the above take/read_to_end (~5x)
        */
        let mut buffer = [0; BUS_ENTRY_SIZE]; // TODO this gets allocated at each next(). We could move it out into the struct: Actually doesnt make a diference, bottleneck is the reading
        match self.reader.read_exact(&mut buffer) {
            Ok(()) => Some(BusRecord::from_bytes(&buffer)),
            Err(e) => match e.kind() {
                // this can happen due to a real EOF, or we're almost at the EOF 
                // and try to read more bytes than whats left (i.e. a truncated file)
                std::io::ErrorKind::UnexpectedEof => None,
                _ => panic!("{e}"),
            }
        }
    }
}

// tag our iterator to be compatible with our framework
impl <'a> CUGIterator for BusReaderPlain<'a> {}


/// actually, also make general iterators of BusRecord amenable
// impl CUGIterator for dyn Iterator<Item = BusRecord> {}
impl CUGIterator for std::vec::IntoIter<BusRecord> {}





/// Main writer for bus records, either uncompressed or compressed.
pub enum BusWriter {
    Plain(BusWriterPlain),
    Compressed(BuszWriter),
}

impl BusWriter {
    pub fn new(filename: &str, params: BusParams) -> Self {
        if filename.ends_with(".busz") {
            Self::Compressed(
                // BuszWriter::new_with_capacity(file_handle, header, DEFAULT_BUF_SIZE)
                BuszWriter::new(filename, params, 100_000)  // TODO: buz_block
            )
        } else {
            Self::Plain(
                BusWriterPlain::new(filename, params)
            )
        }
    }

    pub fn write_iterator(&mut self, iter: impl Iterator<Item=BusRecord>) {
        match self {
            BusWriter::Plain(writer) => writer.write_records(&iter.collect()), // TODO remove allocation
            BusWriter::Compressed(writer) => writer.write_iterator(iter),
        }
    }
    // pub fn new_with_capacity(file_handle: File, header: BusHeader, bufsize: usize) -> Self {


    //     if file_handle.ends_with(".busz") {
    //         Self::Compressed(
    //             BuszWriter::new_with_capacity(file_handle, header, 100_000)
    //         )     
    //     } else {
    //         Self::Plain(
    //             BusWriterPlain::new_with_capacity(file_handle, header, bufsize)
    //         )  
    //     }
    // }
}







// ========================================
/// Writing BusRecords into a File, using Buffers
/// needs the Header of the busfile to be specified (length of CB and UMI)
/// # Example
/// ```
/// # use bustools::io::{BusWriterPlain, BusParams, BusRecord};
///
/// let params = BusParams { cb_len: 16, umi_len: 12};
/// let mut w = BusWriterPlain::new("/tmp/target.bus", params);
/// let record = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 10, FLAG: 0};
/// w.write_record(&record);
/// ```
pub struct BusWriterPlain {
    writer: BufWriter<File>,
    header: BusHeader,
}

impl BusWriterPlain {
    /// create a BusWriter from an open Filehandle
    pub fn from_filehandle(file_handle: File, params: BusParams) -> BusWriterPlain {
        BusWriterPlain::new_with_capacity(file_handle, params, DEFAULT_BUF_SIZE)
    }

    /// create a Buswriter that streams records into a file
    pub fn new(filename: &str, params: BusParams) -> BusWriterPlain {
        let file_handle: File = File::create(filename).expect("FAILED to open");
        BusWriterPlain::from_filehandle(file_handle, params)
    }

    /// Writing a single BusRecord
    pub fn write_record(&mut self, record: &BusRecord) {
        let binrecord = record.to_bytes();
        self.writer
            .write_all(&binrecord)
            .expect("FAILED to write record");
    }

    /// Writing multiple BusRecords at once
    pub fn write_records(&mut self, records: &Vec<BusRecord>) {
        // writes several recordsd and flushes
        for r in records {
            self.write_record(r)
        }
        self.writer.flush().unwrap();
    }

    /// creates a buswriter with specified buffer capacity (after how many records an actual write happens)
    /// dont use , 800KB is the default buffer size and optimal for performance
    pub fn new_with_capacity(file_handle: File, params: BusParams, bufsize: usize) -> Self {
        let mut writer = BufWriter::with_capacity(bufsize, file_handle);

        let custom_header_str = "BUS file produced by kallisto".as_bytes();
        let header = BusHeader::new(params.cb_len, params.umi_len, custom_header_str.len() as u32, false);

        // write the header into the file
        let binheader = header.to_bytes();
        writer
            .write_all(&binheader)
            .expect("FAILED to write header");

        let varheader = custom_header_str;

        writer
            .write_all(varheader)
            .expect("FAILED to write var header");

        BusWriterPlain { writer, header }
    }

    pub fn write_iterator(&mut self, iter: impl Iterator<Item=BusRecord>) {
        for r in iter {
            self.write_record(&r)
        }
        self.writer.flush().unwrap();
    }
}
//=================================================================================

/// Represents the standard file structure of kallisto quantification:
/// - a busfile
/// - and ec.matrix file, containing the mapping between EC and transcript_id
/// - transcripts.txt, containing the ENST transcript id
///
/// makes life easier when dealing with kallisto quantifications and the mapping of EC to Gene.
/// To construct it, one has to supply a transcript_to_gene mapping file. This allows ECs to be mapped to a set of genes via
/// the Ec2GeneMapper
pub struct BusFolder {
    // pub foldername: String,
    busfile: String,
    matrix_ec_file: String,
    transcript_file: String,

}

pub fn parse_ecmatrix(filename: &str) -> HashMap<EC, Vec<TranscriptId>> {
    // parsing an ec.matrix into a Hashmap EC->list_of_transcript_ids
    let file = File::open(filename).unwrap_or_else(|_| panic!("{} not found", filename));
    let reader = BufReader::new(file);
    let mut ec_dict: HashMap<EC, Vec<TranscriptId>> = HashMap::new();

    for line in reader.lines() {
        if let Ok(l) = line {
            // println!("{}", l);
            let mut s = l.split_whitespace();
            let ec = s.next().unwrap().parse::<u32>().unwrap();
            let tmp = s.next().unwrap();

            let transcript_list: Vec<TranscriptId> =
                tmp.split(',').map(|x| TranscriptId(x.parse::<u32>().unwrap())).collect();
            ec_dict.insert(EC(ec), transcript_list);
        } else {
            panic!("Error reading file {}", filename)
        }
    }
    ec_dict
}

fn parse_transcript(filename: &str) -> HashMap<TranscriptId, Transcriptname> {
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut transcript_dict: HashMap<TranscriptId, Transcriptname> = HashMap::new();
    for (i, line) in reader.lines().enumerate() {
        if let Ok(l) = line {
            transcript_dict.insert(TranscriptId(i as u32), Transcriptname(l));
        }
    }
    transcript_dict
}

impl BusFolder {
    pub fn new(foldername: &str) -> BusFolder {

        // the default names
        let busfile = format!("{foldername}/output.corrected.sort.bus");
        let matrix_ec_file = format!("{foldername}/matrix.ec");
        let transcript_file = format!("{foldername}/transcripts.txt");

        BusFolder { busfile, matrix_ec_file, transcript_file}
    }

    pub fn from_files(busfile: &str, matrix_ec_file: &str, transcript_file: &str) -> Self{
        BusFolder { busfile: busfile.to_string(), matrix_ec_file: matrix_ec_file.to_string(), transcript_file: transcript_file.to_string()}
    }
 
    /// returns an iterator of the folder's busfile
    pub fn get_iterator(&self) -> BusReader {
        let bfile = self.get_busfile();
        BusReader::new(&bfile)
    }

    /// return the folders busfile
    pub fn get_busfile(&self) -> String {
        self.busfile.clone()
    }

    /// return the busfiles header
    // pub fn get_bus_header(&self) -> BusHeader {
    //     BusHeader::from_file(&self.get_busfile())
    // }

    pub fn get_bus_params(&self) -> BusParams {
        let h = BusHeader::from_file(&self.get_busfile());
        BusParams {cb_len:h.cb_len, umi_len: h.umi_len}
    }

    /// return the matric.ec file
    pub fn get_ecmatrix_file(&self) -> String {
        self.matrix_ec_file.clone()
    }

    /// return the transcript file
    pub fn get_transcript_file(&self) -> String {
        self.transcript_file.clone()
    }

    /// return the matric.ec file
    pub fn parse_ecmatrix(&self) -> HashMap<EC, Vec<TranscriptId>> {
        let filename = self.get_ecmatrix_file();
        parse_ecmatrix(&filename)
    }

    pub fn parse_transcript(&self) -> HashMap<TranscriptId, Transcriptname> {
        let filename = self.get_transcript_file();
        parse_transcript(&filename)
    }

    pub fn get_cb_size(&self) -> usize {
        let cb_iter_tmp = self.get_iterator().groupby_cb();
        cb_iter_tmp.count()
    }

    pub fn get_cbumi_size(&self) -> usize {
        let cb_iter_tmp = self.get_iterator().groupby_cbumi();
        cb_iter_tmp.count()
    }

    pub fn make_mapper(&self, t2g_file: &str) -> Ec2GeneMapper {
        make_mapper(self, t2g_file)
    }
    pub fn make_mapper_transcript(&self) -> Ec2TranscriptMapper {
        make_mapper_transcript(self)
    }
}

/// group a list of records by their CB/UMI (i.e. a single molecule)
///
/// takes ownership of record_list, since it reorders the elements
pub fn group_record_by_cb_umi(record_list: Vec<BusRecord>) -> HashMap<(u64, u64), Vec<BusRecord>> {
    let mut cbumi_map: HashMap<(u64, u64), Vec<BusRecord>> = HashMap::new();

    for r in record_list {
        let rlist = cbumi_map.entry((r.CB, r.UMI)).or_default(); // inserts a Vec::new if not present
        rlist.push(r);
    }
    cbumi_map
}

pub fn setup_busfile(records: &Vec<BusRecord>) -> (String, TempDir) {
    // just writes the records into a temporay file
    // returns the filename
    // returns TempDir so that it doesnt go out of scope and gets deleted right after tis function returns
    use tempfile::tempdir;
    let dir = tempdir().unwrap();
    let file_path = dir.path().join("output.corrected.sort.bus");
    let tmpfilename = file_path.to_str().unwrap();

    let mut writer = BusWriterPlain::new(tmpfilename, BusParams {cb_len: 16, umi_len: 12});
    writer.write_records(records);

    (tmpfilename.to_string(), dir)
}

pub fn write_partial_busfile(bfile: &str, boutfile: &str, nrecords: usize) {
    // write the first nrecords of the intput file into the output
    let busiter = BusReaderPlain::new(bfile);
    let mut buswriter = BusWriterPlain::new(boutfile, busiter.params.clone());

    for record in busiter.take(nrecords) {
        buswriter.write_record(&record);
    }
}

//=================================================================================
#[cfg(test)]
mod tests {
    use crate::consistent_genes::EC;
    use crate::consistent_transcripts::TranscriptId;
    use crate::io::{setup_busfile, BusHeader, BusReaderPlain, BusRecord, BusWriterPlain, BusParams};
    use crate::iterators::CellGroupIterator;
    use std::io::{BufReader, Read, Write};
    use tempfile::tempdir;

    #[test]
    fn test_read_write_header() {
        let r1 = BusRecord { CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };

        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test_read_write_header.bus");
        let busname = file_path.to_str().unwrap();

        let mut writer = BusWriterPlain::new(busname, BusParams {cb_len: 16, umi_len: 12});
        writer.write_record(&r1);
        writer.writer.flush().unwrap();

        let bheader = BusHeader::from_file(&busname);
        let header = BusHeader::new(16, 12, 20, false);
        // let params = BusParams{cb_len: 16, umi_len:12 };

        assert_eq!(header.magic, bheader.magic);
        assert_eq!(header.cb_len, bheader.cb_len);
        assert_eq!(header.umi_len, bheader.umi_len);
        assert_eq!(header.version, bheader.version);
        // Note: tlen is not preserved!!

    }

    #[test]
    fn test_read_write() {
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 1, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let rlist = vec![r1, r2]; // note: this clones r1, r2!
        let (busname, _dir) = setup_busfile(&rlist);
        let reader = BusReaderPlain::new(&busname);

        let records: Vec<BusRecord> = reader.into_iter().collect();
        assert_eq!(records, rlist)
    }

    use std::fs::File;
    use std::io::BufWriter;

    use super::parse_ecmatrix;
    #[test]
    fn test_ecmatrix() {

        let dir = tempdir().unwrap();
        let file_path = dir.path().join("foo.txt");
        let tmpfilename = file_path.to_str().unwrap();

        let file = File::create(tmpfilename).expect("Unable to create file");
        let mut f = BufWriter::new(file);

        let data = "0 1,2,3\n1 3,4\n2 4";
        f.write_all(data.as_bytes()).expect("Unable to write data");
        f.flush().unwrap();

        let ec = parse_ecmatrix(tmpfilename);
        // println!("{}", ec.len());
        // println!("{:?}", ec);
        let e1 = ec.get(&EC(0)).unwrap();
        assert_eq!(*e1, vec![TranscriptId(1), TranscriptId(2), TranscriptId(3)]);

        let e1 = ec.get(&EC(1)).unwrap();
        assert_eq!(*e1, vec![TranscriptId(3), TranscriptId(4)]);
    }

    use super::group_record_by_cb_umi;
    #[test]
    fn test_grouping() {
        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
        let r3 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r4 = BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
        let r5 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0 };
        let r6 = BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };

        let records = vec![
            r1.clone(),
            r2.clone(),
            r3.clone(),
            r4.clone(),
            r5.clone(),
            r6.clone(),
        ];

        let grouped = group_record_by_cb_umi(records);

        let g1 = grouped.get(&(0 as u64, 1 as u64)).unwrap();
        assert_eq!(*g1, vec![r1, r2]);

        let g2 = grouped.get(&(0 as u64, 2 as u64)).unwrap();
        assert_eq!(*g2, vec![r3]);

        let g3 = grouped.get(&(1 as u64, 1 as u64)).unwrap();
        assert_eq!(*g3, vec![r4, r6]);

        let g4 = grouped.get(&(1 as u64, 2 as u64)).unwrap();
        assert_eq!(*g4, vec![r5]);
    }

    #[test]
    /// make sure the format on disk stays the same when we write it
    fn test_insta_binary_format_write() {
        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.bus");
        let tmpfilename = file_path.to_str().unwrap();

        let mut wri = BusWriterPlain::new(tmpfilename, BusParams {cb_len: 16, umi_len: 12});
        let b = BusRecord {CB: 1, UMI: 1, EC: 1, COUNT: 3, FLAG: 0};
        let c = BusRecord {CB: 0, UMI: 1, EC: 1, COUNT: 3, FLAG: 0};
        let records = vec![b, c];
        wri.write_iterator(records.into_iter());

        let file = File::open(&tmpfilename).unwrap();
        let mut reader = BufReader::new(file);
        let mut buffer = Vec::new();
    
        // Read file into vector.
        reader.read_to_end(&mut buffer).unwrap();
        let digest_str = format!("{:x}", md5::compute(buffer));
        insta::assert_yaml_snapshot!(digest_str, @r###"
        ---
        3b49ecec15ae50a3cf2b7cfd37441ea1
        "###);
    }

    #[test]
    /// make sure the format on disk stays the same when we read it
    fn test_insta_binary_format_read() {
        let external_file = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
        let records: Vec<BusRecord> = BusReaderPlain::new(external_file).step_by(10000).take(100).collect();
        insta::assert_yaml_snapshot!(records)
    }

    #[test]
    fn test_dyn_from_read_file(){
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 1, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let rlist = vec![r1.clone(), r2.clone()]; // note: this clones r1, r2!
        let (busname, _dir) = setup_busfile(&rlist);

        let f = File::open(busname).unwrap();
        let mut r = BusReaderPlain::from_read(f);
        assert_eq!(r.params, BusParams {cb_len: 16, umi_len: 12});

        assert_eq!(r.next(), Some(r1));
        assert_eq!(r.next(), Some(r2));
        assert_eq!(r.next(), None);
    }

    #[test]
    fn test_dyn_from_file(){
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 1, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let rlist = vec![r1.clone(), r2.clone()]; // note: this clones r1, r2!
        let (busname, _dir) = setup_busfile(&rlist);

        let mut r = BusReaderPlain::new(&busname);
        assert_eq!(r.params, BusParams {cb_len: 16, umi_len: 12});

        assert_eq!(r.next(), Some(r1));
        assert_eq!(r.next(), Some(r2));
        assert_eq!(r.next(), None);
    }

    #[test]
    fn test_dyn_from_mem(){
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 1, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let rlist = vec![r1.clone(), r2.clone()]; // note: this clones r1, r2!
        let (busname, _dir) = setup_busfile(&rlist);

        let mut f = File::open(busname).unwrap();

        let mut buffer = Vec::new();
        f.read_to_end(&mut buffer).unwrap();


        let mut r = BusReaderPlain::from_read(buffer.as_slice());
        assert_eq!(r.params, BusParams {cb_len: 16, umi_len: 12});

        assert_eq!(r.next(), Some(r1));
        assert_eq!(r.next(), Some(r2));
        assert_eq!(r.next(), None);
    }

    #[test]
    fn test_dyn_group(){
        let r1 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 3, EC: 0, COUNT: 12, FLAG: 0 };
        let r3 = BusRecord { CB: 1, UMI: 21, EC: 1, COUNT: 2, FLAG: 0 };
        let rlist = vec![r1.clone(), r2.clone(), r3.clone()]; // note: this clones r1, r2!
        let (busname, _dir) = setup_busfile(&rlist);

        let mut f = File::open(busname).unwrap();

        let mut buffer = Vec::new();
        f.read_to_end(&mut buffer).unwrap();


        let r = BusReaderPlain::from_read(buffer.as_slice());

        let mut iter = r.groupby_cb();

        assert_eq!(iter.next(), Some((0, vec![r1, r2])));
        assert_eq!(iter.next(), Some((1, vec![r3])));
    }
}
