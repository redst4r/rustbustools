use crate::consistent_genes::{Ec2GeneMapper, Genename, EC};
use bincode;
use serde::{Deserialize, Serialize};
use std::collections::{HashMap, HashSet};
use std::fs::File;
use std::io::{BufRead, BufReader, BufWriter, Read, Seek, SeekFrom, Write};
use tempfile::TempDir;

const BUS_ENTRY_SIZE: usize = 32;
const BUS_HEADER_SIZE: usize = 20;

// BUS_ENTRY_SIZE = 32  # each entry is 32 bytes!!
//  unpack_str = 'QQiIIxxxx'
// Q: 8byte unsigned long,long int
// i: 4byte int
// I: unsigned int, 4byte
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)] // TODO take away the copy and clone
pub struct BusRecord {
    pub CB: u64,    //8byte
    pub UMI: u64,   // 8byte
    pub EC: u32,    // 4v byte
    pub COUNT: u32, // 4v byte
    pub FLAG: u32,  // 4v byte    including this, we have 28 bytes, missing 4 to fill up
    // PAD: u32     // just padding to fill 32bytes
}
impl BusRecord {
    pub fn to_bytes(&self) -> Vec<u8> {
        let mut binrecord = bincode::serialize(self).expect("FAILED to serialze record");

        // the struct is only 28bytes, so we need 4 padding bytes
        for _i in 0..4 {
            //TODO feel like theres a better way to do this
            binrecord.push(0);
        }
        binrecord
    }
    pub fn from_bytes(bytes: &[u8]) -> Self {
        bincode::deserialize(bytes).expect("deserial error")
    }
}
/// Header of a busfile, as specified by kallisto
/// the only things that are variable are:
/// - cb_len: The number of bases in the cellbarcode (usually 16BP)
/// - umi_len: The number of bases in the cellbarcode (usually 12BP)
/// # Example
/// ```
/// let header = BusHeader::new(16, 12, 1)
/// // can also be obtained from an existing busfile
/// let header = BusHeader::from_file("somefile.bus")
/// ```
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq)]
pub struct BusHeader {
    //4sIIII: 20bytes
    magic: [u8; 4],
    version: u32,
    pub cb_len: u32,
    umi_len: u32,
    tlen: u32,
}

impl BusHeader {
    pub fn new(cb_len: u32, umi_len: u32, tlen: u32) -> BusHeader {
        let magic: [u8; 4] = *b"BUS\x00";
        BusHeader {cb_len, umi_len, tlen, magic, version: 1}
    }

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

    pub fn from_bytes(bytes: &[u8]) -> BusHeader {
        let header_struct: BusHeader =
            bincode::deserialize(bytes).expect("FAILED to deserialze header");
        assert_eq!(
            &header_struct.magic, b"BUS\x00",
            "Header struct not matching; MAGIC is wrong"
        );
        header_struct
    }

    pub fn to_bytes(&self) -> Vec<u8> {
        bincode::serialize(self).expect("FAILED to serialze header")
    }

    pub fn get_tlen(&self) -> u32 {
        self.tlen
    }
}

/// A marker trait for iterators over CB/UMI/gene_EC iterators
/// to create an iterator compatible with this framework,
/// just implement the usual iterator trait (giving the next() function)
/// and also tag the iterator with CUGIterator: `impl CUGIterator for Dummy {}`
/// ## Example
/// ```
/// struct Dummy{
///     f: i64
/// }
///
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

pub struct BusReader {
    pub bus_header: BusHeader,
    reader: BufReader<File>,
}

// benchmarking had a sligh incrase of speed using 800KB instead of *Kb
// further increase buffers dont speed things up more (just need more mem)
// const DEFAULT_BUF_SIZE: usize = 8 * 1024;  //8KB
const DEFAULT_BUF_SIZE: usize = 800 * 1024; // 800  KB
impl BusReader {
    pub fn new(filename: &str) -> BusReader {
        BusReader::new_with_capacity(filename, DEFAULT_BUF_SIZE)
    }
    pub fn new_with_capacity(filename: &str, bufsize: usize) -> BusReader {
        let bus_header = BusHeader::from_file(filename);
        let mut file_handle = File::open(filename).expect("FAIL");

        // advance the file point 20 bytes (pure header) + header.tlen (the variable part)
        let to_seek = BUS_HEADER_SIZE as u64 + bus_header.tlen as u64;
        let _x = file_handle.seek(SeekFrom::Start(to_seek)).unwrap();

        let buf = BufReader::with_capacity(bufsize, file_handle);
        BusReader {bus_header, reader: buf}
    }
}

impl Iterator for BusReader {
    type Item = BusRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buffer = [0; BUS_ENTRY_SIZE];
        let r = match self.reader.read(&mut buffer) {
            Ok(BUS_ENTRY_SIZE) => Some(BusRecord::from_bytes(&buffer)),
            Ok(0) => None,
            Ok(n) => {
                let s: BusRecord = BusRecord::from_bytes(&buffer);
                panic!(
                    "Wrong number of bytes {:?}. Buffer: {:?} record: {:?}",
                    n, buffer, s
                )
            }
            Err(e) => panic!("{:?}", e),
        };
        r
    }
}
// tag our iterator to be compatible with our framework
impl CUGIterator for BusReader {}

// ========================================
pub struct BusWriter {
    pub writer: BufWriter<File>,
    pub header: BusHeader,
}

impl BusWriter {
    pub fn from_filehandle(file_handle: File, header: BusHeader) -> BusWriter {
        BusWriter::new_with_capacity(file_handle, header, DEFAULT_BUF_SIZE)
    }

    pub fn new(filename: &str, header: BusHeader) -> BusWriter {
        let file_handle: File = File::create(filename).expect("FAILED to open");
        BusWriter::from_filehandle(file_handle, header)
    }

    pub fn write_record(&mut self, record: &BusRecord) {
        let binrecord = record.to_bytes();
        self.writer
            .write_all(&binrecord)
            .expect("FAILED to write record");
    }
    pub fn write_records(&mut self, records: &Vec<BusRecord>) {
        // writes several recordsd and flushes
        for r in records {
            self.write_record(r)
        }
        self.writer.flush().unwrap();
    }

    pub fn new_with_capacity(file_handle: File, header: BusHeader, bufsize: usize) -> Self {
        let mut writer = BufWriter::with_capacity(bufsize, file_handle);

        // write the header into the file
        let binheader =header.to_bytes();
        writer.write_all(&binheader).expect("FAILED to write header");

        // write the variable header
        let mut varheader: Vec<u8> = Vec::new();
        for _i in 0..header.tlen {
            varheader.push(0);
        }
        writer.write_all(&varheader).expect("FAILED to write var header");

        BusWriter { writer, header }
    }
}
//=================================================================================

pub struct BusFolder {
    pub foldername: String,
    pub ec2gene: Ec2GeneMapper,
}

pub fn parse_ecmatrix(filename: &str) -> HashMap<EC, Vec<u32>> {
    // parsing an ec.matrix into a Hashmap EC->list_of_transcript_ids
    let file = File::open(filename).unwrap_or_else(|_| panic!("{} not found", filename));
    let reader = BufReader::new(file);
    let mut ec_dict: HashMap<EC, Vec<u32>> = HashMap::new();

    for line in reader.lines() {
        if let Ok(l) = line {
            // println!("{}", l);
            let mut s = l.split_whitespace();
            let ec = s.next().unwrap().parse::<u32>().unwrap();
            let tmp = s.next().unwrap();

            let transcript_list: Vec<u32> =
                tmp.split(',').map(|x| x.parse::<u32>().unwrap()).collect();
            ec_dict.insert(EC(ec), transcript_list);
        } else {
            panic!("Error reading file {}", filename)
        }
    }
    ec_dict
}

fn parse_transcript(filename: &str) -> HashMap<u32, String> {
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut transcript_dict: HashMap<u32, String> = HashMap::new();
    for (i, line) in reader.lines().enumerate() {
        if let Ok(l) = line {
            transcript_dict.insert(i as u32, l);
        }
    }
    transcript_dict
}

fn parse_t2g(t2g_file: &str) -> HashMap<String, Genename> {
    let mut t2g_dict: HashMap<String, Genename> = HashMap::new();
    let file = File::open(t2g_file).unwrap_or_else(|_| panic!("{} not found", t2g_file));
    let reader = BufReader::new(file);
    for (_i, line) in reader.lines().enumerate() {
        if let Ok(l) = line {
            let mut s = l.split_whitespace();
            let transcript_id = s.next().unwrap();
            let ensemble_id = s.next().unwrap();
            let _symbol = s.next().unwrap();

            assert!(!t2g_dict.contains_key(&transcript_id.to_string())); //make sure transcripts dont map to multiple genes
            t2g_dict.insert(transcript_id.to_string(), Genename(ensemble_id.to_string()));
        }
    }
    t2g_dict
}

fn build_ec2gene(
    ec_dict: &HashMap<EC, Vec<u32>>,
    transcript_dict: &HashMap<u32, String>,
    t2g_dict: &HashMap<String, Genename>,
) -> HashMap<EC, HashSet<Genename>> {
    let mut ec2gene: HashMap<EC, HashSet<Genename>> = HashMap::new();

    for (ec, transcript_ints) in ec_dict.iter() {
        let mut genes: HashSet<Genename> = HashSet::new();

        for t_int in transcript_ints {
            let t_name = transcript_dict.get(t_int).unwrap();

            // if we can resolve, put the genename, otherwsise use the transcript name instead
            //
            // actually, turns out that kallisto/bustools treats it differently:
            // if the transcript doenst resolve, just drop tha trasncrip from the EC set
            // TODO what happens if non of the EC transcripts resolve
            match t2g_dict.get(t_name) {
                Some(genename) => {
                    genes.insert(genename.clone());
                }
                // None =>  genes.insert(Genename(t_name.clone())),
                None => {}
            };
        }
        ec2gene.insert(*ec, genes);
    }
    // // sanity check, make sure no Ec set is empty
    // for (ec, gset) in ec2gene.iter(){
    //     if gset.is_empty(){
    //         println!("{ec:?}'s geneset is empty");
    //     }
    // }

    ec2gene
}

impl BusFolder {
    pub fn new(foldername: &str, t2g_file: &str) -> BusFolder {
        println!("building EC->gene for {}", foldername);
        // read EC dict
        let filename = format!("{}/{}", foldername, "matrix.ec");
        let ec_dict = parse_ecmatrix(&filename);

        // read transcript dict
        let filename = format!("{}/{}", foldername, "transcripts.txt");
        let transcript_dict = parse_transcript(&filename);

        // read transcript to gene file
        let t2g_dict = parse_t2g(t2g_file);

        // create the EC->gene mapping
        let ec2gene = build_ec2gene(&ec_dict, &transcript_dict, &t2g_dict);
        let ecmapper = Ec2GeneMapper::new(ec2gene);

        BusFolder {
            foldername: foldername.to_string(),
            ec2gene: ecmapper,
        }
    }

    pub fn get_iterator(&self) -> BusReader {
        let bfile = self.get_busfile();
        BusReader::new(&bfile)
    }

    pub fn get_busfile(&self) -> String {
        format!("{}/output.corrected.sort.bus", self.foldername)
    }

    pub fn get_bus_header(&self) -> BusHeader {
        BusHeader::from_file(&self.get_busfile())
    }

    pub fn get_ecmatrix_file(&self) -> String {
        format!("{}/matrix.ec", self.foldername)
    }

    pub fn get_transcript_file(&self) -> String {
        format!("{}/transcripts.txt", self.foldername)
    }

    pub fn parse_ecmatrix(&self) -> HashMap<EC, Vec<u32>> {
        let filename = self.get_ecmatrix_file();
        parse_ecmatrix(&filename)
    }

    pub fn parse_transcript(&self) -> HashMap<u32, String> {
        let filename = self.get_transcript_file();
        parse_transcript(&filename)
    }
    // pub fn build_ec2gene(&self, ) -> HashMap<u32, HashSet<String>>{
    //     let ec_dict = self.parse_ecmatrix();
    //     let transcript_dict = self.parse_transcript();
    //     let t2g_dict = parse_t2g(&self.t2g_file);
    //     build_ec2gene(&ec_dict, &transcript_dict, &t2g_dict)
    // }
}

pub fn group_record_by_cb_umi(record_list: Vec<BusRecord>) -> HashMap<(u64, u64), Vec<BusRecord>> {
    /*
    group a list of records by their CB/UMI (i.e. a single molecule)

    takes ownership of record_list, since it reorders the elements
    */
    let mut cbumi_map: HashMap<(u64, u64), Vec<BusRecord>> = HashMap::new();

    for r in record_list {
        let rlist = cbumi_map.entry((r.CB, r.UMI)).or_insert(Vec::new());
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

    let header = BusHeader::new(16, 12, 20);
    let mut writer = BusWriter::new(tmpfilename, header);
    writer.write_records(records);

    (tmpfilename.to_string(), dir)
}

pub fn write_partial_busfile(bfile: &str, boutfile: &str, nrecords: usize) {
    // write the first nrecords of the intput file into the output
    let busiter = BusReader::new(bfile);
    let newheader = BusHeader::new(
        busiter.bus_header.cb_len,
        busiter.bus_header.umi_len,
        busiter.bus_header.tlen,
    );
    let mut buswriter = BusWriter::new(boutfile, newheader);

    for record in busiter.take(nrecords) {
        buswriter.write_record(&record);
    }
}

#[test]
fn test_write() {
    // let fname = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
    let fname = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
    // let outname = "/home/michi/bus_testing/bus_output_short/output.corrected.sort.bus";
    // write_partial_busfile(fname, outname, 10_000_000);

    let outname = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";
    write_partial_busfile(fname, outname, 1_500_000);
}

// #[test]
// fn testing2(){
//     use itertools::Itertools;
//     let FOLDERNAME:&str = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
//     let n =100000;
//     let biter= BusIteratorBuffered::new(&FOLDERNAME);
//     let s:Vec<Vec<BusRecord>> = biter
//         .group_by(|r| r.CB)
//         .into_iter().map(|(a, records)|records.collect())
//         .take(n)
//         .collect();
// }

//=================================================================================
#[cfg(test)]
mod tests {
    use crate::consistent_genes::EC;
    use crate::io::{setup_busfile, BusHeader, BusReader, BusRecord, BusWriter};
    use std::io::Write;

    #[test]
    fn test_read_write_header(){
        let r1 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let header = BusHeader::new(16, 12, 20);
        let busname = "/tmp/test_read_write_header.bus";
        let mut writer = BusWriter::new(busname, header);
        writer.write_record(&r1);
        writer.writer.flush().unwrap();

        let bheader = BusHeader::from_file(&busname);
        let header = BusHeader::new(16, 12, 20);
        assert_eq!(header, bheader);
    }

    #[test]
    fn test_read_write(){
        let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 1, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let rlist = vec![r1,r2];   // note: this clones r1, r2!
        let (busname, _dir) = setup_busfile(&rlist);
        let reader = BusReader::new(&busname);

        let records: Vec<BusRecord> = reader.into_iter().collect();
        assert_eq!(records, rlist)
    }

    use std::fs::File;
    use std::io::BufWriter;

    use super::parse_ecmatrix;
    #[test]
    fn test_ecmatrix() {
        let f = File::create("/tmp/foo").expect("Unable to create file");
        let mut f = BufWriter::new(f);

        let data = "0 1,2,3\n1 3,4\n2 4";
        f.write_all(data.as_bytes()).expect("Unable to write data");
        f.flush().unwrap();

        let ec = parse_ecmatrix("/tmp/foo");
        // println!("{}", ec.len());
        // println!("{:?}", ec);
        let e1 = ec.get(&EC(0)).unwrap();
        assert_eq!(*e1, vec![1, 2, 3]);

        let e1 = ec.get(&EC(1)).unwrap();
        assert_eq!(*e1, vec![3, 4]);
    }

    use super::group_record_by_cb_umi;
    #[test]
    fn test_grouping(){
        let r1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r4 = BusRecord{CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r5 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0};
        let r6 = BusRecord{CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

        let records = vec![r1.clone(),r2.clone(),r3.clone(),r4.clone(),r5.clone(), r6.clone()];
        
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
}
