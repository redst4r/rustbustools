use std::collections::HashMap;
use std::fs::File;
use std::io::BufWriter;
use std::io::Read;
use std::io::BufReader;
use std::io::BufRead;
use std::io::Write;
use serde::{Serialize, Deserialize};
use std::io::{Seek, SeekFrom};
use bincode;
// use polars::prelude::*;
use regex::Regex;

const BUS_ENTRY_SIZE: usize = 32;
const BUS_HEADER_SIZE: usize = 20;

// BUS_ENTRY_SIZE = 32  # each entry is 32 bytes!!
//  unpack_str = 'QQiIIxxxx'
// Q: 8byte unsigned long,long int
// i: 4byte int
// I: unsigned int, 4byte
#[allow(non_snake_case)]
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone, Copy)]
pub struct BusRecord {
    pub CB: u64, //8byte
    pub UMI: u64, // 8byte
    pub EC: u32,  // 4v byte
    pub COUNT: u32, // 4v byte
    pub FLAG: u32, // 4v byte    including this, we have 28 bytes, missing 4 to fill up
    // PAD: u32 // just padding to fill 32bytes
}

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq)]
pub struct BusHeader {
//4sIIII: 20bytes
    magic: [u8; 4],
    version: u32,
    cb_len: u32,
    umi_len: u32,
    tlen: u32
}

impl BusHeader {
    pub fn new(cb_len: u32, umi_len: u32, tlen: u32) -> BusHeader{
        let magic: [u8; 4] = *b"BUS\x00";
        BusHeader {cb_len, umi_len, tlen, magic, version: 1}
    }

    pub fn from_file(fname: &String) -> BusHeader{
        // getting 20 bytes out of the file, which is the header
        let file = std::fs::File::open(fname);
        let mut header = Vec::with_capacity(BUS_HEADER_SIZE);
        let _n = file.as_ref().unwrap().take(BUS_HEADER_SIZE as u64).read_to_end(&mut header).unwrap();
        let header_struct: BusHeader = bincode::deserialize(&header).expect("FAILED to deserialze header");
        assert_eq!(&header_struct.magic, b"BUS\x00", "Header struct not matching");
        header_struct
    }
}


pub struct BusWriter{
    pub buf: BufWriter<File>,
    pub header: BusHeader
}
impl BusWriter {
    pub fn new(filename: &String, header: BusHeader) -> BusWriter{
        // let mut file_handle = std::fs::File::open(filename).expect("FAILED to open");
        let file_handle: File = File::create(filename).expect("FAILED to open");

        let mut buf = BufWriter::new(file_handle);

        // write the header into the file
        let binheader = bincode::serialize(&header).expect("FAILED to serialze header");
        buf.write(&binheader).expect("FAILED to write header");

        // write the variable header
        let mut varheader: Vec<u8> = Vec::new();
        for _i in 0..header.tlen{
            varheader.push(0);
        }
        // println!("len of varheader: {}" ,varheader.len());
        buf.write(&varheader).expect("FAILED to write var header");

        BusWriter {buf: buf , header: header}
    }

    pub fn write_record(&mut self, record: &BusRecord){
        let mut binrecord = bincode::serialize(record).expect("FAILED to serialze record");

        // the struct is only 28bytes, so we need 4 padding bytes
        for _i in 0..4{
            binrecord.push(0);
        }

        // println!("Length of record in bytes {}", binrecord.len());
        self.buf.write(&binrecord).expect("FAILED to write record");
    }
    pub fn write_records(&mut self, records: &Vec<BusRecord>){
        // writes several recordsd and flushes
        for r in records{
            self.write_record(r)
        }
        self.buf.flush().unwrap();

    }

}



// ===========Buffered Reader for bus======================
#[derive(Debug)]
pub struct BusIteratorBuffered {
    pub bus_header: BusHeader,
    buf: BufReader<File>
}

impl BusIteratorBuffered {
    pub fn new(filename: &String) -> BusIteratorBuffered{
        let bus_header = BusHeader::from_file(&filename);
        let mut file_handle = std::fs::File::open(filename).expect("FAIL");
        
        // advance the file point 20 bytes (pure header) + header.tlen (the variable part)
        let to_seek:u64 = BUS_HEADER_SIZE.try_into().unwrap();
        let hhh: u64 = bus_header.tlen.into();
        let _x = file_handle.seek(SeekFrom::Start(to_seek + hhh)).unwrap();

        let buf = BufReader::new(file_handle);
        // let mut buf = BufReader::with_capacity(8000, file_handle);
        BusIteratorBuffered {bus_header, buf }
    }
}

impl Iterator for BusIteratorBuffered {
    type Item = BusRecord;

    fn next(&mut self) -> Option<Self::Item> {
        let mut buffer = [0;BUS_ENTRY_SIZE];
        match self.buf.read(&mut buffer){
            Ok(BUS_ENTRY_SIZE) => Some(bincode::deserialize(&buffer).expect("deserial error")),
            Ok(0) => None,
            Ok(n) => {
                let s: BusRecord = bincode::deserialize(&buffer).expect("deserial error");
                panic!("{:?} {:?} {:?}", n, buffer, s)
            }
            Err(e) => panic!("{:?}", e),
        }
    }
}

//=================================================================================
pub struct CbUmiIterator {
    busiter: BusIteratorBuffered,
    last_record: Option<BusRecord>  //option needed to mark the final element of the iteration
}

impl CbUmiIterator {

    pub fn new(fname: &String) ->CbUmiIterator{
        let mut busiter = BusIteratorBuffered::new(fname);
        let last_record = busiter.next(); //initilize with the first record in the file
        CbUmiIterator {busiter, last_record}
    }
}

impl Iterator for CbUmiIterator {
    type Item = ((u64, u64), Vec<BusRecord>);

    fn next(&mut self) -> Option<Self::Item> {

        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted
        
        loop {
            if let Some(last_record) = self.last_record{  //if we're not done with the iteration
                // try to get a new record
                if let Some(new_record) = self.busiter.next(){

                    // determine if we encounter a new cell in this iteration
                    let current_cb = last_record.CB;
                    let current_umi = last_record.UMI;

                    if new_record.CB > current_cb || (new_record.CB == current_cb &&  new_record.UMI > current_umi){  
                        // we ran into a new CB/UMI and its records
                        busrecords.push(last_record); // the stored element from the previous iteration
                        // println!("\tyielding {:?}", (current_cb, &busrecords));
                        self.last_record = Some(new_record);

                        return Some(((current_cb, current_umi), busrecords));
                    }
                    else if (new_record.CB == current_cb) && new_record.UMI == current_umi {
                        busrecords.push(last_record);
                        self.last_record = Some(new_record);

                    }
                    else{  // the new cb is smaller then the current state: this is a bug due to an UNOSORTED busfile
                        panic!("Unsorted busfile: {}/{} -> {}/{}", current_cb, current_umi, new_record.CB, new_record.UMI)
                    }
                }
                else {
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    busrecords.push(last_record);
                    let current_cb = last_record.CB;
                    let current_umi = last_record.UMI;
                    // to mark the end of iteration and all items emitted, set last_item to None
                    self.last_record = None;
                    return Some(((current_cb, current_umi), busrecords));  
                }
            }
            else{  // last_record == None
                // we are done
                return None
            }
        }
    }
}

//=================================================================================
pub struct CellIterator {
    busiter: BusIteratorBuffered,
    last_record: Option<BusRecord>  //option needed to mark the final element of the iteration
}

impl CellIterator {

    pub fn new(fname: &String) ->CellIterator{
        let mut busiter = BusIteratorBuffered::new(fname);
        let last_record = busiter.next(); //initilize with the first record in the file
        CellIterator {busiter, last_record}
    }
}

impl Iterator for CellIterator {
    type Item = (u64, Vec<BusRecord>);

    fn next(&mut self) -> Option<Self::Item> {

        let mut busrecords: Vec<BusRecord> = Vec::new(); // storing the result to be emitted
        
        loop {
            if let Some(last_record) = self.last_record{  //if we're not done with the iteration
                // try to get a new record
                if let Some(new_record) = self.busiter.next(){

                    // determine if we encounter a new cell in this iteration
                    let current_cb = last_record.CB;
                    if new_record.CB > current_cb {  
                        // we ran into a new CB and its records
                        busrecords.push(last_record); // the stored element from the previous iteration
                        // println!("\tyielding {:?}", (current_cb, &busrecords));
                        self.last_record = Some(new_record);

                        return Some((current_cb, busrecords));
                    }
                    else if new_record.CB == current_cb {
                        busrecords.push(last_record);
                        self.last_record = Some(new_record);

                    }
                    else{  // the new cb is smaller then the current state: this is a bug due to an UNOSORTED busfile
                        panic!("Unsorted busfile: {} -> {}", current_cb, new_record.CB)
                    }
                }
                else {
                    // we ran pas the last entry of the file
                    // FINALize the last emit
                    busrecords.push(last_record);
                    let current_cb = last_record.CB;
                    // to mark the end of iteration and all items emitted, set last_item to None
                    self.last_record = None;
                    return Some((current_cb, busrecords));  
                }
            }
            else{  // last_record == None
                // we are done
                return None
            }
        }
    }
}

//=================================================================================

pub struct BusFolder {
    pub foldername: String,
    pub ec_dict: HashMap<u32, Vec<u32>>, // EC-> list of Transcript ids
    pub transcript_dict: HashMap<u32, String>,  // transcript-id -> traqnscript name
    pub transcript_to_gene: HashMap<String, String>,
    pub ec2gene: HashMap<u32, Vec<String>>,
    pub busfile: String
}

pub fn parse_ecmatrix(filename: String) -> HashMap<u32, Vec<u32>>{
    // parsing an ec.matrix into a Hashmap EC->list_of_transcript_ids
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut ec_dict: HashMap<u32, Vec<u32>> = HashMap::new();

    for line in reader.lines(){
        if let Ok(l) = line{
            // println!("{}", l);
            let mut s = l.split_whitespace();
            let ec = s.next().unwrap().parse::<u32>().unwrap();
            let tmp = s.next().unwrap();

            let transcript_list: Vec<u32> = tmp
                                             .split(",")
                                             .map(|x|x.parse::<u32>().unwrap()).collect();
            ec_dict.insert(ec, transcript_list);
        }
    }
    ec_dict
}


pub fn parse_ecmatrix2(filename: String) -> HashMap<u32, Vec<u32>>{

    let file = File::open(filename).unwrap();
    // parsing an ec.matrix into a Hashmap EC->list_of_transcript_ids
    let reader = BufReader::new(file);
    let mut ec_dict: HashMap<u32, Vec<u32>> = HashMap::new();

    let re = Regex::new(r"^(\d+)\s+(.+)$").unwrap();

    for line in reader.lines(){
        if let Ok(l) = line{
            let caps = re.captures(&l).unwrap();
            
            let ecid = caps.get(1).unwrap().as_str().parse::<u32>().unwrap();
            let t_string = caps.get(2).unwrap().as_str();
            let transcript_list: Vec<u32> = t_string
                .split(",")
                .map(|x|x.parse::<u32>().unwrap()).collect();
            ec_dict.insert(ecid, transcript_list);
        }
    }
    ec_dict
}

fn parse_transcript(filename: String) -> HashMap<u32, String>{
    let file = File::open(filename).unwrap();
    let reader = BufReader::new(file);
    let mut transcript_dict: HashMap<u32, String> = HashMap::new();
    for (i, line) in reader.lines().enumerate(){
        if let Ok(l) = line{
            transcript_dict.insert(i as u32, l);
        }
    }
    transcript_dict
}

fn parse_t2g(t2g_file: String) -> HashMap<String, String>{

    let mut t2g_dict: HashMap<String, String> = HashMap::new();        
    let file = File::open(t2g_file).unwrap();
    let reader = BufReader::new(file);
    for (_i, line) in reader.lines().enumerate(){
        if let Ok(l) = line{
            let mut s = l.split_whitespace();
            let transcript_id = s.next().unwrap();
            s.next();
            let symbol = s.next().unwrap();

            assert!(!t2g_dict.contains_key(&transcript_id.to_string()));  //make sure transcripts dont map to multiple genes
            t2g_dict.insert(transcript_id.to_string(), symbol.to_string());
        }
    }
    t2g_dict
}

impl BusFolder {
    pub fn new(foldername: String, t2g_file:String) ->BusFolder{

        // read EC dict
        println!("Reading EC.matrix");
        let filename= format!("{}/{}", foldername, "matrix.ec");
        let ec_dict = parse_ecmatrix(filename);

        // read transcript dict
        println!("Reading transcripts.txt");
        let filename = format!("{}/{}", foldername, "transcripts.txt");
        let transcript_dict = parse_transcript(filename);

        // read transcript to gene file
        println!("Reading t2g_dict");
        let t2g_dict = parse_t2g(t2g_file);

        // create the EC->gene mapping
        println!("building EC->gene");

        let mut ec2gene: HashMap<u32, Vec<String>> = HashMap::new();
        for (ec, transcript_ints) in ec_dict.iter(){

            let mut genes: Vec<String> = Vec::new();
            for t_int in transcript_ints{
                let t_name = transcript_dict.get(&t_int).unwrap();

                // if we can resolve, put the genename, otherwsise use the transcript name instead
                match t2g_dict.get(t_name){
                    Some(genename) => genes.push(genename.clone()),
                    None =>  genes.push(t_name.clone()),
                }
            }
            ec2gene.insert(*ec, genes);
        }


        // let df =  CsvReader::from_path(t2g_file).unwrap()
        //                                         .has_header(false)            
        //                                         .infer_schema(None);
        // let y = df.with_columns(Some(vec!["tid".to_string(),"gid".to_string(), "symbol".to_string()]));
        // v.finish();
        // let mut df = CsvReader::from_path(t2g_file).unwrap()
        // // .with_n_threads(Some(1)) // comment for multithreading
        // // .with_encoding(CsvEncoding::LossyUtf8)
        // .has_header(true)
        // .with_columns(Some(vec!["tid".to_string(),"gid".to_string(), "symbol".to_string()]))
        // .finish().unwrap();
        
        // let x = df["tid"].zip(df["symbol"]);
        let busfile = String::from("output.corrected.sort.bus");
        BusFolder{
            foldername,
            ec_dict,
            transcript_dict,
            transcript_to_gene: t2g_dict,
            ec2gene,
            busfile
        }
    }
}

pub fn group_record_by_cb_umi(record_list: Vec<BusRecord>) -> HashMap<(u64, u64), Vec<BusRecord>>{
    /*
    group a list of records by their CB/UMI (i.e. a single molecule)
    */
    let mut cbumi_map: HashMap<(u64, u64), Vec<BusRecord>> = HashMap::new();

    for r in record_list{
        let empty: Vec<BusRecord> = Vec::new();
        let rlist = cbumi_map.entry((r.CB, r.UMI)).or_insert(empty);
        rlist.push(r);
    }
    cbumi_map
}


//=================================================================================
#[cfg(test)]
mod tests {
    use std::io::Write;
    use crate::io::{BusRecord, BusHeader, CellIterator, BusIteratorBuffered, BusWriter, CbUmiIterator};

    fn setup_busfile(records: Vec<BusRecord>, busname: &String){
        let header = BusHeader::new(16, 12, 20);
        let mut writer = BusWriter::new(&busname, header);
        writer.write_records(&records);
    }

    #[test]
    fn test_read_write_header(){
        let r1 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let header = BusHeader::new(16, 12, 20);
        let busname = "/tmp/test_read_write_header.bus".to_string();
        let mut writer = BusWriter::new(&busname, header);
        writer.write_record(&r1);
        writer.buf.flush().unwrap();

        let bheader = BusHeader::from_file(&busname);
        let header = BusHeader::new(16, 12, 20);
        assert_eq!(header, bheader);
    }

    #[test]
    fn test_read_write(){
        let r1 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};

        let busname = "/tmp/test_read_write.bus".to_string();
        setup_busfile(vec![r1,r2], &busname);

        let mut reader = BusIteratorBuffered::new(&busname);
        let e1 = reader.next().unwrap();
        assert_eq!(e1, r1);
        // println!("{:?} {:?}", r1, e1);

        let e2 = reader.next().unwrap();
        assert_eq!(e2, r2);
        assert_eq!(reader.next(), None);

        // let records: Vec<BusRecord> = reader.into_iter().collect();
        // assert_eq!(records, vec![r1, r2])

    }

    #[test]
    #[should_panic(expected = "Unsorted busfile: 2 -> 0")]
    fn test_panic_on_unsorted(){  
        let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 2, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r4 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

        let records = vec![r1,r2,r3,r4];

        let busname = "/tmp/test_panic_on_unsorted.bus".to_string();
        setup_busfile(records, &busname);

        let cb_iter = CellIterator::new(&busname);
        let _n: Vec<(u64, Vec<BusRecord>)> = cb_iter.collect();
    }

    #[test]
    fn test_cb_iter(){   
        let r1 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r2 = BusRecord{CB: 0, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r3 = BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
        let r4 = BusRecord{CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
        let r5 = BusRecord{CB: 2, UMI: 21, EC: 1, COUNT: 2, FLAG: 0};
        let r6 = BusRecord{CB: 3, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};

        let records = vec![r1,r2,r3,r4,r5, r6];

        let busname = "/tmp/test_iter.bus".to_string();
        setup_busfile(records, &busname);


        let cb_iter = CellIterator::new(&busname);
        // let n: Vec<Vec<BusRecord>> = cb_iter.map(|(_cb, rec)| rec).collect();
        let n: Vec<(u64, Vec<BusRecord>)> = cb_iter.collect();
        // println!("{:?}", n);

        assert_eq!(n.len(), 4);
        // println!("{:?}", n);
        // println!("First");
        let c1 = &n[0];
        assert_eq!(*c1, (0, vec![r1,r2]));

        // println!("Second");
        let c2 = &n[1];
        assert_eq!(*c2, (1,vec![r3]));

        // println!("Third");
        let c3 = &n[2];
        assert_eq!(*c3, (2, vec![r4,r5]));

        let c4 = &n[3];
        assert_eq!(*c4, (3, vec![r6]));

        // assert_eq!(n, vec![vec![r1,r2], vec![r3], vec![r4,r5]])

     }


     #[test]
     fn test_cbumi_iter(){   
         let r1 = BusRecord{CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0};
         let r2 = BusRecord{CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
         let r3 = BusRecord{CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0};
         let r4 = BusRecord{CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
         let r5 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0};
         let r6 = BusRecord{CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0};
 
         let records = vec![r1,r2,r3,r4,r5, r6];
 
         let busname = "/tmp/test_cbumi_iter.bus".to_string();
         setup_busfile(records, &busname);
 
 
         let cb_iter = CbUmiIterator::new(&busname);
         // let n: Vec<Vec<BusRecord>> = cb_iter.map(|(_cb, rec)| rec).collect();
         let n: Vec<((u64, u64), Vec<BusRecord>)> = cb_iter.collect();
         println!("{:?}", n);
 
         assert_eq!(n.len(), 5);
         // println!("{:?}", n);
         // println!("First");
         let c1 = &n[0];
         assert_eq!(*c1, ((0,1 ), vec![r1,r2]));
 
         // println!("Second");
         let c2 = &n[1];
         assert_eq!(*c2, ((0, 2), vec![r3]));
 
         // println!("Third");
         let c3 = &n[2];
         assert_eq!(*c3, ((1,1), vec![r4]));
 
         let c4 = &n[3];
         assert_eq!(*c4, ((1,2), vec![r5]));
 
         let c5 = &n[4];
         assert_eq!(*c5, ((2,1), vec![r6]));
         // assert_eq!(n, vec![vec![r1,r2], vec![r3], vec![r4,r5]])
 
      }

    //   #[test]
    //   fn dummy(){
    //       let b= crate::BusFolder::new(
    //           // "/home/michi/mounts/TB4drive/ISB_data/LT_pilot/LT_pilot/kallisto_quant/Fresh1/kallisto/sort_bus/bus_output".to_string(),
    //           // "/home/michi/mounts/TB4drive/kallisto_resources/transcripts_to_genes.txt".to_string()
    //           "/home/michi/bus_output".to_string(),
    //           "/home/michi/transcripts_to_genes.txt".to_string()
    //       );
      
    //       println!("{}", b.ec_dict.len());
    //       println!("{:?}", b.ec_dict.get(&812069).unwrap());
    //       println!("{:?}", b.ec2gene.get(&613).unwrap());
          
    //   }
    use std::fs::File;
    use std::io::{BufWriter};

    // use super::parse_ecmatrix;
    use super::parse_ecmatrix2;
    #[test]
    fn test_ecmatrix(){
        let f = File::create("/tmp/foo").expect("Unable to create file");
        let mut f = BufWriter::new(f);

        let data = "0 1,2,3\n1 3,4\n2 4";
        f.write_all(data.as_bytes()).expect("Unable to write data");
        f.flush().unwrap();

        let ec= parse_ecmatrix2("/tmp/foo".to_string());
        // println!("{}", ec.len());
        println!("{:?}", ec);
        let e1 = ec.get(&0).unwrap();
        assert_eq!(*e1, vec![1,2,3]);

        let e1 = ec.get(&1).unwrap();
        assert_eq!(*e1, vec![3 ,4]);
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

        let records = vec![r1,r2,r3,r4,r5, r6];
        
        let grouped = group_record_by_cb_umi(records);


        let g1 = grouped.get(&(0 as u64,1 as u64)).unwrap();
        assert_eq!(*g1, vec![r1,r2]);

        let g2 = grouped.get(&(0 as u64,2 as u64)).unwrap();
        assert_eq!(*g2, vec![r3]);

        let g3 = grouped.get(&(1 as u64,1 as u64)).unwrap();
        assert_eq!(*g3, vec![r4, r6]);

        let g4 = grouped.get(&(1 as u64,2 as u64)).unwrap();
        assert_eq!(*g4, vec![r5]);

    }

}