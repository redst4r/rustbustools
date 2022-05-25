use std::fs::File;
use std::io::Read;
use std::io::BufReader;
use serde::{Serialize, Deserialize};
use bincode;


const BUS_ENTRY_SIZE: usize = 32;
const BUS_HEADER_SIZE: usize = 20;

// BUS_ENTRY_SIZE = 32  # each entry is 32 bytes!!
//  unpack_str = 'QQiIIxxxx'
// Q: 8byte unsigned long,long int
// i: 4byte int
// I: unsigned int, 4byte
#[derive(Serialize, Deserialize, Debug)]
pub struct BusRecord {
    pub CB: u64, //8byte
    pub UMI: u64, // 8byte
    pub EC: i32,  // 4v byte
    pub COUNT: u32, // 4v byte
    pub FLAG: u32, // 4v byte    including this, we have 28 bytes, missing 4 to fill up
    // PAD: u32 // just padding to fill 32bytes
}

#[derive(Serialize, Deserialize, Debug)]
pub struct BusHeader {
//4sIIII: 20bytes
    magic: [u8; 4],
    version: u32,
    cb_len: u32,
    umi_len: u32,
    tlen: u32
}


pub fn read_bus_header(fname: &str) -> BusHeader {
    // getting 20 bytes out of the file, which is the header
    let file = std::fs::File::open(fname);
    let mut header = Vec::with_capacity(BUS_HEADER_SIZE);
    let _n = file.as_ref().unwrap().take(BUS_HEADER_SIZE as u64).read_to_end(&mut header).unwrap();
    let header_struct: BusHeader = bincode::deserialize(&header).unwrap();

    assert_eq!(&header_struct.magic, b"BUS\x00");
    header_struct
}

pub fn read_bus_records(fname: &str) -> Vec<BusRecord>{

    let mut file = std::fs::File::open(fname).expect("FAIL");

    let header = read_bus_header(fname);
    // move the cursor across the header to the first entry
    let to_seek:u64 = BUS_HEADER_SIZE.try_into().unwrap();
    let hhh: u64 = header.tlen.into();
    let _x = file.seek(SeekFrom::Start(to_seek + hhh)).unwrap();


    let mut buf = BufReader::new(file);
    let mut records: Vec<BusRecord> = Vec::new();
    let mut record_bytes = [0; BUS_ENTRY_SIZE];
    loop {
        match buf.read(&mut record_bytes) {
            Ok(0) => break,
            Ok(BUS_ENTRY_SIZE) => records.push(bincode::deserialize(&record_bytes).unwrap()),
            Ok(n) => panic!("{:?}", n),
            Err(e) => panic!("{:?}", e),
        };
    }
    records
}

#[derive(Debug)]
pub struct BusIteratorBuffered {
    filename: String,
    pub bus_header: BusHeader,
    buf: BufReader<File>
}
use std::io::{Seek, SeekFrom};

impl  BusIteratorBuffered {
    pub fn new(filename: &String) -> BusIteratorBuffered{
        let bus_header = read_bus_header(filename.as_str());
        let mut file_handle = std::fs::File::open(filename).expect("FAIL");
        
        // advance the file point 20 bytes (pure header) + header.tlen (the variable part)
        let to_seek:u64 = BUS_HEADER_SIZE.try_into().unwrap();
        let hhh: u64 = bus_header.tlen.into();
        let _x = file_handle.seek(SeekFrom::Start(to_seek + hhh)).unwrap();

        let mut buf = BufReader::new(file_handle);
        // let mut buf = BufReader::with_capacity(8000, file_handle);
        BusIteratorBuffered {filename: filename.to_string(), bus_header, buf }
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
