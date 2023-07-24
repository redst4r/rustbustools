//! Dealing with the [busz compression format](https://github.com/BUStools/BUSZ-format)
//! 
//! ## About Bitvec
//! This code relies heavily on BitVec. It uses BitVec to encode/decode
//! the bits of the busz records, in particular Fibbonnaci encoding and NewPFD encoding
//! **A certain peculiarity though**:
//! To turn bytes (e.g from a u64 or read from the file) into Bitvec we use `BitVec::from_bytes(byte_Array)`
//! This takes the bytes literally in the order of the array.
//! Yet `bustools` writes busz in little endian format, i.e. the byte order is reversed.
//! In particular, each busz block contains entries for CB,UMI... each PADDED with zeros afterwards(to a multiple of 64)
//! On disk this is how it looks like:
//! ```
//! 0000000...00000000[CBs in Fibbonnaci]
//! 0000000...00000000[UMIs in Fibbonnaci]
//! ```
//! 
//! Even more, the fibbonacci encoding must be done with little endian byte order, if on disk it looks like
//! ```bash,no_run
//! aaaaaaaabbbbbbbbccccccccddddddddeeeeeeeeffffffffgggggggghhhhhhhh  //bits
//! ```
//! the correct fibonacci stream to decode is
//! ```
//! ddddddddccccccccbbbbbbbbaaaaaaaahhhhhhhhgggggggg....
//! ``` 
use std::{io::{BufWriter, Write, BufReader, Read, SeekFrom, Seek}, fs::File, collections::VecDeque};

use bitvec::{prelude as bv, field::BitField};
// use fibonacci_codec::Encode;
use itertools::{Itertools, izip};
use serde::{Serialize, Deserialize};
use crate::io::{BusRecord, BusReader, BusHeader, BUS_HEADER_SIZE, DEFAULT_BUF_SIZE, BusWriter};
use crate::runlength_codec::RunlengthCodec;


#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
pub struct BuszSpecificHeader {
    block_size: u32,
    pfd_block_size: u32,
    lossy_umi: u32
}
impl BuszSpecificHeader {
    fn to_bytes(&self) -> Vec<u8>{
        bincode::serialize(self).expect("FAILED to serialze header")
    }
}

const PFD_BLOCKSIZE: usize = 512;

pub struct CompressedBlock {

}

pub struct CompressedBlockHeader {
    // the 34 most significant bits denote the size of the compressed block in bytes. 
    // The 30 least significant bits denote the number of BUS records in the block.
    header_bytes: u64  
}

pub fn setbits_u32(x: u8) -> u32 {
    u32::MAX >> (32 - x)
}

pub fn setbits_u64(x: u8) -> u64 {
    u64::MAX >> (64 - x)
}

impl CompressedBlockHeader {
    pub fn new(block_size_bytes: u64, block_size_records:u64) -> Self {
        let header_bytes = (block_size_bytes << 30) | block_size_records ;

        // we only have 30 bits to store the nrecords
        if (setbits_u32(30) as u64) <= block_size_records {
            panic!("Cant store more than {} records, trying {}", setbits_u32(30), block_size_records)
        }
        // we only have 34 bits to store the blocksize
        if setbits_u64(34) <= block_size_bytes {
            panic!("Cant store more than {} records, trying {}", setbits_u32(34), block_size_bytes)
        }        
        CompressedBlockHeader { header_bytes }
    }

    /// decodes the header bytes into blocksize and number of records
    /// folllowing the header
    pub fn get_blocksize_and_nrecords(&self) -> (u64, u64) {
        let bit_length = 30; // encoding scheme imposed by bustools

        // dbg!(display_u64_in_bits(self.header_bytes));

        let block_size_bytes = self.header_bytes >> bit_length;
        // dbg!(display_u64_in_bits(block_size_bytes));

        // let bitmask: u32 = (!0) << (32 - bit_length);
        let bitmask = setbits_u32(bit_length);
        // dbg!(display_u32_in_bits(bitmask));
        
        let bitmask_64 = bitmask as u64;
        let block_size_records = self.header_bytes & bitmask_64;
        // dbg!(display_u64_in_bits(block_size_records));
        (block_size_bytes, block_size_records)
    }
}

pub fn display_u64_in_bits(x: u64) -> String{
    let s: Vec<u64> = (0..64).rev().map (|n| (x >> n) & 1).collect();
    s.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join("")
}

pub fn display_u32_in_bits(x: u32) -> String{
    let s: Vec<u32> = (0..32).rev().map (|n| (x >> n) & 1).collect();
    s.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join("")
}

/// round an integer to the next bigger multiple
/// ```rust
///  use bustools::busz::round_to_multiple;
///  assert_eq!(round_to_multiple(10,10), 10);
///  assert_eq!(round_to_multiple(11,10), 20);
///  assert_eq!(round_to_multiple(6,5), 10);
/// ```
pub fn round_to_multiple(i: usize, multiple: usize) -> usize {
    ((i+multiple-1)/multiple)*multiple
}

fn compress_barcodes2(records: &[BusRecord]) -> bv::BitVec<u8,bv::Msb0> {
    let runlength_codec = RunlengthCodec {RLE_VAL: 0, shift_up_1: true};

    let mut cb_iter = records.iter().map(|r| r.CB);
    let mut last_el = cb_iter.next().unwrap();
    let mut delta_encoded = Vec::new();
    delta_encoded.push(last_el);
    for el in cb_iter{
        delta_encoded.push(el-last_el);
        last_el=el
    }
    // println!("Delta enc: {:?}",delta_encoded);
    let runlen_encoded = runlength_codec.encode(delta_encoded.into_iter());
    // println!("run enc: {:?}",runlen_encoded);

    //fibbonaci encoding
    // let mut enc =runlen_encoded.fib_encode().unwrap();
    let mut enc: bv::BitVec<u8, bv::Msb0> = bv::BitVec::new();
    for fib_encoded in runlen_encoded.into_iter().map(newpfd::fibonacci::fib_enc) {
        enc.extend(fib_encoded);
    }


    // pad to next multiple of 64
    let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    // enc.grow(n_pad, false);
    for _ in 0..n_pad {
        enc.push(false);
    }
    
    // pad to next multiple of 64
    // let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    // enc.grow(n_pad, false);
    enc
}

// TODO: the first value encoded is a little funny, it gets incremented 2x!!
// periodic runlength encoding since the UMI resets to smaller values (when changing CB)
// and we cant handle negative differences!
fn compress_umis(records: &[BusRecord]) -> bv::BitVec<u8, bv::Msb0> {
    let runlength_codec = RunlengthCodec {RLE_VAL: 0,shift_up_1: true};
    let mut periodic_delta_encoded = Vec::new();
    let iii = records.iter();
    let last_record = &records[0];
    let mut last_umi = 0;
    let mut last_bc = last_record.CB + 1;
    let mut bc: u64;
    let mut umi: u64;
    for current_record in iii {
        bc = current_record.CB;
        umi = current_record.UMI + 1;
        if bc != last_bc {
            last_umi= 0;
        }
        let diff = umi - last_umi;

        periodic_delta_encoded.push(diff);
        last_umi = umi;
        last_bc = bc;

    };
    let runlen_encoded = runlength_codec.encode(periodic_delta_encoded.into_iter());
    //fibbonaci encoding
    let mut enc: bv::BitVec<u8, bv::Msb0> = bv::BitVec::with_capacity(2*runlen_encoded.len()); // the capacity is a minimum, assuming each RLE is a 1, i.e. `11` in fib encoding
    for mut fib_encoded in runlen_encoded.into_iter().map(newpfd::fibonacci::fib_enc) {
        enc.append(&mut fib_encoded);
    }

    // pad to next multiple of 64
    let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    // enc.grow(n_pad, false);
    for _ in 0..n_pad {
        enc.push(false);
    }
    enc
}


/// Compress ECs with NewPFD encoding
fn compress_ecs(records: &[BusRecord]) -> bv::BitVec<u8, bv::Msb0> {
    let ecs = records.iter().map(|r|r.EC as u64);
    let (mut encoded, _n_el) = newpfd::newpfd_bitvec::encode(ecs, PFD_BLOCKSIZE);

    // let n_pad = round_to_multiple(encoded.len(), 64) - encoded.len();
    let n_pad = round_to_multiple(encoded.len(), 32) - encoded.len();
    for _ in 0..n_pad {
        encoded.push(false);
    }
    assert_eq!(encoded.len() % 32, 0,  "newPFD block size needs to be a mutiple of 64, but is {}", encoded.len());
    
    // pad to next multiple of 64
    // encoded.grow(n_pad, false);

    // the rather strange swapping around of endianess, see parse_ec() too
    let bytes = bitslice_to_bytes(encoded.as_bitslice()); 
    let a = swap_endian(&bytes, 8);
    let swapped = swap_endian(&a, 4);
    let swapped_bv: bv::BitVec<u8, bv::Msb0> =  bv::BitVec::from_slice(&swapped);

    swapped_bv
}

/// Compress counts with RunLength(1) encoding
fn compress_counts(records: &[BusRecord]) -> bv::BitVec<u8, bv::Msb0> {
    let runlength_codec = RunlengthCodec {RLE_VAL: 1, shift_up_1: false};
    let count_iter = records.iter().map(|r| r.COUNT as u64);

    let runlen_encoded = runlength_codec.encode(count_iter);
    // println!("Counts: run enc: {:?}",runlen_encoded);

    //fibbonaci encoding
    // if true {
    let mut enc: bv::BitVec<u8, bv::Msb0> = bv::BitVec::new();
    for fib_encoded in runlen_encoded.into_iter().map(newpfd::fibonacci::fib_enc) {
        enc.extend(fib_encoded);
    }

    // } else {
        // let mut enc = runlen_encoded.fib_encode().unwrap();

    // // pad to next multiple of 64
    let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    for _ in 0..n_pad {
        enc.push(false);
    }
    enc
}

/// Compress flags with RunLength(0) encoding
fn compress_flags(records: &[BusRecord]) -> bv::BitVec<u8, bv::Msb0> {
    let runlength_codec = RunlengthCodec {RLE_VAL: 0, shift_up_1: true};
    let flag_iter = records.iter().map(|r| r.FLAG as u64);

    let runlen_encoded = runlength_codec.encode(flag_iter);
    // println!("run enc: {:?}",runlen_encoded);

    //fibbonaci encoding
    // let mut enc = runlen_encoded.fib_encode().unwrap();
    let mut enc: bv::BitVec<u8, bv::Msb0> = bv::BitVec::new();
    for fib_encoded in runlen_encoded.into_iter().map(newpfd::fibonacci::fib_enc) {
        enc.extend(fib_encoded);
    }

    // pad to next multiple of 64
    let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    // enc.grow(n_pad, false);
    for _ in 0..n_pad {
        enc.push(false);
    }

    enc
}

/// puts all given records into a single busz-block, including header
fn compress_busrecords_into_block(records: &[BusRecord]) -> Vec<u8> {//bv::BitVec<u8, bv::Msb0> {
    let bcs = compress_barcodes2(records);
    let umis = compress_umis(records);
    let ecs = compress_ecs(records);
    let counts = compress_counts(records);
    let flags = compress_flags(records);

    // bytes
    let bitsize: usize = bcs.len()+umis.len()+ecs.len()+counts.len()+flags.len();
    assert_eq!(bitsize%8, 0);
    let mut body_bytes = Vec::with_capacity(bitsize/8);

    body_bytes.extend(swap_endian(&bitslice_to_bytes(&bcs), 8));
    body_bytes.extend(swap_endian(&bitslice_to_bytes(&umis), 8));
    body_bytes.extend(swap_endian(&bitslice_to_bytes(&ecs), 8));
    body_bytes.extend(swap_endian(&bitslice_to_bytes(&counts), 8));
    body_bytes.extend(swap_endian(&bitslice_to_bytes(&flags), 8));

    let nbytes = body_bytes.len();
    let header = CompressedBlockHeader::new(
        nbytes.try_into().unwrap(), 
        records.len().try_into().unwrap());
    
    let mut header_bytes = header.header_bytes.to_le_bytes().to_vec();
    header_bytes.extend(body_bytes); 

    header_bytes
    // let _a = bcs.as_bitslice().load_be();
    // assert_eq!(
        // _a,
        // bitslice_to_bytes(&bcs));
  

    // println!("encode CB pos:{}", 0);
    // println!("encode UMI pos:{}", bcs.len());
    // //concat
    // bcs.append(&mut umis);
    // println!("encode ecs pos:{}", bcs.len());

    // bcs.append(&mut ecs);
    // println!("encode counts pos:{}", bcs.len());

    // bcs.append(&mut counts);
    // println!("encode flag pos:{}", bcs.len());
    // bcs.append(&mut flags);
    // println!("encode end pos:{}", bcs.len());

    // // need to be a multiple of 8
    // assert_eq!(bcs.len() % 8 , 0);

    // // need to be a multiple of 864
    // assert_eq!(bcs.len() % 64 , 0);

    // let nbytes = bcs.len() / 8;
    // let header = CompressedBlockHeader::new(
    //     nbytes.try_into().unwrap(), 
    //     records.len().try_into().unwrap());

    // // let mut header_bits = bit_vec::BitVec::from_bytes(&header.header_bytes.to_be_bytes()) ; 
    // let mut header_bits: bv::BitVec<u8, bv::Msb0> = bv::BitVec::from_slice(&header.header_bytes.to_be_bytes()) ; // definitely BE here

    // header_bits.append(&mut bcs);
    // header_bits
    // // header.header_bytes
}

pub fn compress_busfile(input: &str, output: &str, blocksize: usize) {

    let reader = BusReader::new(input);
    let out_fh: File = File::create(output).expect("FAILED to create output file");
    
    // write BusZ header
    let mut writer = BufWriter::with_capacity(DEFAULT_BUF_SIZE, out_fh);
    // write the header into the file
    let magic: [u8; 4] = *b"BUS\x01";
    let busheader = BusHeader{
        magic, 
        version: reader.bus_header.version, 
        cb_len: reader.bus_header.cb_len, 
        umi_len: reader.bus_header.umi_len, 
        tlen: 0 //reader.bus_header.tlen
    };
    let binheader = busheader.to_bytes();
    writer
        .write_all(&binheader)
        .expect("FAILED to write header");

    // write the variable header
    let mut varheader: Vec<u8> = Vec::new();
    for _i in 0..busheader.tlen {
        varheader.push(0);
    }
    writer
        .write_all(&varheader)
        .expect("FAILED to write var header");

    // BusZ header
    let busz_header = BuszSpecificHeader {
        block_size: blocksize.try_into().unwrap(),
        pfd_block_size: 512,
        lossy_umi: 0
    };
    let binzheader = busz_header.to_bytes();
    writer
        .write_all(&binzheader)
        .expect("FAILED to write header");

    for chunk in &reader.chunks(blocksize) {
        let records: Vec<BusRecord> = chunk.collect();
        let compressed_block = compress_busrecords_into_block(&records);
        let little_endian = compressed_block;
        // BitVec.to_bytes() implicitly does bytes in big-endian (first bit in the stream is the high-order bit of the low order byte),
        // yet in the busz file its stored in little endian and we need to convert this
        // let little_endian = swap_endian(&compressed_block.to_bytes(), 8);
        // let little_endian = swap_endian(&bitslice_to_bytes(&compressed_block), 8);

        writer.write_all(&little_endian).unwrap();
    }

    // write a final, emtpy block header (0_u64) that signals the EOF
    writer.write_all(&[0;8]).unwrap();

}

/// swaps endianness of the byte-vector
/// assuming 8byte (u64) words
/// simply grabs each 8byte section and reverses the order within
/// # params
/// * wordsize: 4 (bytes) for u32, 8(bytes) for u64
fn swap_endian(bytes: &[u8], wordsize: usize) -> Vec<u8>{
    let mut swapped_endian: Vec<u8> = Vec::with_capacity(bytes.len());
    for bytes in bytes.chunks(wordsize){
        swapped_endian.extend(bytes.iter().rev());
    }
    swapped_endian
}

pub fn decompress_busfile(input: &str, output: &str) {

    let reader = BuszReader::new(input);
    let mut writer = BusWriter::new(
        output,
        reader.bus_header.clone()
    );
    println!("{:?}", reader);

    for r in reader {
        writer.write_record(&r);
    }
}

// all the encoded parts (cb,umi,ec...) are padded with zeros until a mutiple of 64
// which we need to remove
fn calc_n_trailing_bits(bits_processed: usize) -> usize {
    let padded_size = round_to_multiple(bits_processed, 64);
    let zeros_toremoved = padded_size - bits_processed;
    zeros_toremoved
}

// // turn a bitvector (64 elements) into a u64
pub fn bits64_to_u64(x: bit_vec::BitVec) -> u64{
    assert_eq!(x.len(), 64);
    u64::from_le_bytes(x.to_bytes()[..8].try_into().unwrap())
}


const BUSZ_HEADER_SIZE: usize = 4+4+4;
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
pub struct BuszHeader {
    block_size: u32,
    pfd_block_size: u32,
    lossy_umi: u32,
}
impl BuszHeader {
    /// desearializes a BusHeader from Bytes; when reading busfiles
    /// assumes Little-Endian! (https://docs.rs/bincode/latest/bincode/config/index.html#options-struct-vs-bincode-functions)
    pub fn from_bytes(bytes: &[u8]) -> BuszHeader {
        let header_struct: BuszHeader =
            // this interprets the bytes in Little Endian!, i.e bytes=[1,0,0,0,0,0,0,0] = 1_u64
            bincode::deserialize(bytes).expect("FAILED to deserialze busz header");
        assert_eq!(
            header_struct.lossy_umi, 0,
            "lossy_umi != 0 not supported"
        );
        header_struct
    }
    /// seialize the header to bytes
    /// assumes Little-Endian! (https://docs.rs/bincode/latest/bincode/config/index.html#options-struct-vs-bincode-functions)
    pub fn to_bytes(&self) -> Vec<u8> {
        bincode::serialize(self).expect("FAILED to serialze header")
    }
}

#[derive(Debug)]
pub struct BuszReader {
    pub bus_header: BusHeader,
    pub busz_header: BuszHeader,
    pub reader: BufReader<File>,
    buffer: VecDeque<BusRecord>
}

impl BuszReader {
    /// main constructor for busreader, buffersize is set to best performance
    pub fn new(filename: &str) -> Self {
        BuszReader::new_with_capacity(filename, DEFAULT_BUF_SIZE)
    }

    /// Creates a buffered reader over busfiles, with specific buffersize
    pub fn new_with_capacity(filename: &str, bufsize: usize) -> Self {

        // read the regular busheader
        let bus_header = BusHeader::from_file(filename);
        println!("Reader: bus header {:?}", bus_header);
        let mut file_handle = File::open(filename).expect("FAIL");

        // advance the file point 20 bytes (pure header) + header.tlen (the variable part)
        let to_seek = BUS_HEADER_SIZE as u64 + bus_header.get_tlen() as u64;
        let _x = file_handle.seek(SeekFrom::Start(to_seek)).unwrap();

        // now we're at the start of the BusZHeader
        let mut buszheader_bytes = [0_u8; BUSZ_HEADER_SIZE];
        file_handle.read_exact(&mut buszheader_bytes).unwrap();
        // println!("{:?}", buszheader_bytes);

        let bzheader = BuszHeader::from_bytes(&buszheader_bytes);
        println!("Reader: busz header {:?}", bzheader);
        //FH is now at the start of the actual contents
        let buf = BufReader::with_capacity(bufsize, file_handle);

        let buffer = VecDeque::with_capacity(bzheader.block_size as usize);
        BuszReader { bus_header, busz_header:bzheader, reader: buf, buffer }
    }

    /// takes the next 8 bytes (u64) out of the stream and interprets it as a busheader
    /// if we encounter the zero byte [00000000], this indicates the EOF
    fn load_busz_header(&mut self) -> Option<CompressedBlockHeader>{
        // get block header
        let mut blockheader_bytes = [0_u8;8];
        self.reader.read_exact(&mut blockheader_bytes).unwrap();
        if blockheader_bytes == [0,0,0,0,0,0,0,0] {  // EOF
            return None
        }
        let h = CompressedBlockHeader { header_bytes: u64::from_le_bytes(blockheader_bytes)};
        // println!("H bytes {}, H records {}", H.get_blocksize_and_nrecords().0, H.get_blocksize_and_nrecords().1);
        Some(h)
    }

    pub fn load_busz_block_faster(&mut self) -> Option<Vec<BusRecord>>{

        let h: Option<CompressedBlockHeader> = self.load_busz_header();
        if h.is_none(){
            return None
        }
        let header = h.unwrap();
        let (block_size_bytes, nrecords) = header.get_blocksize_and_nrecords();

        println!("BusZ block-header bytes:{block_size_bytes} #records{nrecords}");
        // read the busZ-block body
        let mut block_buffer: Vec<u8> = Vec::from_iter((0..block_size_bytes).map(|x| x as u8));
        self.reader.read_exact(&mut block_buffer).unwrap();

        
        // conversion to big-endian to make the reading work right
        let bigendian_buf = swap_endian(&block_buffer, 8);
        let theblock: &bv::BitVec<u8, bv::Msb0> = &bv::BitVec::from_slice(&bigendian_buf);
        let mut block = BuszBlock::new(theblock.as_bitslice(),nrecords as usize);
        
        let records =block.parse_block();
    
        Some(records)
    }
}

impl Iterator for BuszReader {
    type Item=BusRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer.is_empty(){
            println!("buffer empty, loading new block");
            // pull in another block
            let block =self.load_busz_block_faster(); 
            match block {
                Some(records) => {self.buffer = records.into()},
                None => {return None}
            };
        }
        match self.buffer.pop_front() {
            Some(record) => Some(record),
            None => panic!("cant happen")
        }
    }
}


/// BusZ format has the problem that it can only write to disk
/// once `busz_blocksize` elements have been encountered.
/// The last block might be truncated (but not intermediates!)
/// 
/// Due to this, we have to buffer the data. Yet when done writing,
/// we need to ensure the buffer gets emptied -> truncated block
/// Also we need to add an EOF tag.
/// 
/// Now we want to avoid writing anything else after we've flushed/closed the writer
/// Hence the State to keep track of
#[derive(Eq, PartialEq, Debug)]
enum BuszWriterState {
    Open,
    FlushedAndClosed,
}

/// Writing BusRecords into compressed .busz format
/// 
/// # Important:
/// One **ALWAYS** has to call `terminal_flush()` in the end to properly terminate the busz file.
/// Otherwise the .busz won't be readable and some entries might not make it to disk
/// 
/// The safest option is to rely on the `write_iterator` method. 
/// This writes the entire iterator to file with proper termination.
/// 
pub struct BuszWriter {
    writer: BufWriter<File>,
    internal_buffer: Vec<BusRecord>, // store items until we got enough to write a new block
    busz_blocksize: usize,
    header: BusHeader,
    state: BuszWriterState
}

impl BuszWriter {
    /// create a BuszWriter from an open Filehandle
    pub fn from_filehandle(file_handle: File, header: BusHeader, busz_blocksize: usize) -> BuszWriter {
        BuszWriter::new_with_capacity(file_handle, header, busz_blocksize)
    }

    /// create a Buszwriter that streams records into a file
    pub fn new(filename: &str, header: BusHeader, busz_blocksize: usize) -> BuszWriter {
        let file_handle: File = File::create(filename).expect("FAILED to open");
        BuszWriter::from_filehandle(file_handle, header, busz_blocksize)
    }

    pub fn write_iterator(&mut self, iter: impl Iterator<Item=BusRecord>) {
        for chunk in &iter.chunks(self.busz_blocksize) {
            let records: Vec<BusRecord> = chunk.collect();
            let compressed_block = compress_busrecords_into_block(&records);
            self.writer.write_all(&compressed_block).unwrap();
        }
    
        // write a final, emtpy block header (0_u64) that signals the EOF
        self.writer.write_all(&[0;8]).unwrap();
        self.writer.flush().unwrap();
        self.state = BuszWriterState::FlushedAndClosed;
    }


    /// Writing a single BusRecord
    pub fn write_record(&mut self, record: BusRecord) {

        if self.state == BuszWriterState::FlushedAndClosed {
            panic!("Buffer has been flushed and closed!")
        }

        if self.internal_buffer.len() < self.busz_blocksize - 1{  // adding an element, but it doesnt fill the buffer
            self.internal_buffer.push(record);
            println!("adding to buffer: {}", self.internal_buffer.len());
        }
        else { // buffer is full including this element
            self.internal_buffer.push(record);
            println!("buffer full: {}", self.internal_buffer.len());
            self.write_buffer_to_disk()
        }
    }

    fn write_buffer_to_disk(&mut self) {
        if self.state == BuszWriterState::FlushedAndClosed {
            panic!("Buffer has been flushed and closed!")
        }
        println!("Writing to disk {} entires", self.internal_buffer.len());
        let compressed_block = compress_busrecords_into_block(&self.internal_buffer);
        self.writer.write_all(&compressed_block).unwrap();
        self.internal_buffer.clear();

        // f;lush the writer too
        self.writer.flush().unwrap();
    }

    /// Writing multiple BusRecords at once
    pub fn write_records(&mut self, records: Vec<BusRecord>) {
        if self.state == BuszWriterState::FlushedAndClosed {
            panic!("Buffer has been flushed and closed!")
        }

        // writes several recordsd and flushes
        for r in records {
            self.write_record(r)
        }
    }

    /// this needs to be called when we're done writing to the buszfile
    /// empties the internal buffer, adds the EOF
    pub fn terminal_flush(&mut self) {

        if !self.internal_buffer.is_empty() {
            self.write_buffer_to_disk();
        }

        // write a final, emtpy block header (0_u64) that signals the EOF
        self.writer.write_all(&[0;8]).unwrap();
        self.writer.flush().unwrap();

        self.state = BuszWriterState::FlushedAndClosed;

    }

    /// creates a buszwriter with specified buffer capacity (after how many records an actual write happens)
    /// dont use , 800KB is the default buffer size and optimal for performance
    pub fn new_with_capacity(file_handle: File, header: BusHeader, busz_blocksize: usize) -> Self {
        let mut writer = BufWriter::with_capacity(DEFAULT_BUF_SIZE, file_handle);

        // write the header into the file
        // write the header into the file
        let magic: [u8; 4] = *b"BUS\x01";
        let busheader = BusHeader{
            magic, 
            version: header.version, 
            cb_len: header.cb_len, 
            umi_len: header.umi_len, 
            tlen: 0 //header.tlen
        };
        let binheader = busheader.to_bytes();
        writer
            .write_all(&binheader)
            .expect("FAILED to write header");

        // write the variable header
        let mut varheader: Vec<u8> = Vec::new();
        for _i in 0..busheader.tlen {
            varheader.push(0);
        }
        writer
            .write_all(&varheader)
            .expect("FAILED to write var header");
        
        // BusZ header
        let busz_header = BuszSpecificHeader {
            block_size: busz_blocksize.try_into().unwrap(),
            pfd_block_size: 512,
            lossy_umi: 0
        };
        let binzheader = busz_header.to_bytes();
        writer
            .write_all(&binzheader)
            .expect("FAILED to write header");

        let internal_buffer: Vec<BusRecord> = Vec::with_capacity(busz_blocksize);

        BuszWriter { writer, internal_buffer, busz_blocksize, header, state: BuszWriterState::Open }
    
    }
}


/// to keep track in which section we are in the busz-block
#[derive(Debug, PartialEq, Eq)]
enum BuszBlockState {
    // Header,
    Cb,
    Umi,
    Ec,
    Count,
    Flag,
    Finished
}

/// A single block of a busz-file, header has already been parsed
/// Contains all bytes from the block (we know how many bytes from the header)
/// and parses those bytes into BusRecords
#[derive(Debug)]
struct BuszBlock <'a> {
    buffer: &'a bv::BitSlice<u8, bv::Msb0>,
    pos: usize, // where we are currently in the buffer
    n_elements: usize ,// how many busrecords are stored in the block
    state: BuszBlockState,
    debug: bool
}

impl <'a> BuszBlock <'a> {
    // pub fn new(buffer: &'a BitSlice<u8, Msb0>, n_elements: usize) -> Self {
    pub fn new(buffer: &'a bv::BitSlice<u8, bv::Msb0>, n_elements: usize) -> Self {
        // TODO: warning, buffer must be conveted in a special way if the bytes come out of a file
        // see BuszReader::load_busz_block_faster
        // cant do it in herer due to lifetime issues
        BuszBlock { 
            buffer, 
            pos: 0, 
            n_elements, 
            state: BuszBlockState::Cb,
            debug: false
        }
    }

    fn parse_cb(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Cb);
        assert_eq!(self.pos, 0, "must be at beginning of buffer to parse CBs");

        //decode the CBs
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
        // println!("Decoding CBs {:?}", block_body);

        let cb_buffer = &self.buffer[self.pos..];
        // let bits_before = self.buffer.len();

        let mut fibdec = newpfd::fibonacci::FibonacciDecoder::new(cb_buffer);
        // let mut fibdec = FibbonacciDecoder { bitstream: block_body};

        let CB_RLE_VAL = 0 + 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut cb_delta_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            // println!("Counter {counter}");
            // println!("Bits {:?}", fibdec.bitstream);
            // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
            let item = fibdec.next().unwrap();
            // println!("Item {item}");
            if item == CB_RLE_VAL {
                let runlength = fibdec.next().unwrap();
                // let (runlength, remainder) = bitvec_single_fibbonacci_decode(&mut block_body);
                // block_body = remainder;  //weird scope thing: if we call the above with return named block_body, it'll be its own var and not progagate
                // println!("Decoding run of  {runlength}");
                // println!("Body after   {:?}", fibdec.bitstream);

                for _ in 0..runlength{
                    cb_delta_encoded.push(CB_RLE_VAL -1 );  //due to shift+1
                    counter+= 1;
                }
            } else {
                // println!("Decoding sigle element");
                cb_delta_encoded.push(item - 1);//due to shift+1
                counter+= 1;
            }
        }
        // undo delta encoding, i,e. sumsum
        cb_delta_encoded.iter_mut().fold(0, |acc, x| {
            *x += acc;
            *x
        });
        // println!("cb_delta_encoded: {:?}", cb_delta_encoded);
        // println!("decoded {:?} CBs", cb_delta_encoded.len());
        // println!("Bitstream after CB {:?}", fibdec.bitstream);

        // need to remove the trailing 0 which come from padding to the next mutiple of u64 in th CB block
        // we know that the CB block has a length of n x 64
        let bits_processed = fibdec.get_bits_processed();

        // let padded_size = round_to_multiple(bits_processed, 64);
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
        // println!("before {bits_before}, after {bits_after}; zeros to remove: {zeros_toremoved}");
        // println!("CBs: {}", cb_delta_encoded.len());
        // println!("Cbs: bits proc {}, bits padded {}", bits_processed, zeros_toremoved);
        // println!("CBs: buffer at {}", self.pos + bits_processed + zeros_toremoved);

        assert!(
            !self.buffer[self.pos + bits_processed..self.pos + bits_processed + zeros_toremoved].any()
        );

        // adcance the pointer to the start of the next section, i.e. the umis
        self.pos = self.pos + bits_processed + zeros_toremoved;
        self.state = BuszBlockState::Umi;

        cb_delta_encoded
    }

    fn parse_umi(&mut self, cbs: &[u64]) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Umi);

        let umi_buffer = &self.buffer[self.pos..];
        let mut fibdec = newpfd::fibonacci::FibonacciDecoder::new(umi_buffer);

        // =====================================
        // UMIs

        let UMI_RLE_VAL = 0 + 1 ; // since shifted
        let mut counter = 0;

        // guarantee that last_barcode != rows[0].barcode
        let mut last_cb = cbs[0] + 1;
        let mut umi =0_u64;
        let mut umis: Vec<u64> = Vec::with_capacity(self.n_elements);
        while counter < self.n_elements {
            let diff = fibdec.next().unwrap() - 1;
            // println!("Decoded item {diff}");
            let current_cb =  cbs[counter];
            if last_cb !=current_cb {
                umi=0;
            }
            // println!("Item {item}");
            if diff == UMI_RLE_VAL - 1 {
                let runlength = fibdec.next().unwrap();
                // let (runlength, remainder) = bitvec_single_fibbonacci_decode(&mut block_body);
                // block_body = remainder;  //weird scope thing: if we call the above with return named block_body, it'll be its own var and not progagate
                // println!("Decoding run of  {runlength}");
                // println!("Body after   {:?}", fibdec.bitstream);

                for _ in 0..runlength{
                    umis.push(umi -1 );  //due to shift+1
                    counter+= 1;
                }
            } else {
                umi+= diff;
                // println!("Decoding sigle element");
                umis.push(umi - 1);//due to shift+1
                counter+= 1;
            }
            last_cb = current_cb;
        }
        // println!("UMIs {:?}", umis);
        // println!("decoded {:?} UMIs", umis.len());

        // need to remove the trailing 0 which come from padding to the next mutiple of u64
        // we know that the CB block has a length of n x 64
        let bits_processed =  fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
        // println!("UMIs: {}", umis.len());
        // println!("UMIs: bits proc {}, bits padded {}", bits_processed, zeros_toremoved);
        // println!("UMIs: buffer at {}", self.pos + bits_processed + zeros_toremoved);

        assert!(
            !self.buffer[self.pos + bits_processed..self.pos + bits_processed + zeros_toremoved].any()
        );

        // adcance the pointer to the start of the next section, i.e. the umis
        self.pos = self.pos + bits_processed + zeros_toremoved;
        self.state = BuszBlockState::Ec;

        umis

    }

    fn parse_ec(&mut self) -> Vec<u64> {
        // here it gets tricky: previously we swapped the entire stream
        // to u64-little endian
        //
        // however, ECs uses u32 as a basis of operatiors
        // and in order to get the bits appear in the right order we need u32-little endian
        // 1. revert the u64 endian swap
        // 2. swap to u32 little endian
        // 3. do decoding
        // 4 swap back to u64-endian

        let ec_buffer = &self.buffer[self.pos..];

        // println!("converting EC bitstream");
        let bytes = bitslice_to_bytes(ec_buffer);
        let _original = swap_endian(&bytes, 8);
        let little_endian_32_bytes = swap_endian(&_original, 4);
       
        let remainder_little_endian_32: &bv::BitSlice<u8, bv::Msb0> =  bv::BitSlice::from_slice(&little_endian_32_bytes);

        // println!("len remainder_little_endian_32 {}", remainder_little_endian_32.len());
        // println!("remainder_little_endian_32\n{}", bitstream_to_string(remainder_little_endian_32));
        // println!("main buf: remainder_little_endian_32\n{}", bitstream_to_string(ec_buffer));

        let (ecs, bits_consumed) = newpfd::newpfd_bitvec::decode(remainder_little_endian_32, self.n_elements, PFD_BLOCKSIZE);

        // println!("EC buffer used:\n{}", bitstream_to_string(&remainder_little_endian_32[..bits_consumed]));

        self.pos += bits_consumed;
        // println!("after EC main buffer: {}", bitstream_to_string(&self.buffer[self.pos..]));
        // println!("after EC main buffer: {}", bitstream_to_string(&ec_buffer[bits_consumed..]));

        // println!("ECs: {}", ecs.len());
        // println!("ECs: bits proc {}, bits padded 0", bits_consumed);
        // println!("ECs: buffer at {}", self.pos);


        self.state = BuszBlockState::Count;

        ecs
    }

    fn parse_counts(&mut self) ->Vec<u64> {
        // ===============================
        // count decoding
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
        // println!("Decoding Counts {:?}", remainder_little_endian_64);
        assert_eq!(self.state, BuszBlockState::Count);

        // things get complicated here: in EC-we had to swap_endian and accomodate u32 encoding
        // this was done within `parse_ec` and did not affect the self.buffer
        // Now we're back to u64. However somethings off. In particular, trying to same procedure 
        // as with the old `decompress_busz_block` (undoing the u32 swap from ec) doesnt work
        // Its an alignment issue!! The EC block might end in the middle of a u64:
        //       ----------------------------u64-----------------------------|--------------------------u64------------
        // `aaaaaaaabbbbbbbbccccccccdddddddd|eeeeeeeeffffffffgggggggghhhhhhhh|iiiiiiiijjjjjjjjkkkkkkkkllllllllmmmmmmmmnnnnnnnnoooooooopppppppp`
        // `EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE|CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC|CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC`  // E= EC, C=COUNT
        // 
        //initially stuf has been swapped around using the whole buffer
        // now, trying to undo the swapping on only the COUNT-part wont work
        // The original byte order was `hgfedcba|ponmlkji` -> endian-swap-64 -> abcdefgh|ijklmnop
        // undoing that swap on the count part only: efgh|ijklmnop (without taking into accont the alignment) will lead to lkjihgfe|....ponm|
        //
        // Solution: operate swaping on the entire buffer, THEN select the COUNTs part

        // undoing initial flip
        let bytes = bitslice_to_bytes(self.buffer);
        // this is what happens in EC: undo u64-swap, apply u32-swap
        let a = swap_endian(&bytes, 8);
        let b = swap_endian(&a, 4);
        // move the bufferposition after EC
        assert_eq!(self.pos % 8, 0);
        let ec_pos = self.pos / 8;

        let bytes_behind_ec = &b[ec_pos..];
        let c = swap_endian(bytes_behind_ec, 4);
        let d = swap_endian(&c, 8);
        // println!("bytes\n{:?}", bytes);
        // println!("after bytes\n{:?}", d);
        let count_buffer: &bv::BitSlice<u8, bv::Msb0> =  bv::BitSlice::from_slice(&d);

        // println!("Count buffer, reswapped\n{}", bitstream_to_string(count_buffer));


        // let mut count_buffer = &self.buffer[self.pos..];
        // println!("before transform count_buffer\n{}", bitstream_to_string(count_buffer));

        let mut fibdec = newpfd::fibonacci::FibonacciDecoder::new(count_buffer);

        let COUNT_RLE_VAL = 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut counts_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            let item = fibdec.next().unwrap();
            // println!("Item {item}");
            if item == COUNT_RLE_VAL {
                let runlength = fibdec.next().unwrap();
                for _ in 0..runlength{
                    counts_encoded.push(COUNT_RLE_VAL);  //due to shift+1
                    counter+= 1;
                }
            } else {
                counts_encoded.push(item);//due to shift+1
                counter+= 1;
            }
        }

        let bits_processed =  fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;

        // println!("COUNT: {}", counts_encoded.len());
        // println!("COUNT: bits proc {}, bits padded {}", bits_processed, zeros_toremoved);
        // println!("COUNT: buffer at {}", self.pos + bits_processed + zeros_toremoved);

        assert!(
            !count_buffer[bits_processed..bits_processed + zeros_toremoved].any()
        );

        // adcance the pointer to the start of the next section, i.e. the umis
        self.pos = self.pos + bits_processed + zeros_toremoved;
        self.state = BuszBlockState::Flag;

        counts_encoded
    }

    fn parse_flags(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Flag);

        // same issue as in count: some misalginemtn and endianess issues. 
        // do the same thing: transformations on the entire buffer, the move position to FLAG section
        // undoing initial flip

        let bytes = bitslice_to_bytes(self.buffer);
        // this is what happens in EC: undo u64-swap, apply u32-swap
        let a = swap_endian(&bytes, 8);
        let b = swap_endian(&a, 4);
        // move the bufferposition after EC
        assert_eq!(self.pos % 8, 0);
        let count_pos = self.pos / 8;

        let bytes_behind_count = &b[count_pos..];
        let c = swap_endian(bytes_behind_count, 4);
        let d = swap_endian(&c, 8);
        // println!("after bytes\n{:?}", d);
        let flag_buffer: &bv::BitSlice<u8, bv::Msb0> =  bv::BitSlice::from_slice(&d);
        // println!("flag_buffer\n{}", bitstream_to_string(flag_buffer));


        let mut fibdec = newpfd::fibonacci::FibonacciDecoder::new(flag_buffer);

        let FLAG_RLE_VAL = 0+1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut flag_decoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            // println!("Counter {counter}");
            // println!("Bits {:?}", fibdec.bitstream);
            // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
            let item = fibdec.next().unwrap();
            // println!("Item {item}");
            if item == FLAG_RLE_VAL {
                let runlength = fibdec.next().unwrap();
                // let (runlength, remainder) = bitvec_single_fibbonacci_decode(&mut block_body);
                // block_body = remainder;  //weird scope thing: if we call the above with return named block_body, it'll be its own var and not progagate
                // println!("Decoding run of  {runlength}");
                // println!("Body after   {:?}", fibdec.bitstream);
                for _ in 0..runlength{
                    flag_decoded.push(FLAG_RLE_VAL - 1 );  //due to shift+1
                    counter+= 1;
                }
            } else {
                // println!("Decoding sigle element {}", item-1 );
                flag_decoded.push(item-1);//due to shift+1
                counter+= 1;
            }
        }
        // println!("Falgs: {:?}", flag_decoded);
        // println!("Reamining: {:?}", fibdec.bitstream);

        let bits_processed = fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
        
        assert!(
            !flag_buffer[bits_processed..bits_processed + zeros_toremoved].any()
        );
        // adcance the pointer to the start of the next section, i.e. the umis
        self.pos = self.pos + bits_processed + zeros_toremoved;
        self.state = BuszBlockState::Finished;
        
        assert_eq!(self.pos, self.buffer.len(), "still leftover bits in the buffer!"); 

        flag_decoded
    }

    pub fn parse_block(&mut self) -> Vec<BusRecord>{
        if self.debug {
            println!("Block-bits:\n{}", bitstream_to_string(&self.buffer[self.pos..]));

            println!("CB pos {}", self.pos);
            }
        let cbs = self.parse_cb();

        if self.debug {
            println!("CBs: {:?}", cbs);

            println!("UMI pos {}", self.pos);
            }
        // println!("umi: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let umis = self.parse_umi(&cbs);


        if self.debug {
            println!("UMIs: {:?}", umis);

            println!("EC pos {}", self.pos);
            }
        // println!("ec: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let ec = self.parse_ec();


        if self.debug {
            println!("ECs: {:?}", ec);

            println!("COUNT pos {}", self.pos);
            println!("COUNT buf\n{}", bitstream_to_string(&self.buffer[self.pos..]));
            }
        // println!("count: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let count = self.parse_counts();


        if self.debug {
            println!("count: {:?}", count);

            println!("FLAG pos {}", self.pos);
            }
        // println!("flag: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let flag = self.parse_flags();
        if self.debug {
            println!("flag: {:?}", flag);
        }
        // println!("{:?}", cbs);
        // println!("{:?}", umis);
        // println!("{:?}", ec);
        // println!("{:?}", count);
        // println!("{:?}", flag);

        let mut decoded_records = Vec::with_capacity(self.n_elements);
        for (cb,umi,ec,count,flag) in izip!(cbs, umis, ec, count, flag) {
            decoded_records.push(
                BusRecord {
                    CB: cb, 
                    UMI: umi, 
                    EC:ec.try_into().unwrap(), 
                    COUNT:count.try_into().unwrap(), 
                    FLAG:flag.try_into().unwrap()
                }
            );
        }
        decoded_records
    }
}

/// turn a bitslice into an array of bytes
/// the first 8 bits (bits[..8]) will become the first byte in the result
/// i.e. a sort of BigEndian encoding
fn bitslice_to_bytes(bits: &bv::BitSlice<u8, bv::Msb0>) -> Vec<u8>{

    assert_eq!(bits.len() % 8,  0, "cant covnert to bytes if Bitsclie is not a multiple of 8");

    let nbytes = bits.len() / 8;
    let mut bytes = Vec::with_capacity(nbytes);
    for i in 0..nbytes {
        let b = &bits[i*8.. (i+1)*8];
        let a: u8 = b.load_be(); // doesnt matter if be/le here since its only 1byte anyway
        bytes.push(a);
    }
    bytes
}


fn bitstream_to_string(buffer: &bv::BitSlice<u8, bv::Msb0>) -> String{
    let mut s = String::new();
    let x = buffer.iter().map(|x| if *x{"1"} else {"0"});
    for bit64 in &x.into_iter().chunks(64){
        let concat = bit64.collect::<Vec<_>>().join("");
        s.push_str(&concat);
        s.push('\n');
    }
    s
}

#[cfg(test)]
mod test {
    use crate::{busz::{CompressedBlockHeader,decompress_busfile, swap_endian}};
    use super::{BuszReader, compress_busfile};

    #[test]
    fn endian_swapping() {
        let v = vec![0_u8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
        let a = swap_endian(&v, 8);
        let b = swap_endian(&a, 4);
        let c = swap_endian(&b, 4);
        let d = swap_endian(&c, 8);
        assert_eq!(v,d);
    }

    #[test]
    fn test_header_encode_decode() {
        let nbytes = 20;
        let nrecords = 10;
        let h = CompressedBlockHeader::new(nbytes, nrecords);

        assert_eq!(h.get_blocksize_and_nrecords().0, nbytes);
        assert_eq!(h.get_blocksize_and_nrecords().1, nrecords);
    }
    mod encode {
        use crate::{io::BusRecord, busz::{compress_barcodes2, compress_umis, compress_ecs}};

        #[test]
        fn test_cb_encode(){
            let v = vec![ 
                BusRecord {CB:0,UMI:0,EC:0,COUNT:1, FLAG: 0 },
                BusRecord {CB:0,UMI:0,EC:1,COUNT:1, FLAG: 0 },
                BusRecord {CB:1,UMI:0,EC:2,COUNT:1, FLAG: 0 },
                BusRecord {CB:1,UMI:0,EC:3,COUNT:1, FLAG: 0 },
                BusRecord {CB:1,UMI:0,EC:4,COUNT:1, FLAG: 0 },
                BusRecord {CB:1,UMI:0,EC:5,COUNT:1, FLAG: 0 },
            ];
            let enc = compress_barcodes2(&v);
            let decoded: Vec<_> = fibonacci_codec::fib_decode_u64(enc).map(|x|x.unwrap()).collect();
            
            assert_eq!(decoded, vec![
                1,2,  // two zers
                2,    // a single 1
                1,3   // three zero
                ]);
        }

        #[test]
        fn test_umi_encode(){
            let v = vec![ 
                BusRecord {CB:0,UMI:0,EC:0,COUNT:1, FLAG: 0 },
                BusRecord {CB:0,UMI:0,EC:1,COUNT:1, FLAG: 0 },
            ];
            let enc = compress_umis(&v);
            let decoded: Vec<_> = fibonacci_codec::fib_decode_u64(enc).map(|x|x.unwrap()).collect();
            
            assert_eq!(decoded, vec![
                2,  // a one (the refence umi is zero, but all are incemrented by one)
                1,1 //a 1-run
                ]);
        }
        #[test]
        fn test_umi_encode_wraparound(){
            let v = vec![ 
                BusRecord {CB:0,UMI:10,EC:0,COUNT:1, FLAG: 0 },   // 10
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:1,UMI:0,EC:1,COUNT:1, FLAG: 0 },    // 1
            ];
            let enc = compress_umis(&v);
            let decoded: Vec<_> = fibonacci_codec::fib_decode_u64(enc).map(|x|x.unwrap()).collect();
            
            assert_eq!(decoded, vec![
                12,  // a 10 (the refence umi is zero, but all are incemrented by one)
                1,2, //a 2-run
                2
                ]);
        }

        #[test]
        fn test_ec_encode(){
            let v = vec![ 
                BusRecord {CB:0,UMI:10,EC:0,COUNT:1, FLAG: 0 },   // 10
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:1,UMI:0,EC:1,COUNT:1, FLAG: 0 },    // 1
            ];
            let enc = compress_ecs(&v);
            println!("enc size {}", enc.len());
        }   
    }
    
    mod read_write {
        use crate::{io::{BusRecord, BusHeader}, busz::{BuszWriter, BuszReader, BuszWriterState}};
        #[test]
        fn test_busz_writer_good_blocksize() {
            let blocksize = 2;  // blocksize wont fit all elements

            let v = vec![ 
                BusRecord {CB:0,UMI:10,EC:0,COUNT:1, FLAG: 0 },   // 10
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:1,UMI:0,EC:1,COUNT:1, FLAG: 0 },    // 1
            ];

            let buszfile = "/tmp/busz_writer_test.buzs";
            let header = BusHeader::new(16,12,0);
            let mut writer = BuszWriter::new(buszfile, header , blocksize);
            writer.write_records(v.clone());
            writer.terminal_flush();

            assert_eq!(writer.state, BuszWriterState::FlushedAndClosed);
            assert_eq!(writer.internal_buffer.len(), 0);

            let reader = BuszReader::new(buszfile);
            assert_eq!(reader.collect::<Vec<_>>(), v);
        }

        #[test]
        fn test_busz_writer_crooked_blocksize() {
            let blocksize = 3;  // blocksize wont fit all elements

            let v = vec![ 
                BusRecord {CB:0,UMI:10,EC:0,COUNT:1, FLAG: 0 },   // 10
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:1,UMI:0,EC:1,COUNT:1, FLAG: 0 },    // 1
            ];

            let buszfile = "/tmp/busz_writer_test.buzs";
            let header = BusHeader::new(16,12,0);
            let mut writer = BuszWriter::new(buszfile, header , blocksize);
            writer.write_records(v.clone());
            writer.terminal_flush();

            assert_eq!(writer.state, BuszWriterState::FlushedAndClosed);
            assert_eq!(writer.internal_buffer.len(), 0);

            let reader = BuszReader::new(buszfile);
            assert_eq!(reader.collect::<Vec<_>>(), v);
        }
        #[test]
        fn test_busz_write_iterator() {
            let blocksize = 3;  // blocksize wont fit all elements

            let v = vec![ 
                BusRecord {CB:0,UMI:10,EC:0,COUNT:1, FLAG: 0 },   // 10
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:1,UMI:0,EC:1,COUNT:1, FLAG: 0 },    // 1
            ];
            let buszfile = "/tmp/busz_writer_test.buzs";
            let header = BusHeader::new(16,12,0);
            let mut writer = BuszWriter::new(buszfile, header , blocksize);
            writer.write_iterator(v.iter().cloned());


            assert_eq!(writer.state, BuszWriterState::FlushedAndClosed);
            assert_eq!(writer.internal_buffer.len(), 0);
            
            let reader = BuszReader::new(buszfile);

            assert_eq!(reader.collect::<Vec<_>>(), v);
        }
    }
    
    mod external {
        use crate::{io::{BusRecord, BusWriter, BusHeader, BusReader}, busz::{BuszReader, compress_busfile, decompress_busfile}};

        #[test]
        fn test_external(){
            let v = vec![ 
                BusRecord {CB:10,UMI:11,EC:10,COUNT:13, FLAG: 14 },   // 10
                BusRecord {CB:11,UMI:11,EC:10,COUNT:13, FLAG: 14 },   // 0
                BusRecord {CB:22,UMI:10,EC:10,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:22,UMI:11,EC:10,COUNT:1, FLAG: 0 },    // 1
            ];
            let mut  writer = BusWriter::new(
                "/tmp/buscompress.bus", 
                BusHeader::new(16, 12, 1)
            );
            writer.write_records(&v);
        }

        #[test]
        fn test_encode_decode_busz(){

            println!("Writing plain busfile");
            let v = vec![ 
                BusRecord {CB:10,UMI:11,EC:10,COUNT:13, FLAG: 20 },   // 10
                BusRecord {CB:11,UMI:11,EC:10,COUNT:13, FLAG: 20 },   // 0
                BusRecord {CB:22,UMI:10,EC:10,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:22,UMI:11,EC:10,COUNT:1, FLAG: 0 },    // 1
            ];
            let input_plain = "/tmp/buscompress.bus";
            let mut  writer = BusWriter::new(
                input_plain, 
                BusHeader::new(16, 12, 1)
            );
            writer.write_records(&v);
            drop(writer);

            let copmressed_output = "/tmp/lalalala.busz";
            println!("copmressing busfile");
            compress_busfile(
                input_plain,
                copmressed_output,
                100
            );

            println!("decoding busfile");

            // // decode it
            let reader = BuszReader::new("/tmp/lalalala.busz");
            let recs: Vec<_> = reader.collect();
            assert_eq!(v, recs);

        }

        #[test]
        fn test_encode_decode_busz_biggerfile(){

            let input_plain = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";
            let input_compressed_true = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";

            let copmressed_output = "/tmp/output.corrected.sort.busz";
            println!("copmressing busfile");
            compress_busfile(
                input_plain,
                copmressed_output,
                10000
            );
            println!("decoding busfile");

            // // decode it
            let reader = BuszReader::new("/tmp/output.corrected.sort.busz");
            // reader.load_busz_block_faster();
            let recs: Vec<_> = reader.collect();

            let x = BusReader::new(input_compressed_true);
            assert_eq!(x.collect::<Vec<_>>(), recs);

        }


        #[test]
        fn test_compress1() {
            let input_compressed = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.busz"; 
            let input_plain = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";
            let copmressed_output = "/tmp/buscompress_testing.busz";
            compress_busfile(
                input_plain,
                copmressed_output,
                100
            );
        }


        #[test]
        fn test_compress_full() {
            let input_compressed = "/home/michi/bus_testing/bus_output/output.corrected.sort.busz"; 
            let input_plain = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";
            let copmressed_output = "/tmp/buscompress_testing_full.busz";
            compress_busfile(
                input_plain,
                copmressed_output,
                10000
            );
        }

        #[test]
        fn test_compress2() {
            let mut reader = BuszReader::new("/tmp/lalalala_true.busz");
            println!("{:?}", reader.bus_header);
            reader.next();
            println!("==========================================");
            println!("==========================================");
            let mut reader = BuszReader::new("/tmp/lalalala.busz");
            println!("{:?}", reader.bus_header);

            reader.next();
            // let records:Vec<_> = reader.collect();
            // println!("{}", records.len())
        }
        #[test]
        fn test_decompress(){
            // decompress a busfile, check that the contents match the true (uncompressed version)
            let input_compressed = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.busz"; 
            let input_plain = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";
            decompress_busfile(
                input_compressed,
                "/tmp/buscompress_lala.bus");
            let r = BusReader::new("/tmp/buscompress_lala.bus");
            let records:Vec<_> = r.collect();
            println!("{} compressed records read", records.len());

            let r_original = BusReader::new(input_plain);
            let records_original:Vec<_> = r_original.collect();

            assert_eq!(records.len(), records_original.len());
            assert_eq!(records, records_original);
        }

        #[test]
        fn test_iterator(){
            let reader = BuszReader::new("/home/michi/bus_testing/bus_output_shortest/output.corrected.sort.busz");
            let records:Vec<_> = reader.collect();

            let r_original = BusReader::new("/home/michi/bus_testing/bus_output_shortest/output.corrected.sort.bus");
            let records_original:Vec<_> = r_original.collect();

            assert_eq!(records.len(), records_original.len());
            assert_eq!(records, records_original);
        }   
    }

}
