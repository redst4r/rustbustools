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
//! ```
//! aaaaaaaabbbbbbbbccccccccddddddddeeeeeeeeffffffffgggggggghhhhhhhh  //bits
//! ```
//! the correct fibonacci stream to decode is
//! ```
//! ddddddddccccccccbbbbbbbbaaaaaaaahhhhhhhhgggggggg....
//! ``` 
use std::{io::{BufWriter, Write, BufReader, Read, SeekFrom, Seek}, fs::File, collections::VecDeque};

use bit_vec::BitVec;
use itertools::{Itertools, izip};
use serde::{Serialize, Deserialize};
use fibonacci_codec::Encode;
use crate::{io::{BusRecord, BusReader, BusHeader, BUS_HEADER_SIZE, DEFAULT_BUF_SIZE, BusWriter}, newpfd::{NewPFDCodec, bitvec_single_fibbonacci_decode}, runlength_codec::RunlengthCodec};

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
pub struct BuszSpecificHeader {
    block_size: u32,
    pfd_block_size: u32,
    lossy_umi: u32
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
        if (setbits_u64(34) as u64) <= block_size_bytes {
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
    let t = s.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join("");
    t
}

pub fn display_u32_in_bits(x: u32) -> String{
    let s: Vec<u32> = (0..32).rev().map (|n| (x >> n) & 1).collect();
    let t = s.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join("");
    t
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


fn compress_barcodes2(records: &[BusRecord]) -> BitVec {
    let runlength_codec = RunlengthCodec {RLE_VAL: 0};

    let mut cb_iter = records.iter().map(|r| r.CB);
    let mut last_el = cb_iter.next().unwrap();
    let mut delta_encoded = Vec::new();
    delta_encoded.push(last_el);
    for el in cb_iter{
        delta_encoded.push(el-last_el);
        last_el=el;
    }
    // println!("Delta enc: {:?}",delta_encoded);
    let runlen_encoded = runlength_codec.encode(delta_encoded.into_iter());
    // println!("run enc: {:?}",runlen_encoded);

    //fibbonaci encoding
    let mut enc =runlen_encoded.fib_encode().unwrap();

    // pad to next multiple of 64
    let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    enc.grow(n_pad, false);
    enc
}

// TODO: the first value encoded is a little funny, it gets incremented 2x!!
// periodic runlength encoding since the UMI resets to smaller values (when changing CB)
// and we cant handle negative differences!
fn compress_umis(records: &[BusRecord]) -> BitVec {
    let runlength_codec = RunlengthCodec {RLE_VAL: 0};
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
    // println!("Delta enc: {:?}",periodic_delta_encoded);
    let runlen_encoded = runlength_codec.encode(periodic_delta_encoded.into_iter());
    // println!("run enc: {:?}",runlen_encoded);
    //fibbonaci encoding
    let mut enc =runlen_encoded.fib_encode().unwrap();

    // pad to next multiple of 64
    let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    enc.grow(n_pad, false);

    enc
}


/// Compress ECs with NewPFD encoding
fn compress_ecs(records: &[BusRecord]) -> BitVec {
    let newpfd = NewPFDCodec::new(PFD_BLOCKSIZE);

    let ecs = records.iter().map(|r|r.EC as u64);
    let mut encoded = newpfd.encode(ecs);

    // pad to next multiple of 64
    let n_pad = round_to_multiple(encoded.len(), 64) - encoded.len();
    encoded.grow(n_pad, false);

    encoded
}

/// Compress counts with RunLength(1) encoding
fn compress_counts(records: &[BusRecord]) -> BitVec {
    let runlength_codec = RunlengthCodec {RLE_VAL: 1};
    let count_iter = records.iter().map(|r| r.COUNT as u64);

    let runlen_encoded = runlength_codec.encode(count_iter);
    // println!("run enc: {:?}",runlen_encoded);

    //fibbonaci encoding
    let mut enc = runlen_encoded.fib_encode().unwrap();

    // pad to next multiple of 64
    let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    enc.grow(n_pad, false);

    enc
}

/// Compress flags with RunLength(0) encoding
fn compress_flags(records: &[BusRecord]) -> BitVec {
    let runlength_codec = RunlengthCodec {RLE_VAL: 0};
    let flag_iter = records.iter().map(|r| r.FLAG as u64);

    let runlen_encoded = runlength_codec.encode(flag_iter);
    // println!("run enc: {:?}",runlen_encoded);

    //fibbonaci encoding
    let mut enc = runlen_encoded.fib_encode().unwrap();

    // pad to next multiple of 64
    let n_pad =  round_to_multiple(enc.len(), 64) - enc.len();
    enc.grow(n_pad, false);

    enc
}

/// puts all given records into a single busz-block, including header
fn compress_busrecords_into_block(records: &[BusRecord]) -> BitVec {
    let mut bcs = dbg!(compress_barcodes2(records));
    let mut umis = dbg!(compress_umis(records));
    let mut ecs = dbg!(compress_ecs(records));
    let mut counts = dbg!(compress_counts(records));
    let mut flags = dbg!(compress_flags(records));

    //concat
    bcs.append(&mut umis);
    bcs.append(&mut ecs);
    bcs.append(&mut counts);
    bcs.append(&mut flags);

    // need to be a multiple of 8
    assert_eq!(bcs.len() % 8 , 0);

    // need to be a multiple of 864
    assert_eq!(bcs.len() % 64 , 0);

    let nbytes = bcs.len() / 8;
    let header = CompressedBlockHeader::new(
        nbytes.try_into().unwrap(), 
        records.len().try_into().unwrap());

    let mut header_bits = bit_vec::BitVec::from_bytes(&header.header_bytes.to_le_bytes()) ; //TODO endian-ness

    header_bits.append(&mut bcs);
    header_bits
    // header.header_bytes
}

pub fn compress_busfile(input: &str, output: &str, blocksize: usize) {

    let reader = BusReader::new(input);
    let mut writer = BufWriter::new(File::open(output).unwrap_or_else(|_| panic!("file not found: {output}")));
    
    // write BusZ header
    
    for chunk in &reader.chunks(blocksize) {
        let records: Vec<BusRecord> = chunk.collect();
        let compressed_block = compress_busrecords_into_block(&records);

        // BitVec.to_bytes() implicitly does bytes in big-endian (first bit in the stream is the high-order bit of the low order byte),
        // yet in the busz file its stored in little endian and we need to convert this
        let little_endian = swap_endian(&compressed_block.to_bytes(), 8);

        writer.write_all(&little_endian).unwrap();
    }
}

/// swaps endianness of the byte-vector
/// assuming 8byte (u64) words
/// simply grabs each 8byte section and reverses the order within
/// # params
/// * wordsize: 4 (bytes) for u32, 8(bytes) for u64
fn swap_endian(bytes: &[u8], wordsize: usize) -> Vec<u8>{
    let mut swapped_endian: Vec<u8> = Vec::with_capacity(bytes.len());
    for bytes in bytes.chunks(wordsize){
        let mut rev_bytes = Vec::new();
        for c in bytes{
            rev_bytes.push(c);
        }
        rev_bytes.reverse();
        swapped_endian.extend(rev_bytes);
    }
    swapped_endian
}

pub fn decompress_busfile(input: &str, output: &str) {

    // let mut reader = BufReader::new(File::open(input).unwrap_or_else(|_| panic!("file not found: {input}")));
    //let mut writer = BufWriter::new(File::open(output).unwrap_or_else(|_| panic!("file not found: {output}")));
    let mut reader = BuszReader::new(input);
    let mut writer = BusWriter::new(
        output,
        reader.bus_header.clone()
    );
    println!("{:?}", reader);

    // read BusZ header
    let mut blockcounter = 0;
    // iterate over the busz-blocks
    loop {
        blockcounter+=1;
        println!("Parsing block {}", blockcounter);

        // get block header
        let mut blockheader_bytes = [0_u8;8];
        reader.reader.read_exact(&mut blockheader_bytes).unwrap();
        if blockheader_bytes == [0,0,0,0,0,0,0,0] {  // EOF
            break;
        }
        let H = CompressedBlockHeader { header_bytes: u64::from_le_bytes(blockheader_bytes)};
        let (block_size_bytes, nrecords) = H.get_blocksize_and_nrecords();
        // println!("H bytes {}, H records {}", block_size_bytes, nrecords);

        // read the busZ-block body
        let mut block_buffer: Vec<u8> = Vec::from_iter((0..block_size_bytes).map(|x| x as u8));
        reader.reader.read_exact(&mut block_buffer).unwrap();

        let bigendian_buf = swap_endian(&block_buffer, 8);
        let the_block = BitVec::from_bytes(&bigendian_buf);
        
        let records = decompress_busz_block( the_block, nrecords.try_into().unwrap());
    
        writer.write_records(&records);
    }
}


// all the encoded parts (cb,umi,ec...) are padded with zeros until a mutiple of 64
// which we need to remove
fn calc_n_trailing_bits(bits_processed: usize) -> usize {
    let padded_size = round_to_multiple(bits_processed, 64);
    let zeros_toremoved = padded_size - bits_processed;
    zeros_toremoved
}

/// Decompress a single BusZ block
fn decompress_busz_block(block_body: BitVec, n_elements: usize) -> Vec<BusRecord>{

    //decode the CBs
    // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
    // println!("Decoding CBs {:?}", block_body);
    let bits_before = block_body.len();
    let mut fibdec = FibbonacciDecoder { bitstream: block_body};

    let CB_RLE_VAL = 0 + 1;  //since everhthing is shifted + 1, the RLE element is also +1
    let mut cb_delta_encoded: Vec<u64> = Vec::with_capacity(n_elements);
    let mut counter = 0;
    while counter < n_elements {
        // println!("Counter {counter}");
        // println!("Bits {:?}", fibdec.bitstream);
        // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
        let item = fibdec.decode_single();
        // println!("Item {item}");
        if item == CB_RLE_VAL {
            let runlength = fibdec.decode_single();
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

    // need to remove the trailing 0 which come from padding to the next mutiple of u64
    // we know that the CB block has a length of n x 64
    let bits_after = fibdec.bitstream.len();
    let bits_processed =  bits_before - bits_after;
    // let padded_size = round_to_multiple(bits_processed, 64);
    let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
    // println!("before {bits_before}, after {bits_after}; zeros to remove: {zeros_toremoved}");

    let remainder = fibdec.bitstream.split_off(zeros_toremoved);
    // println!("truncated: {:?}", fibdec.bitstream);

    assert!(!fibdec.bitstream.any());

    // println!("Bitstream after truncate: {:?}", remainder);

    let mut fibdec = FibbonacciDecoder { bitstream: remainder};

    // =====================================
    // UMIs
    let bits_before = fibdec.bitstream.len();

    let UMI_RLE_VAL = 0 + 1 ; // since shifted
    let mut counter = 0;

    // guarantee that last_barcode != rows[0].barcode
    let mut last_cb = cb_delta_encoded[0] + 1;
    let mut umi =0_u64;
    let mut umis: Vec<u64> = Vec::with_capacity(n_elements);
    while counter < n_elements {
        let diff = fibdec.decode_single() - 1;
        // println!("Decoded item {diff}");
        let current_cb =  cb_delta_encoded[counter];
        if last_cb !=current_cb {
            umi=0;
        }
        // println!("Item {item}");
        if diff == UMI_RLE_VAL - 1 {
            let runlength = fibdec.decode_single();
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
    let bits_after = fibdec.bitstream.len();
    let bits_processed =  bits_before - bits_after;
    let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;

    // println!("before {bits_before}, after {bits_after}; zeros to remove: {zeros_toremoved}");

    let remainder = fibdec.bitstream.split_off(zeros_toremoved);
    // println!("truncated: {:?}", fibdec.bitstream);
    assert!(!fibdec.bitstream.any());

    // println!("Bitstream after truncate: {:?}", remainder);

    // println!("Decoding ECS");
    // here it gets tricky: previously we swapped the entire stream
    // to u64-little endian
    //
    // however, ECs uses u32 as a basis of operatiors
    // and in order to get the bits appear in the right order we need u32-little endian
    // 1. revert the u64 endian swap
    // 2. swap to u32 little endian
    // 3. do decoding
    // 4 swap back to u64-endian
    let _original = swap_endian(&remainder.to_bytes(), 8);
    // println!("Bitstream  original: {:?}", BitVec::from_bytes(&_original));
    let little_endian_32 = swap_endian(&_original, 4);
    // println!("Bitstream  original: {:?}", BitVec::from_bytes(&little_endian_32));

    let remainder_little_endian_32 = BitVec::from_bytes(&little_endian_32);
    drop(remainder);

    let newpfd = NewPFDCodec::new(PFD_BLOCKSIZE);
    let bits_before = remainder_little_endian_32.len();
    let (ecs, re) = newpfd.decode(remainder_little_endian_32, umis.len());
    let bits_after = re.len();
    // println!("before bits {bits_before}, after {bits_after}");
    // println!("ECS {ecs:?}");

    //undo the u64 little endian
    let _original = swap_endian(&re.to_bytes(), 4);
    // println!("Bitstream  original: {:?}", BitVec::from_bytes(&_original));
    let little_endian_64 = swap_endian(&_original, 8);
    // println!("Bitstream  little: {:?}", BitVec::from_bytes(&little_endian_64));

    let remainder_little_endian_64 = BitVec::from_bytes(&little_endian_64);
    
    // ===============================
    // count decoding
    // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
    // println!("Decoding Counts {:?}", remainder_little_endian_64);
    let bits_before: usize = remainder_little_endian_64.len();
    let mut fibdec = FibbonacciDecoder { bitstream: remainder_little_endian_64};

    let COUNT_RLE_VAL = 1;  //since everhthing is shifted + 1, the RLE element is also +1
    let mut counts_encoded: Vec<u64> = Vec::with_capacity(n_elements);
    let mut counter = 0;
    while counter < n_elements {
        // println!("Counter {counter}");
        // println!("Bits {:?}", fibdec.bitstream);
        // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
        let item = fibdec.decode_single();
        // println!("Item {item}");
        if item == COUNT_RLE_VAL {
            let runlength = fibdec.decode_single();
            // let (runlength, remainder) = bitvec_single_fibbonacci_decode(&mut block_body);
            // block_body = remainder;  //weird scope thing: if we call the above with return named block_body, it'll be its own var and not progagate
            // println!("Decoding run of  {runlength}");
            // println!("Body after   {:?}", fibdec.bitstream);
            for _ in 0..runlength{
                counts_encoded.push(COUNT_RLE_VAL);  //due to shift+1
                counter+= 1;
            }
        } else {
            // println!("Decoding sigle element {}", item);
            counts_encoded.push(item);//due to shift+1
            counter+= 1;
        }
    }

    // println!("counts: {:?}", counts_encoded);
    let bits_after = fibdec.bitstream.len();
    let bits_processed =  bits_before - bits_after;
    let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
    // println!("before {bits_before}, after {bits_after}; zeros to remove: {zeros_toremoved}");
    let remainder = fibdec.bitstream.split_off(zeros_toremoved);
    // println!("truncated: {:?}", fibdec.bitstream);
    assert!(!fibdec.bitstream.any());


    // flag decoding
    // println!("Decoding Flags {:?}", remainder);
    let bits_before: usize = remainder.len();
    let mut fibdec = FibbonacciDecoder { bitstream: remainder};

    let FLAG_RLE_VAL = 0+1;  //since everhthing is shifted + 1, the RLE element is also +1
    let mut flag_decoded: Vec<u64> = Vec::with_capacity(n_elements);
    let mut counter = 0;
    while counter < n_elements {
        // println!("Counter {counter}");
        // println!("Bits {:?}", fibdec.bitstream);
        // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
        let item = fibdec.decode_single();
        // println!("Item {item}");
        if item == FLAG_RLE_VAL {
            let runlength = fibdec.decode_single();
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

    let bits_after = fibdec.bitstream.len();
    let bits_processed =  bits_before - bits_after;
    let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
    // println!("before {bits_before}, after {bits_after}; zeros to remove: {zeros_toremoved}");
    let remainder = fibdec.bitstream.split_off(zeros_toremoved);
    // println!("truncated: {:?}", fibdec.bitstream);
    assert!(!fibdec.bitstream.any());
    
    let mut decoded_records = Vec::with_capacity(n_elements);
    for (cb,umi,ec,count,flag) in izip!(cb_delta_encoded, umis, ecs, counts_encoded, flag_decoded) {
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
    // println!("Records: {:?}", decoded_records);
    // println!("reamining stream: {:?}", remainder);

    decoded_records

}

struct FibbonacciDecoder {
    bitstream: BitVec
}

impl FibbonacciDecoder {
    pub fn decode_single(&mut self) -> u64{

        // find the location of the first occurance of 11
        let mut lastbit = false;
        let mut ix: Option<usize> = None;
        for (i, b) in self.bitstream.iter().enumerate() {
            if lastbit & b { // we ran into a 11
                ix = Some(i);
                break;
            }
            lastbit = b
        }
        // split off anything before the delimiter `11`
        if let Some(i) = ix{
            let remainder = self.bitstream.split_off(i + 1); // i marks the last `1`: 00011, but the  
            // silly API: front|11|remainder: it return remainder, stores from in the self.bitstream
            // i.e. we need to swap self.bitstream <-> remainder
            let head = std::mem::replace(&mut self.bitstream, remainder);

            // println!("{x:?}");
            // println!("{:?}", remainder);
    
            let dec: Vec<_> = fibonacci_codec::fib_decode_u64(head.iter())
                .map(|a|a.unwrap())
                .collect();

            // println!("{:?}", dec);
            assert_eq!(dec.len(), 1);
            // println!("single dec {x:?} -> {}", dec[0]);
    
            dec[0]
            // a
        }
        else {
            let error = format!("error in decoding {} {:?}",self.bitstream.len(), self.bitstream);
            panic!("{}", &error);
        }
    
    }
}

#[cfg(test)]
mod test {
    use std::io::Read;
    use bit_vec::BitVec;
    use crate::{busz::{CompressedBlockHeader, compress_barcodes2, compress_umis, compress_ecs, BuszSpecificHeader, swap_endian}, io::{BusRecord, setup_busfile, BusWriter, BusHeader, BusReader}, newpfd::bits64_to_u64};
    use super::{compress_busrecords_into_block, decompress_busz_block, FibbonacciDecoder, BuszReader, decompress_busfile};

    #[test]
    fn test_fibonacci_decoder() {
        let b = BitVec::from_iter(vec![
            true, true, // 1
            false, true, true, // 2
            true, false, true, true, // 4
            true, false, false, false, false, true, true, //14
            ]);
        let mut dc = FibbonacciDecoder {bitstream: b};
        assert_eq!(dc.decode_single(), 1);
        assert_eq!(dc.decode_single(), 2);
        assert_eq!(dc.decode_single(), 4);
        assert_eq!(dc.decode_single(), 14);
    }

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
    fn test_header_encode_decode() {
        let nbytes = 20;
        let nrecords = 10;
        let h = CompressedBlockHeader::new(nbytes, nrecords);

        assert_eq!(h.get_blocksize_and_nrecords().0, nbytes);
        assert_eq!(h.get_blocksize_and_nrecords().1, nrecords);
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

    #[test]
    fn test_encode_decode_block(){
        let v = vec![ 
            BusRecord {CB:0,UMI:10,EC:0,COUNT:1, FLAG: 0 },   // 10
            BusRecord {CB:0,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
            BusRecord {CB:2,UMI:10,EC:1,COUNT:1, FLAG: 0 },   // 0
            BusRecord {CB:2,UMI:11,EC:1,COUNT:1, FLAG: 0 },    // 1
        ];
        let mut comp_block = compress_busrecords_into_block(&v);

        //pop off the header, which is u64
        let mut body = comp_block.split_off(64);

        let header = bits64_to_u64(comp_block);
        let (bytes, nrecords) = CompressedBlockHeader {header_bytes: header}.get_blocksize_and_nrecords();
        assert_eq!(nrecords, v.len() as u64);

        println!("Byts: {bytes}, Elements: {nrecords}, Body: {:?}", body);
        decompress_busz_block( body, nrecords.try_into().unwrap());
    }

    // #[test]
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
    fn test_external2(){
        let mut zreader = BuszReader::new(
            "/tmp/buscompress.busz");
        println!("{:?}", zreader.bus_header);
        println!("{:?}", zreader.busz_header);

        // parse the block header
        let mut blockheader = [0_u8; 8];
        zreader.reader.read_exact(&mut blockheader).unwrap();
        println!("Block header {}, {:?}", blockheader.len(), blockheader);
        let H = CompressedBlockHeader { header_bytes: u64::from_le_bytes(blockheader)};
        let (nbytes, nrecords) = H.get_blocksize_and_nrecords();
        println!("H bytes {}, H records {}", nbytes, nrecords);

        // put the bytes of the block
        let mut buf: Vec<u8> = Vec::from_iter((0..nbytes).map(|x| x as u8));
        zreader.reader.read_exact(&mut buf).unwrap();

        let bigendian_buf = swap_endian(&buf, 8);
        println!("{:?}", buf);
        println!("{:?}", bigendian_buf);

        let b = BitVec::from_bytes(&buf);
        println!("{}, {:?}", b.len(), b);

        let b = BitVec::from_bytes(&bigendian_buf);
        println!("{}, {:?}", b.len(), b);


        let mut fibdec = FibbonacciDecoder { bitstream: b};
        println!("Decode {}", fibdec.decode_single());
        println!("Decode {}", fibdec.decode_single());
        println!("Decode {}", fibdec.decode_single());
        println!("Decode {}", fibdec.decode_single());
        
    }
    #[test]
    fn test3(){

        /*
        the bitstream coming from the busz is in unit of u64, ie. 8bytes
        and in little endian ( least significant byte first), ie a byte array of [1,0,0,0,0,0,0,0] corresponds to the number1
         */

        // order in which the bytes come out of the file
        // let mut b : Vec<bool> = vec![
        //     1,1,1,0,0,0,0,0, // 1
        //     1,1,0,1,0,1,1,1, // 2
        //     0,0,1,0,1,1,0,1  // 3
        //     ].iter().map(|x| *x>0).collect();
        let mut b : Vec<bool> = vec![
            0,0,1,0,1,1,0,1,  // 3
            1,1,0,1,0,1,1,1, // 2
            1,1,1,0,0,0,0,0, // 1
            ].iter().map(|x| *x>0).collect();
            
        // b.reverse();
        println!("{:?}", b);
        let bv = BitVec::from_fn(b.len(), |i| b[i]);
        println!("{:?}", bv);
        let mut fibdec = FibbonacciDecoder { bitstream: bv};
        println!("{}", fibdec.decode_single());
        println!("{}", fibdec.decode_single());
        println!("{}", fibdec.decode_single());
        println!("{}", fibdec.decode_single());
        // println!("{}", fibdec.decode_single());

    }

    #[test]
    fn test_3(){
        let mut zreader = BuszReader::new(
            "/tmp/buscompress.busz");
        println!("{:?}", zreader.bus_header);
        println!("{:?}", zreader.busz_header);

        // parse the block header
        let mut blockheader = [0_u8; 8];
        zreader.reader.read_exact(&mut blockheader).unwrap();
        println!("Block header {}, {:?}", blockheader.len(), blockheader);
        let H = CompressedBlockHeader { header_bytes: u64::from_le_bytes(blockheader)};
        let (nbytes, nrecords) = H.get_blocksize_and_nrecords();
        println!("H bytes {}, H records {}", nbytes, nrecords);

        // put the bytes of the block
        let mut buf: Vec<u8> = Vec::from_iter((0..nbytes).map(|x| x as u8));
        zreader.reader.read_exact(&mut buf).unwrap();

        let bigendian_buf = swap_endian(&buf, 8);
        let b = BitVec::from_bytes(&bigendian_buf);
        let records = decompress_busz_block(b, nrecords as usize);
        assert_eq!(records.len(), 4);

        let mut buf = Vec::new();
        zreader.reader.read_to_end(&mut buf).unwrap();
        println!("{:?}", buf)
    }

    #[test]
    fn test_4(){
        decompress_busfile("/tmp/buscompress.busz", "/tmp/buscompress_lala.bus");
        let r = BusReader::new("/tmp/buscompress_lala.bus");
        let records:Vec<_> = r.collect();
        assert_eq!(records.len(), 4);
        println!("{:?}", records);
    }

    #[test]
    fn test_5(){
        decompress_busfile(
            "/home/michi/bus_testing/bus_output_shortest/output.corrected.sort.busz",
            "/tmp/buscompress_lala.bus");
        let r = BusReader::new("/tmp/buscompress_lala.bus");
        let records:Vec<_> = r.collect();

        let r_original = BusReader::new("/home/michi/bus_testing/bus_output_shortest/output.corrected.sort.bus");
        let records_original:Vec<_> = r_original.collect();

        assert_eq!(records.len(), records_original.len());
        assert_eq!(records, records_original);

        // assert_eq!(
        //     records.iter().map(|r|r.CB).collect::<Vec<_>>(),
        //     records_original.iter().map(|r|r.CB).collect::<Vec<_>>(),
        // );
        // assert_eq!(
        //     records.iter().map(|r|r.UMI).collect::<Vec<_>>(),
        //     records_original.iter().map(|r|r.UMI).collect::<Vec<_>>(),
        // );
        // assert_eq!(
        //     records.iter().map(|r|r.COUNT).collect::<Vec<_>>(),
        //     records_original.iter().map(|r|r.COUNT).collect::<Vec<_>>(),
        // );
        // assert_eq!(
        //     records.iter().map(|r|r.FLAG).collect::<Vec<_>>(),
        //     records_original.iter().map(|r|r.FLAG).collect::<Vec<_>>(),
        // );
        // assert_eq!(
        //     records.iter().map(|r|r.EC).collect::<Vec<_>>(),
        //     records_original.iter().map(|r|r.EC).collect::<Vec<_>>(),
        // );

        // println!("{:?}", records);
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


const BUSZ_HEADER_SIZE: usize = 4+4+4;
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
pub struct BuszHeader {
    block_size: u32,
    pfd_block_size: u32,
    lossy_umi: u32,
}
impl BuszHeader {
    /// desearializes a BusHeader from Bytes; when reading busfiles
    pub fn from_bytes(bytes: &[u8]) -> BuszHeader {
        let header_struct: BuszHeader =
            bincode::deserialize(bytes).expect("FAILED to deserialze busz header");
        assert_eq!(
            header_struct.lossy_umi, 0,
            "lossy_umi != 0 not supported"
        );
        header_struct
    }
    /// seialize the header to bytes
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
        let mut file_handle = File::open(filename).expect("FAIL");

        // advance the file point 20 bytes (pure header) + header.tlen (the variable part)
        let to_seek = BUS_HEADER_SIZE as u64 + bus_header.get_tlen() as u64;
        let _x = file_handle.seek(SeekFrom::Start(to_seek)).unwrap();

        // now we're at the start of the BusZHeader
        let mut buszheader_bytes = [0_u8; BUSZ_HEADER_SIZE];
        file_handle.read_exact(&mut buszheader_bytes).unwrap();
        println!("{:?}", buszheader_bytes);

        let bzheader = BuszHeader::from_bytes(&buszheader_bytes);
    
        //FH is now at the start of the actual contents
        let buf = BufReader::with_capacity(bufsize, file_handle);

        let buffer = VecDeque::with_capacity(bzheader.block_size as usize);
        BuszReader { bus_header, busz_header:bzheader, reader: buf, buffer:buffer }
    }

    pub fn load_busz_block(&mut self) -> Option<Vec<BusRecord>>{

        // get block header
        let mut blockheader_bytes = [0_u8;8];
        self.reader.read_exact(&mut blockheader_bytes).unwrap();
        if blockheader_bytes == [0,0,0,0,0,0,0,0] {  // EOF
            // break;
            return None
        }
        let H = CompressedBlockHeader { header_bytes: u64::from_le_bytes(blockheader_bytes)};
        let (block_size_bytes, nrecords) = H.get_blocksize_and_nrecords();
        // println!("H bytes {}, H records {}", block_size_bytes, nrecords);

        // read the busZ-block body
        let mut block_buffer: Vec<u8> = Vec::from_iter((0..block_size_bytes).map(|x| x as u8));
        self.reader.read_exact(&mut block_buffer).unwrap();

        let bigendian_buf = swap_endian(&block_buffer, 8);
        let the_block = BitVec::from_bytes(&bigendian_buf);
        
        let records = decompress_busz_block( the_block, nrecords.try_into().unwrap());
    
        Some(records)
    }
}

impl Iterator for BuszReader {
    type Item=BusRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer.is_empty(){
            println!("buffer empty, loading new block");

            let block =self.load_busz_block(); 
            // pull in another block
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
    // fn next_old(&mut self) -> Option<Self::Item> {
    //     let mut buffer = [0; BUS_ENTRY_SIZE]; // TODO this gets allocated at each next(). We could move it out into the struct: Actually doesnt make a diference, bottleneck is the reading
    //     match self.reader.read(&mut buffer) {
    //         Ok(BUS_ENTRY_SIZE) => Some(BusRecord::from_bytes(&buffer)),
    //         Ok(0) => None,
    //         Ok(n) => {
    //             let s: BusRecord = BusRecord::from_bytes(&buffer);
    //             panic!(
    //                 "Wrong number of bytes {:?}. Buffer: {:?} record: {:?}",
    //                 n, buffer, s
    //             )
    //         }
    //         Err(e) => panic!("{:?}", e),
    //     }
    // }
}


#[test]
fn bit_vec_testing() {

    let mut b = BitVec::from_elem(33, false);
    b.set(15, true);

    println!("{:?}",b);

    let x = b.storage();
    println!("{:?}",x)
}

#[test]
fn bitvec_testing() {

    use bitvec::prelude::*;

    // fixed size array
    let arr = bitarr![u32, Lsb0; 0; 80];

    // a slice
    let bits = bits![u16, Msb0; 0; 40];

    let mut data = [0u8; 2];
    let bits = data.view_bits_mut::<Lsb0>();

    bits.set(9, true);
    // bits.set(9, true);
    // // Unsigned integers (scalar, array,
    // // and slice) can be borrowed.
    // let data = 0x2021u16;
    // let bits = data.view_bits::<Msb0>();
    // let data = [0xA5u8, 0x3C];
    // let bits = data.view_bits::<Lsb0>();

    // // Bit-slices can split anywhere.
    // let (head, rest) = bits.split_at(4);


    let x = BitVec::from_bitslice(bits);
    println!("{:?}",x);
    let y = &x[..10];
    println!("{:?}",y);


}