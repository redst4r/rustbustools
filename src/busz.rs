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
use bitvec::{slice::BitSlice, prelude::Msb0, field::BitField};
use bitvec::prelude as bv;
use itertools::{Itertools, izip};
use serde::{Serialize, Deserialize};
use fibonacci_codec::Encode;
use crate::{io::{BusRecord, BusReader, BusHeader, BUS_HEADER_SIZE, DEFAULT_BUF_SIZE, BusWriter}, newpfd::{NewPFDCodec, bitvec_single_fibbonacci_decode}, runlength_codec::RunlengthCodec, myfibonacci::{self, bitslice_to_fibonacci}, newpfd_bitvec};

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

#[deprecated(since="0.7.0", note="please use `decompress_busfile` instead")]
pub fn decompress_busfile_old(input: &str, output: &str) {

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

/// Decompress a single BusZ block
fn decompress_busz_block(block_body: BitVec, n_elements: usize) -> Vec<BusRecord>{

    //decode the CBs
    // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
    // println!("Decoding CBs\n{:?}", block_body);
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

    let after_cb_pos = bits_processed+zeros_toremoved;
    // println!("CBs: bits proc {}, bits padded {}", bits_processed, zeros_toremoved);

    // println!("CBs buffer pos: {}", after_cb_pos);
    // println!("Bitstream after truncate: {:?}", remainder);

    // println!("Decoding UMIs\n{:?}", remainder);

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
    let after_umi_pos = after_cb_pos + bits_processed+zeros_toremoved;
    // println!("UMIs: bits proc {}, bits padded {}", bits_processed, zeros_toremoved);
    
    // println!("UMIs buffer pos: {}", after_umi_pos);



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
    // println!("Decoding ECs\n{:?}", remainder);


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

    let after_ec_pos = after_umi_pos + bits_before - bits_after;
    // println!("ECs buffer pos: {}", after_ec_pos);
    // println!("before bits {bits_before}, after {bits_after}");
    // println!("ECS {ecs:?}");


    // println!("Decoding Counts before backtransform\n{:?}", re);

    //undo the u64 little endian
    let _original = swap_endian(&re.to_bytes(), 4);
    // println!("Decoding Counts original\n{:?}", BitVec::from_bytes(&_original));

    let little_endian_64 = swap_endian(&_original, 8);
    // println!("Bitstream  little:\n{:?}", BitVec::from_bytes(&little_endian_64));

    let remainder_little_endian_64 = BitVec::from_bytes(&little_endian_64);
    
    // ===============================
    // count decoding
    // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
    // println!("Decoding Counts {:?}", remainder_little_endian_64);
    let bits_before: usize = remainder_little_endian_64.len();

    // println!("Decoding Counts after backtransform\n{:?}", remainder_little_endian_64);


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

    let after_count_pos = after_ec_pos + bits_processed+zeros_toremoved;
    // println!("Counts buffer pos: {}", after_count_pos);

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
    
    let after_flag_pos = after_count_pos + bits_processed+zeros_toremoved;
    // println!("Flags buffer pos: {}", after_flag_pos);

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
    use crate::{busz::{CompressedBlockHeader, compress_barcodes2, compress_umis, compress_ecs, BuszSpecificHeader, swap_endian, decompress_busfile}, io::{BusRecord, setup_busfile, BusWriter, BusHeader, BusReader}, newpfd::bits64_to_u64};
    use super::{compress_busrecords_into_block, decompress_busz_block, FibbonacciDecoder, BuszReader, decompress_busfile_old};

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
        decompress_busfile_old("/tmp/buscompress.busz", "/tmp/buscompress_lala.bus");
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

    /// takes the next 8 bytes (u64) out of the stream and interprets it as a busheader
    /// if we encounter the zero byte [00000000], this indicates the EOF
    fn load_busz_header(&mut self) -> Option<CompressedBlockHeader>{
        // get block header
        let mut blockheader_bytes = [0_u8;8];
        self.reader.read_exact(&mut blockheader_bytes).unwrap();
        if blockheader_bytes == [0,0,0,0,0,0,0,0] {  // EOF
            // break;
            return None
        }
        let H = CompressedBlockHeader { header_bytes: u64::from_le_bytes(blockheader_bytes)};
        // println!("H bytes {}, H records {}", H.get_blocksize_and_nrecords().0, H.get_blocksize_and_nrecords().1);
        Some(H)
    }

    pub fn load_busz_block(&mut self) -> Option<Vec<BusRecord>>{

        let h = self.load_busz_header();
        if h.is_none(){
            return None
        }
        let header = h.unwrap();
        let (block_size_bytes, nrecords) = header.get_blocksize_and_nrecords();

        // read the busZ-block body
        let mut block_buffer: Vec<u8> = Vec::from_iter((0..block_size_bytes).map(|x| x as u8));
        self.reader.read_exact(&mut block_buffer).unwrap();

        let bigendian_buf = swap_endian(&block_buffer, 8);
        let the_block = BitVec::from_bytes(&bigendian_buf);
        
        let records = decompress_busz_block( the_block, nrecords.try_into().unwrap());
    
        Some(records)
    }

    pub fn load_busz_block_faster(&mut self) -> Option<Vec<BusRecord>>{

        let h: Option<CompressedBlockHeader> = self.load_busz_header();
        if h.is_none(){
            return None
        }
        let header = h.unwrap();
        let (block_size_bytes, nrecords) = header.get_blocksize_and_nrecords();

        // read the busZ-block body
        let mut block_buffer: Vec<u8> = Vec::from_iter((0..block_size_bytes).map(|x| x as u8));
        self.reader.read_exact(&mut block_buffer).unwrap();

        // conversion to big-endian to make the reading work right
        let bigendian_buf = swap_endian(&block_buffer, 8);
        let theblock: &bv::BitVec<u8, bv::Msb0> = &bv::BitVec::from_slice(&bigendian_buf);
        // println!("the_block: {}, {:?}", theblock.len(), theblock);
        let slice: &BitSlice<u8, bv::Msb0> = &theblock.as_bitslice();


        let mut block = BuszBlock::new(&slice,nrecords as usize);
        
        let records =block.parse_block();
    
        Some(records)
    }
}

impl Iterator for BuszReader {
    type Item=BusRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer.is_empty(){
            println!("buffer empty, loading new block");

            // let block =self.load_busz_block(); 
            let block =self.load_busz_block_faster(); 
// 
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
}


/// to keep track in which section we are in the busz-block
#[derive(Debug, PartialEq, Eq)]
enum BuszBlockState {
    Header,
    CB,
    UMI,
    EC,
    COUNT,
    FLAG,
    FINISHED
}

/// A single block of a busz-file, header has already been parsed
/// Contains all bytes from the block (we know how many bytes from the header)
/// and parses those bytes into BusRecords
#[derive(Debug)]
struct BuszBlock <'a> {
    buffer: &'a BitSlice<u8, Msb0>,
    pos: usize, // where we are currently in the buffer
    n_elements: usize ,// how many busrecords are stored in the block
    state: BuszBlockState
}


impl <'a> BuszBlock <'a>{
    // pub fn new(buffer: &'a BitSlice<u8, Msb0>, n_elements: usize) -> Self {
    pub fn new(buffer: &'a BitSlice<u8, bv::Msb0>, n_elements: usize) -> Self {
        // TODO: warning, buffer must be conveted in a special way if the bytes come out of a file
        // see BuszReader::load_busz_block_faster
        // cant do it in herer due to lifetime issues
        BuszBlock { 
            buffer: buffer, 
            pos: 0, 
            n_elements: n_elements, 
            state: BuszBlockState::CB 
        }
    }

    fn parse_cb(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::CB);
        assert_eq!(self.pos, 0, "must be at beginning of buffer to parse CBs");

        //decode the CBs
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
        // println!("Decoding CBs {:?}", block_body);

        let cb_buffer = &self.buffer[self.pos..];
        // let bits_before = self.buffer.len();

        let mut fibdec = myfibonacci::MyFibDecoder::new(cb_buffer);
        // let mut fibdec = FibbonacciDecoder { bitstream: block_body};

        let CB_RLE_VAL = 0 + 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut cb_delta_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            // println!("Counter {counter}");
            // println!("Bits {:?}", fibdec.bitstream);
            // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
            let item = bitslice_to_fibonacci( fibdec.next().unwrap());
            // println!("Item {item}");
            if item == CB_RLE_VAL {
                let runlength = bitslice_to_fibonacci( fibdec.next().unwrap());
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
        self.state = BuszBlockState::UMI;

        cb_delta_encoded
    }

    fn parse_umi(&mut self, cbs: &[u64]) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::UMI);

        let umi_buffer = &self.buffer[self.pos..];
        let mut fibdec = myfibonacci::MyFibDecoder::new(umi_buffer);

        // =====================================
        // UMIs

        let UMI_RLE_VAL = 0 + 1 ; // since shifted
        let mut counter = 0;

        // guarantee that last_barcode != rows[0].barcode
        let mut last_cb = cbs[0] + 1;
        let mut umi =0_u64;
        let mut umis: Vec<u64> = Vec::with_capacity(self.n_elements);
        while counter < self.n_elements {
            let diff = bitslice_to_fibonacci( fibdec.next().unwrap()) - 1;
            // println!("Decoded item {diff}");
            let current_cb =  cbs[counter];
            if last_cb !=current_cb {
                umi=0;
            }
            // println!("Item {item}");
            if diff == UMI_RLE_VAL - 1 {
                let runlength = bitslice_to_fibonacci( fibdec.next().unwrap());
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
        self.state = BuszBlockState::EC;

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
       
        // let remainder_little_endian_32: bv::BitVec<u8, Msb0> =  bv::BitVec::from_vec(little_endian_32_bytes);
        let remainder_little_endian_32: &bv::BitSlice<u8, Msb0> =  bv::BitSlice::from_slice(&little_endian_32_bytes);

        // println!("len remainder_little_endian_32 {}", remainder_little_endian_32.len());
        // println!("remainder_little_endian_32\n{}", bitstream_to_string(remainder_little_endian_32));
        // println!("main buf: remainder_little_endian_32\n{}", bitstream_to_string(ec_buffer));

        let (ecs, bits_consumed) = newpfd_bitvec::decode(remainder_little_endian_32, self.n_elements, PFD_BLOCKSIZE);

        // println!("after EC: {}", bitstream_to_string(&remainder_little_endian_32[bits_consumed..]));

        self.pos = self.pos + bits_consumed;
        // println!("after EC main buffer: {}", bitstream_to_string(&self.buffer[self.pos..]));
        // println!("after EC main buffer: {}", bitstream_to_string(&ec_buffer[bits_consumed..]));

        // println!("ECs: {}", ecs.len());
        // println!("ECs: bits proc {}, bits padded 0", bits_consumed);
        // println!("ECs: buffer at {}", self.pos);


        self.state = BuszBlockState::COUNT;

        ecs
    }

    fn parse_counts(&mut self) ->Vec<u64> {
        // ===============================
        // count decoding
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
        // println!("Decoding Counts {:?}", remainder_little_endian_64);
        assert_eq!(self.state, BuszBlockState::COUNT);

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
        let count_buffer: &bv::BitSlice<u8, Msb0> =  bv::BitSlice::from_slice(&d);

        // let mut count_buffer = &self.buffer[self.pos..];
        // println!("before transform count_buffer\n{}", bitstream_to_string(count_buffer));

        let mut fibdec = myfibonacci::MyFibDecoder::new(count_buffer);

        let COUNT_RLE_VAL = 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut counts_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            let item = bitslice_to_fibonacci( fibdec.next().unwrap());
            // println!("Item {item}");
            if item == COUNT_RLE_VAL {
                let runlength = bitslice_to_fibonacci( fibdec.next().unwrap());
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
        self.state = BuszBlockState::FLAG;

        counts_encoded
    }

    fn parse_flags(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::FLAG);

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
        let flag_buffer: &bv::BitSlice<u8, Msb0> =  bv::BitSlice::from_slice(&d);
        // println!("flag_buffer\n{}", bitstream_to_string(flag_buffer));


        let mut fibdec = myfibonacci::MyFibDecoder::new(flag_buffer);

        let FLAG_RLE_VAL = 0+1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut flag_decoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            // println!("Counter {counter}");
            // println!("Bits {:?}", fibdec.bitstream);
            // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
            let item = bitslice_to_fibonacci( fibdec.next().unwrap());
            // println!("Item {item}");
            if item == FLAG_RLE_VAL {
                let runlength = bitslice_to_fibonacci( fibdec.next().unwrap());
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
        self.state = BuszBlockState::FINISHED;
        
        assert_eq!(self.pos, self.buffer.len(), "still leftover bits in the buffer!"); 

        flag_decoded
    }

    pub fn parse_block(&mut self) -> Vec<BusRecord>{
        // println!("Cb: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let cbs = self.parse_cb();
        // println!("CBs: {:?}", cbs);


        // println!("umi: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let umis = self.parse_umi(&cbs);

        // println!("ec: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let ec = self.parse_ec();

        // println!("count: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let count = self.parse_counts();

        // println!("flag: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let flag = self.parse_flags();
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
fn bitslice_to_bytes(bits: &BitSlice<u8, Msb0>) -> Vec<u8>{

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

#[test]
fn test_busblock(){
    use bitvec::prelude as bv;

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
    // println!("H bytes {}, H records {}", nbytes, nrecords);

    // put the bytes of the block
    let mut buf: Vec<u8> = Vec::from_iter((0..nbytes).map(|x| x as u8));
    zreader.reader.read_exact(&mut buf).unwrap();

    // conversion to big-endian to make the reading work right
    let bigendian_buf = swap_endian(&buf, 8);
    let theblock: &bv::BitVec<u8, bv::Msb0> = &bv::BitVec::from_slice(&bigendian_buf);
    // println!("the_block: {}, {:?}", theblock.len(), theblock);
    let slice: &BitSlice<u8, bv::Msb0> = &theblock.as_bitslice();


    let mut block = BuszBlock::new(slice, nrecords as usize);

    
    // // println!("Cb: {}", bitstream_to_string(&block.buffer[block.pos..]));
    // let cbs = block.parse_cb();

    // // println!("umi: {}", bitstream_to_string(&block.buffer[block.pos..]));
    // let umis = block.parse_umi(&cbs);

    // // println!("ec: {}", bitstream_to_string(&block.buffer[block.pos..]));
    // let ec = block.parse_ec();

    // // println!("count: {}", bitstream_to_string(&block.buffer[block.pos..]));
    // let count = block.parse_counts();

    // // println!("flag: {}", bitstream_to_string(&block.buffer[block.pos..]));
    // let flag = block.parse_flags();
    // println!("{:?}", cbs);
    // println!("{:?}", umis);
    // println!("{:?}", ec);
    // println!("{:?}", count);
    // println!("{:?}", flag);

    // let mut decoded_records = Vec::with_capacity(nrecords.try_into().unwrap());
    // for (cb,umi,ec,count,flag) in izip!(cbs, umis, ec, count, flag) {
    //     decoded_records.push(
    //         BusRecord {
    //             CB: cb, 
    //             UMI: umi, 
    //             EC:ec.try_into().unwrap(), 
    //             COUNT:count.try_into().unwrap(), 
    //             FLAG:flag.try_into().unwrap()
    //         }
    //     );
    // }
    let records = block.parse_block();
    println!("Records: {:?}", records);
    
}

fn bitstream_to_string(buffer: &bv::BitSlice<u8, Msb0>) -> String{
    let s = buffer.iter().map(|x| if *x==true{"1"} else {"0"}).join("");
    s
}

#[test]
fn endian_Swapping() {
    let v = vec![0_u8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
    let a = swap_endian(&v, 8);
    let b = swap_endian(&a, 4);
    let c = swap_endian(&b, 4);
    let d = swap_endian(&c, 8);
    assert_eq!(a,d);
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

    let mut data = [0u8; 16];
    let mut data = [0x00u8,0x11u8,0x22u8,0x33u8,0x44,0x55,0x66,0x77,0x88,0x99,0xAAu8, 0xBB, 0xCC, 0xDD, 0xEE, 0xFF];

    let mut bits = data.view_bits_mut::<Lsb0>();


    bits.set(0, true);
    bits.set(1, true);
    bits.set(4, true);
    bits.set(100, true);


    let mut bytes = Vec::new();
    for i in 0..16 {
        let b = &bits[i*8.. (i+1)*8];
        let a: u8 = b.load_be();
        bytes.push(a);
    }
    println!("{:?}", bytes);

    let s = swap_endian(&bytes, 8);
    println!("{:?}", s);

    let t = swap_endian(&s, 4);
    println!("{:?}", t);

    
    // let a_be: u32 = x.load_be();
    // let a_le: u32 = x.load_le();
    // println!("Bits {:?} Int_be: {}, Int_le {}",x, a_be, a_le);
    // u64::from_be_bytes(bytes[..8].try_into().unwrap());

    // x.
}