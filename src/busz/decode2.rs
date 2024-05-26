use std::{collections::VecDeque, fs::File, io::{BufReader, Read}};
use crate::{busz::{encode::compress_busrecords_into_block, utils::{bitstream_to_string, calc_n_trailing_bits, swap_endian}, BUSZ_HEADER_SIZE}, io::{BusHeader, BusParams, BusRecord, BusWriterPlain, CUGIterator, BUS_HEADER_SIZE, DEFAULT_BUF_SIZE}};
use bitvec::prelude as bv;
use itertools::izip;
use fastfibonacci::{byte_decode::{faster::{FB_LOOKUP_NEW_U16, FB_LOOKUP_NEW_U8}, u64_fibdecoder::U64Decoder}, FbDec};
use super::{BuszBitSlice, BuszHeader, CompressedBlockHeader, PFD_BLOCKSIZE};

pub struct BuszReader <'a> {
    params: BusParams,
    busz_header: BuszHeader,
    reader:  Box<dyn Read+ 'a>, //ugly way to store any Read-like object in here, BufferedReader, File, Cursor, or just a vec<u8>!
    buffer: VecDeque<BusRecord>,
    is_done: bool,  // weird case where .groupby() keeps calling .next() even after it recieved None. We need to keep emitting None once we're done with the file
}

impl <'a>BuszReader<'a> {
    /// main constructor for busreader, buffersize is set to best performance
    pub fn new(filename: &str) -> Self {
        let file_handle = File::open(filename).expect("FAIL");       
        let buf = BufReader::with_capacity(DEFAULT_BUF_SIZE, file_handle);
        Self::from_read(buf)
    }

    /// Construct a decoder from a `Read`-able object. Make sure that the supplied
    /// reader is buffered, otherwise each call will result in an io/operation.
    /// (proably not super-critical, as BuszReader reads entire blocks of bytes anyways, i.e. its internally buffered)
    pub fn from_read(mut reader: impl Read + 'a) -> Self {
        // parse header
        let mut header_bytes = [0_u8; BUS_HEADER_SIZE];
        reader.read_exact(&mut header_bytes).expect("failed to read header");
        let header = BusHeader::from_bytes(&header_bytes);
        let params = BusParams { cb_len: header.cb_len, umi_len: header.umi_len };
        
        assert_eq!(
            &header.magic, b"BUS\x01",
            "Header struct not matching; MAGIC is wrong"
        );

        // the variable header
        let mut var_buffer = Vec::with_capacity(header.tlen as usize);
        for _i in 0..header.tlen {
            var_buffer.push(0_u8);
        }
        reader.read_exact(&mut var_buffer).expect("failed to read variable header");
        
        // BusZHeader
        let mut buszheader_bytes = [0_u8; BUSZ_HEADER_SIZE];
        reader.read_exact(&mut buszheader_bytes).unwrap();
        let busz_header = BuszHeader::from_bytes(&buszheader_bytes);
        
        // done, create the BuszReader
        let buffer = VecDeque::with_capacity(busz_header.block_size as usize);
        BuszReader { params, busz_header, reader: Box::new(reader), buffer, is_done: false }
    }

    /// takes the next 8 bytes (u64) out of the stream and interprets it as a buszheader
    /// if we encounter the zero byte [00000000], this indicates the EOF
    fn load_block_header(&mut self) -> Option<CompressedBlockHeader>{
        // get block header
        let mut blockheader_bytes = [0_u8;8];
        self.reader.read_exact(&mut blockheader_bytes).expect("couldnt read the busz block header");
        if blockheader_bytes == [0,0,0,0,0,0,0,0] {  // EOF
            return None
        }
        let h = CompressedBlockHeader { header_bytes: u64::from_le_bytes(blockheader_bytes)};
        // println!("H bytes {}, H records {}", h.get_blocksize_and_nrecords().0, h.get_blocksize_and_nrecords().1);
        Some(h)
    }

    /// loads a BusZ-block from the stream.
    /// We parse the header, which tells us how many bytes and records the block has
    /// with that we pull that amount of bytes out of the stream and parse them into Records.
    fn load_busz_block_faster(&mut self) -> Option<Vec<BusRecord>>{

        let h: Option<CompressedBlockHeader> = self.load_block_header();
        if h.is_none(){
            return None
        }
        let header = h.unwrap();
        let (block_size_bytes, nrecords) = header.get_blocksize_and_nrecords();

        // println!("BusZ block-header bytes:{block_size_bytes} #records{nrecords}");
        // theres usually ~50k bytes in a block, 10000 records
        // read the busZ-block body
        let mut block_buffer: Vec<u8> = vec![0;block_size_bytes as usize];
        self.reader.read_exact(&mut block_buffer).unwrap();

        
        // HERES THE NEW STUFF
        // create a block from the bytes, parse it
        let mut block = BuszBlock::new(block_buffer.as_slice(), nrecords as usize);
        
        let records =block.parse_block();
    
        Some(records)
    }

    pub fn get_params(&self) -> &BusParams {
        &self.params
    }

    // pub fn get_busz_header(&self) -> BuszHeader{
    //     self.busz_header.clone()
    // }
}

impl <'a>Iterator for BuszReader<'a> {
    type Item=BusRecord;

    fn next(&mut self) -> Option<Self::Item> {

        if self.is_done {  // weird little workaround such that the iterator keeps emitting None when were done without trying to read more blocks from disk
            return None
        }

        if self.buffer.is_empty(){
            // pull in another block
            // println!("Pulling in another block");
            let block =self.load_busz_block_faster(); 
            match block {
                Some(records) => {
                    // println!("Pulled in another {} records", records.len());
                    self.buffer = records.into()
                },
                None => {
                    // println!("no more records to pull");
                    self.is_done = true;
                    return None
                }
            };
        }
        match self.buffer.pop_front() {
            Some(record) => {
                // println!("Emitting {:?}\nRemainingBuffer: {:?}\n", record, self.buffer);
                Some(record)
            },
            None => panic!("cant happen")
        }
    }
}

// To make `grouby_cb()` etc work
impl <'a> CUGIterator for BuszReader<'a> {}

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
    buffer: &'a [u8],
    pos: usize, // where we are currently in the buffer
    n_elements: usize ,// how many busrecords are stored in the block
    state: BuszBlockState,
}

impl <'a> BuszBlock <'a> {
    fn new(buffer: &'a [u8], n_elements: usize) -> Self {
        // TODO: warning, buffer must be conveted in a special way if the bytes come out of a file
        // see BuszReader::load_busz_block_faster
        // cant do it in herer due to lifetime issues
        BuszBlock { 
            buffer, 
            pos: 0, 
            n_elements, 
            state: BuszBlockState::Cb,
        }
    }

    // fn fibonacci_factory<R:Read>(stream: R) -> U64Decoder<R>{
    //      U64Decoder::new(stream)
    // }

    fn fibonacci_factory<R:Read+'a>(stream: R) -> fastfibonacci::byte_decode::faster::FastFibonacciDecoderNewU16<'a, R>{
        fastfibonacci::byte_decode::faster::FastFibonacciDecoderNewU16::new(stream, &FB_LOOKUP_NEW_U16, false, fastfibonacci::byte_decode::faster::StreamType::U64)
    }

    fn parse_cb(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Cb);
        assert_eq!(self.pos, 0, "must be at beginning of buffer to parse CBs");

        // The CBs are in a u64 block, hence use a decoder for u64s
        let mut dd = Self::fibonacci_factory(&self.buffer[self.pos..]);
        // let mut dd = U64Decoder::new(&self.buffer[self.pos..]);

        //decode the CBs
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb

        const CB_RLE_VAL:u64 = 0 + 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut cb_delta_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        while cb_delta_encoded.len() < self.n_elements {
            let item = dd.next().unwrap();
        
            if item == CB_RLE_VAL {
                let runlength = dd.next().unwrap() as usize;
                // add the same element over and over
                cb_delta_encoded.resize(cb_delta_encoded.len() + runlength, CB_RLE_VAL -1);
            } else {
                cb_delta_encoded.push(item - 1);//due to shift+1
            }
        }
        // undo delta encoding, i,e. sumsum
        cb_delta_encoded.iter_mut().fold(0, |acc, x| {
            *x += acc;
            *x
        });

        // we're done decoding the CBs;
        // lets make sure we're flush with the u64 array
        // ie whatever is left in the current u64 is just zero-bits
        assert!(dd.is_clean());

        // update the position in the bytes array
        // we've read 8bytes per u64!
        let bytes_consumed = dd.get_consumed_bytes();
        self.pos += bytes_consumed;

        self.state = BuszBlockState::Umi;

        cb_delta_encoded
    }
    
    fn parse_umi(&mut self, cbs: &[u64]) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Umi);

        let mut dd = Self::fibonacci_factory(&self.buffer[self.pos..]);
        // let mut dd = U64Decoder::new(&self.buffer[self.pos..]);

        // let umi_buffer = &self.buffer[self.pos..];
        // let mut fibdec = Self::fibonacci_factory(umi_buffer);

        const UMI_RLE_VAL: u64 = 0 + 1 ; // since shifted

        // guarantee that last_barcode != rows[0].barcode
        let mut last_cb = cbs[0] + 1;
        let mut umi =0_u64;
        let mut umis: Vec<u64> = Vec::with_capacity(self.n_elements);
        while umis.len() < self.n_elements {
            let diff = dd.next().unwrap() - 1;

            let current_cb =  cbs[umis.len()];
            if last_cb !=current_cb {
                umi=0;
            }

            if diff == UMI_RLE_VAL - 1 {
                let runlength = dd.next().unwrap() as usize; 
                // adding the same element repeatedly
                umis.resize(umis.len() + runlength, umi - 1);
            } else {
                umi+= diff;
                // decoding single element
                umis.push(umi - 1);//due to shift+1
            }
            last_cb = current_cb;
        }

        // we're done decoding the UMIs;
        // lets make sure we're flush with the u64 array
        // ie whatever is left in the current u64 is just zero-bits
        assert!(dd.is_clean());

        // update the position in the bytes array
        // we've read 8bytes per u64!
        let bytes_consumed = dd.get_consumed_bytes();
        self.pos += bytes_consumed;

        self.state = BuszBlockState::Ec;

        umis
    }

    fn parse_ec(&mut self) -> Vec<u64> {


        let bytes = &self.buffer[self.pos..];       

        // annoyingly we have to swap the endianness for BitVec to work properly
        // currently the first bits to consider are stored in byte[3]
        // Bitstream: 
        // 00000000000000000000000001111110|00000000000000000000000001011100|
        // instead of 
        // 01111110000000000000000000000000|01011100000000000000000000000000

        // this swapping copies, i.e. it DOES NOT affect the self.buffer!
        let little_endian_32_bytes = swap_endian(&bytes, 4);  //TODO performance: could be done in-place
        let remainder_little_endian_32: &BuszBitSlice =  bv::BitSlice::from_slice(&little_endian_32_bytes);
        
        // IMPORTANT NOTE: the next line copies self.buffer
        // which is subsequently swapped!
        // if WE DONT COPY, the entire self.buffer would be swapped
        // let mut bytes = bitslice_to_bytes(ec_buffer);       
        // swap_endian8_swap_endian4_inplace(&mut bytes);
        // let remainder_little_endian_32: &BuszBitSlice =  bv::BitSlice::from_slice(&bytes);

        // println!("Decoding EC: with {} els", self.n_elements);
        // println!("Bitstream: \n{}", bitstream_to_string_pretty(remainder_little_endian_32, 32));
        let (ecs, bits_consumed) = newpfd::newpfd_bitvec::decode(remainder_little_endian_32, self.n_elements, PFD_BLOCKSIZE);

        assert_eq!(bits_consumed % 32, 0);
        let bytes_consumed = bits_consumed / 8;

        self.pos += bytes_consumed;
        self.state = BuszBlockState::Count;

        ecs
    }

    fn parse_counts(&mut self) -> Vec<u64> {
        // ===============================
        // count decoding
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
        // println!("Decoding Counts {:?}", remainder_little_endian_64);
        assert_eq!(self.state, BuszBlockState::Count);

        let mut dd = Self::fibonacci_factory(&self.buffer[self.pos..]);
        // let mut dd = U64Decoder::new(&self.buffer[self.pos..]);

        const COUNT_RLE_VAL: u64 = 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut counts_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        while counts_encoded.len() < self.n_elements {
            let item = dd.next().unwrap();

            if item == COUNT_RLE_VAL {
                let runlength = dd.next().unwrap() as usize;
                // just append the RLE repeatedly:
                counts_encoded.resize(counts_encoded.len() + runlength, COUNT_RLE_VAL);
            } else {
                //decode single element
                counts_encoded.push(item);//due to shift+1
            }
        }

        // we're done decoding the counts;
        // lets make sure we're flush with the u64 array
        // ie whatever is left in the current u64 is just zero-bits
        assert!(dd.is_clean());

        // update the position in the bytes array
        // we've read 8bytes per u64!
        let bytes_consumed = dd.get_consumed_bytes();
        self.pos += bytes_consumed;

        self.state = BuszBlockState::Flag;

        counts_encoded
    }

    fn parse_flags(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Flag);

        let mut dd = Self::fibonacci_factory(&self.buffer[self.pos..]);
        // let mut dd = U64Decoder::new(&self.buffer[self.pos..]);

        const FLAG_RLE_VAL: u64 = 0+1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut flag_decoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        while flag_decoded.len() < self.n_elements {

            let item = dd.next().unwrap();
            if item == FLAG_RLE_VAL {
                let runlength = dd.next().unwrap() as usize;
                // repeatedly add RLE
                flag_decoded.resize(
                    flag_decoded.len() + runlength,
                    FLAG_RLE_VAL - 1 //due to shift+1
                );
            } else {
                //decode single element
                flag_decoded.push(item-1); //due to shift+1
            }
        }

        assert!(dd.is_clean());

          // update the position in the bytes array
        // we've read 8bytes per u64!
        let bytes_consumed = dd.get_consumed_bytes();        
        self.pos += bytes_consumed;

        self.state = BuszBlockState::Finished;
        
        assert_eq!(self.pos, self.buffer.len(), "still leftover bytes in the buffer!"); 

        flag_decoded
    }

    /// parses the raw bytes in this block
    /// into a list of busrecords
    fn parse_block(&mut self) -> Vec<BusRecord>{

        let cbs = self.parse_cb();
        let umis = self.parse_umi(&cbs);
        let ecs = self.parse_ec();
        let count = self.parse_counts();
        let flag = self.parse_flags();

        let mut decoded_records = Vec::with_capacity(self.n_elements);
        for (cb,umi,ec,count,flag) in izip!(cbs, umis, ecs, count, flag) {
            decoded_records.push(
                BusRecord {
                    CB: cb, 
                    UMI: umi, 
                    EC: ec as u32,  
                    COUNT:count as u32,
                    FLAG: flag as u32
                }
            );
        }
        decoded_records
    }
}


#[test]
fn testing(){
    let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0 };
    let r2 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
    let r3 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 };
    let r4 = BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
    let r5 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0 };
    let r6 = BusRecord { CB: 1, UMI: 3, EC: 1, COUNT: 2, FLAG: 1 };

    let records = vec![
        r1.clone(),
        r2.clone(),
        r3.clone(),
        r4.clone(),
        r5.clone(),
        r6.clone(),
    ];

    // let mut w = BusWriter::new("/tmp/test.bus", BusParams {cb_len: 16, umi_len: 12});
    // w.write_iterator(records.iter().cloned());

    // let mut r = BuszReader::new("/tmp/test.busz");
    // let x = r.load_block_header().unwrap();
    // println!("ACTUAL FILE {:?}", x.get_blocksize_and_nrecords());

    let mut block_bytes = compress_busrecords_into_block(&records);

    let mut header_bytes =[0_u8; 8];
    for i in 0..8 {
        header_bytes[i] = block_bytes[i];
    }
    let h = CompressedBlockHeader { header_bytes: u64::from_le_bytes(header_bytes)};
    // println!("{:?}",h.get_blocksize_and_nrecords());

    // get rid of the header
    block_bytes = block_bytes[8..].to_vec();
    // println!("Bytes {:?}",block_bytes);
    
    // println!("Before Swap {:?}",block_bytes);
    // block_bytes = swap_endian(&block_bytes, 8);

    let (_blksize, nrecords) = h.get_blocksize_and_nrecords();

    println!("{}", block_bytes.len());
    // aggr to u64s
    assert_eq!(block_bytes.len() % 4, 0);


    let mut block = BuszBlock::new(&block_bytes, nrecords as usize);

    let records_dec = block.parse_block();

    assert_eq!(records, records_dec)
}


#[cfg(test)]
mod testing{
    use tempfile::tempdir;
    use crate::{busz::BuszWriter, iterators::CellGroupIterator};
    use super::*;

    #[test]
    /// this tests some weird assumption of the grouping iterators (.groupby_cb(), .groupby_cbumi()):
    /// They call .next() once more even thoguh they already recieved a `None`
    /// 
    fn test_groupby() {
        // odd! groupby_cbumi() calls next() once more after we're already done with the file
        // which leads to an exception
        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0 };
        let r2 = BusRecord { CB: 0, UMI: 2, EC: 1, COUNT: 2, FLAG: 0 };
        let r3 = BusRecord { CB: 0, UMI: 3, EC: 0, COUNT: 12, FLAG: 0 };
        let r4 = BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };
        let r5 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0 };
        let r6 = BusRecord { CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 };

        let records = vec![
            r1.clone(),
            r2.clone(),
            r3.clone(),
            r4.clone(),
            r5.clone(),
            r6.clone(),
        ];

        let dir = tempdir().unwrap();
        let file_path = dir.path().join("test.bus");
        let fname = file_path.to_str().unwrap();
        let mut w = BuszWriter::new(
            fname, 
            BusParams {cb_len:16, umi_len:12}, 
            4
        );
        w.write_iterator(records.into_iter());


        let mut reader = BuszReader::new(fname).groupby_cb();


        // lets step through it
        assert_eq!(reader.next(), Some((0, vec![r1,r2,r3])));
        assert_eq!(reader.next(), Some((1, vec![r4,r5])));
        assert_eq!(reader.next(), Some((2, vec![r6])));
        assert_eq!(reader.next(), None);

        // just to make sure it keeps yielding None
        assert_eq!(reader.next(), None);


    }

    #[test]
    fn test_inmem() {
        let fname = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.busz";
        let mut fh = BufReader::new(
            File::open(fname).unwrap()
        );
    
        let mut bytes: Vec<u8> = Vec::new();
        let bytes_read = fh.read_to_end(&mut bytes).unwrap();
        println!("bytes_read : {bytes_read}");

        let reader_inmem = BuszReader::from_read(bytes.as_slice());
        let reader_ondisk = BuszReader::new(fname);

        for (r1, r2) in izip!(reader_inmem, reader_ondisk) {
            assert_eq!(r1, r2)
        }
    }

    /// Make sure all BusRecord fields are decoded properly (we dont usualy spend time one FLAG)
    #[test]
    fn test_all_busrecord_fields() {

        let r1 = BusRecord { CB: 0, UMI: 0, EC: 1, COUNT: 2, FLAG: 5 };
        let r2 = BusRecord { CB: 100, UMI: 10, EC: 2, COUNT: 4, FLAG: 10 };
        let r3 = BusRecord { CB: 200, UMI: 20, EC: 3, COUNT: 8, FLAG: 15 };
        let r4 = BusRecord { CB: 300, UMI: 30, EC: 4, COUNT: 16, FLAG: 20 };
        let records = vec![r1, r2, r3, r4];


        let mut bytes: Vec<u8> = Vec::new();
        let mut w = BuszWriter::new(
            "/tmp/allfields.busz", 
            BusParams{ cb_len:16, umi_len: 12}, 
            3
        );
        w.write_iterator(records.clone().into_iter());

        let mut fh = File::open("/tmp/allfields.busz").unwrap();
        let bytes_read = fh.read_to_end(&mut bytes).unwrap();
        println!("bytes_read : {bytes_read}");


        let reader_inmem = BuszReader::from_read(bytes.as_slice());

        for (r1, r2) in izip!(reader_inmem, records) {
            assert_eq!(r1, r2)
        }
    }

    #[test]
    fn confirm_bus_oder() {
        let fname = "/tmp/tmp.bus";

        let r1 = BusRecord { CB: 6, UMI: 4, EC: 1, COUNT: 2, FLAG: 5 };
        let r2 = BusRecord { CB: 6+6, UMI: 4, EC: 1, COUNT: 2, FLAG: 5 };
        let mut w = BuszWriter::new(fname, BusParams { cb_len: 16, umi_len: 12 }, 100);
        w.write_records(vec![r1, r2]);
        w.terminal_flush();


        let mut r = BuszReader::new(fname);
        
        // the reader is now after the busheader, at the first blockheader

        let mut buf = [0_u8;8];
        r.reader.read_exact(&mut buf).unwrap();

        let h = CompressedBlockHeader { header_bytes: u64::from_le_bytes(buf)};
        println!("{:?}", h.get_blocksize_and_nrecords());

        let mut buf = [0_u8;8];
        r.reader.read_exact(&mut buf).unwrap();
        println!("CB {buf:?}, {:064b} -> u64 {}", u64::from_le_bytes(buf), u64::from_le_bytes(buf));


        let mut buf = [0_u8;8];
        r.reader.read_exact(&mut buf).unwrap();
        println!("UMI {buf:?}, {:064b} -> u64 {}", u64::from_le_bytes(buf), u64::from_le_bytes(buf));


        let mut buf = [0_u8;4];
        r.reader.read_exact(&mut buf).unwrap();
        println!("EC {buf:?} {:032b}", u32::from_le_bytes(buf));

        let mut buf = [0_u8;8];
        r.reader.read_exact(&mut buf).unwrap();
        println!("count {buf:?}");

    }
}