use std::{fs::File, io::{BufWriter, Write}};

use crate::{io::{BusRecord, DEFAULT_BUF_SIZE, BusReader, BusHeader}, busz::{utils::{bitslice_to_bytes, swap_endian}, PFD_BLOCKSIZE, CompressedBlockHeader}};
use bitvec::prelude as bv;
use itertools::Itertools;
use newpfd::fibonacci::fib_enc_multiple_fast;
use super::{runlength_codec::RunlengthCodec, utils::round_to_multiple, BuszSpecificHeader};


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
    let mut enc = fib_enc_multiple_fast(&runlen_encoded);


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
    let mut enc= fib_enc_multiple_fast(&runlen_encoded);

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
    let mut enc = fib_enc_multiple_fast(&runlen_encoded);

    // pad to next multiple of 64
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
    let mut enc = fib_enc_multiple_fast(&runlen_encoded);


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

/// Compress `input` busfile into `output` busz-file using `blocksize`
/// 
/// # Parameters
/// * blocksize: How many elements are grouped together and compressed together
pub fn compress_busfile(input: &str, output: &str, blocksize: usize) {

    let reader = BusReader::new(input);
    let mut  writer = BuszWriter::new(output, reader.bus_header.clone(), blocksize);
    writer.write_iterator(reader.into_iter());
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

    /// Writes an iterator of Busrecords into a compressed busfile on disk
    /// This is the preferred way of using BusZWriter as it guarantees
    /// proper EOF and closing the compressed file
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

#[cfg(test)]
mod test {
    use crate::{io::BusRecord, busz::{encode::{compress_barcodes2, compress_umis, compress_ecs, BuszWriter, BuszWriterState}, decode::BuszReader}};

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
        // let decoded: Vec<_> = fibonacci_codec::fib_decode_u64(enc).map(|x|x.unwrap()).collect();
        let decoder = newpfd::fibonacci::FibonacciDecoder::new(&enc);
        let decoded: Vec<_> = decoder.collect();
        
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
        // let decoded: Vec<_> = fibonacci_codec::fib_decode_u64(enc).map(|x|x.unwrap()).collect();
        let decoder = newpfd::fibonacci::FibonacciDecoder::new(&enc);
        let decoded: Vec<_> = decoder.collect();


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
        // let decoded: Vec<_> = fibonacci_codec::fib_decode_u64(enc).map(|x|x.unwrap()).collect();
        let decoder = newpfd::fibonacci::FibonacciDecoder::new(&enc);
        let decoded: Vec<_> = decoder.collect();
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

    mod buszwriter {
        use crate::{io::{BusRecord, BusHeader}, busz::{encode::{BuszWriter, BuszWriterState}, decode::BuszReader}};

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
            assert_eq!(reader.collect::<Vec<BusRecord>>(), v);
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

            let buszfile = "/tmp/busz_writer_test2.buzs";
            let header = BusHeader::new(16,12,0);
            let mut writer = BuszWriter::new(buszfile, header , blocksize);
            writer.write_records(v.clone());
            writer.terminal_flush();

            assert_eq!(writer.state, BuszWriterState::FlushedAndClosed);
            assert_eq!(writer.internal_buffer.len(), 0);

            let reader = BuszReader::new(buszfile);
            assert_eq!(reader.collect::<Vec<BusRecord>>(), v);
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
            let buszfile = "/tmp/busz_writer_test3.buzs";
            let header = BusHeader::new(16,12,0);
            let mut writer = BuszWriter::new(buszfile, header , blocksize);
            writer.write_iterator(v.iter().cloned());


            assert_eq!(writer.state, BuszWriterState::FlushedAndClosed);
            assert_eq!(writer.internal_buffer.len(), 0);
            
            let reader = BuszReader::new(buszfile);

            assert_eq!(reader.collect::<Vec<BusRecord>>(), v);
        }
    }
}

