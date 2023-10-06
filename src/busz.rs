//! Dealing with the [busz compression format](https://github.com/BUStools/BUSZ-format)
//! 
//! # Examples
//! ## Reading a compressed bus file
//! ```rust, no_run
//! use bustools::busz::BuszReader;
//! let reader = BuszReader::new("/some/file.busz");
//! for record in reader {
//!     // ...
//! }
//! ```
//! ## Writing to a compressed bus file
//! ```rust, no_run
//! use bustools::busz::BuszWriter;
//! use bustools::io::{BusRecord, BusHeader};
//! let blocksize = 10000;
//! let mut writer = BuszWriter::new("/some/file.busz", BusHeader::new(16,12,0), blocksize);
//! let records = vec![
//!     BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0 },
//!     BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 },
//!     BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0 },
//!     BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 },
//!     BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0 },
//!     BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0 },
//! ];
//! writer.write_iterator(records.into_iter());
//! ```
//! 
//! # About Bitvec and Memory layout
//! This code relies heavily on BitVec. It uses [`bitvec`] to encode/decode
//! the bits of the busz records, in particular Fibbonnaci encoding and NewPFD encoding.
//! 
//! **A certain peculiarity though**:
//! To turn bytes (e.g from a u64 or read from the file) into [`bitvec::vec::BitVec`] we use `BitVec::from_bytes(byte_Array)`
//! This takes the bytes literally in the order of the array.
//! Yet `bustools` writes busz in little endian format, i.e. the byte order is reversed.
//! In particular, each busz block contains entries for CB,UMI... each PADDED with zeros afterwards(to a multiple of 64)
//! On disk this is how it looks like:
//! ```bash, no_run
//! 0000000...00000000[CBs in Fibbonnaci]
//! 0000000...00000000[UMIs in Fibbonnaci]
//! ```
//! 
//! Even more, the fibbonacci encoding must be done with little endian byte order, if on disk it looks like
//! ```bash,no_run
//! aaaaaaaabbbbbbbbccccccccddddddddeeeeeeeeffffffffgggggggghhhhhhhh  //bits
//! ```
//! the correct fibonacci stream to decode is
//! ```bash, no_run
//! ddddddddccccccccbbbbbbbbaaaaaaaahhhhhhhhgggggggg....
//! ``` 

use serde::{Serialize, Deserialize};
use self::utils::{setbits_u32, setbits_u64};

mod decode;
mod encode;
mod utils;
mod runlength_codec;


// exposing some core classes/functions to the public API
pub use encode::BuszWriter;
pub use decode::BuszReader;
pub use encode::compress_busfile;
pub use decode::decompress_busfile;

const PFD_BLOCKSIZE: usize = 512; // size of a PFD block within busz (this many ECs get encoded together)


#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
struct BuszSpecificHeader {
    block_size: u32,
    pfd_block_size: u32,
    lossy_umi: u32
}
impl BuszSpecificHeader {
    fn to_bytes(&self) -> Vec<u8>{
        bincode::serialize(self).expect("FAILED to serialze header")
    }
}

struct CompressedBlockHeader {
    // the 34 most significant bits denote the size of the compressed block in bytes. 
    // The 30 least significant bits denote the number of BUS records in the block.
    header_bytes: u64  
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
        let block_size_bytes = self.header_bytes >> bit_length;
        let bitmask_64 = setbits_u32(bit_length) as u64;
        let block_size_records = self.header_bytes & bitmask_64;
        (block_size_bytes, block_size_records)
    }
}


const BUSZ_HEADER_SIZE: usize = 4+4+4;
#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
struct BuszHeader {
    block_size: u32,
    pfd_block_size: u32,
    lossy_umi: u32,
}
impl BuszHeader {
    /// desearializes a `BusHeader` from Bytes; when reading busfiles
    /// assumes Little-Endian! [see here](https://docs.rs/bincode/latest/bincode/config/index.html#options-struct-vs-bincode-functions)
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
    /// assumes Little-Endian! [see here](https://docs.rs/bincode/latest/bincode/config/index.html#options-struct-vs-bincode-functions)
    pub fn to_bytes(&self) -> Vec<u8> {
        bincode::serialize(self).expect("FAILED to serialze header")
    }
}


#[cfg(test)]
mod test {
    use crate::busz::CompressedBlockHeader;
    
    #[test]
    fn test_header_encode_decode() {
        let nbytes = 20;
        let nrecords = 10;
        let h = CompressedBlockHeader::new(nbytes, nrecords);

        assert_eq!(h.get_blocksize_and_nrecords().0, nbytes);
        assert_eq!(h.get_blocksize_and_nrecords().1, nrecords);
    }

    mod external {
        use tempfile::tempdir;
        use pretty_assertions::assert_eq;
        use crate::io::{BusRecord, BusWriter, BusHeader, BusReader};
        use crate::busz::decode::{BuszReader, decompress_busfile};
        use crate::busz::encode::compress_busfile;

        #[test]
        fn test_external(){
            let v = vec![ 
                BusRecord {CB:10,UMI:11,EC:10,COUNT:13, FLAG: 14 },   // 10
                BusRecord {CB:11,UMI:11,EC:10,COUNT:13, FLAG: 14 },   // 0
                BusRecord {CB:22,UMI:10,EC:10,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:22,UMI:11,EC:10,COUNT:1, FLAG: 0 },    // 1
            ];

            let dir = tempdir().unwrap();
            let file_path= dir.path().join("buscompress.bus");
            let filename = file_path.to_str().unwrap();
            
            let mut  writer = BusWriter::new(
                filename, 
                BusHeader::new(16, 12, 1)
            );
            writer.write_records(&v);
        }

        #[test]
        fn test_encode_decode_busz(){
            let v = vec![ 
                BusRecord {CB:10,UMI:11,EC:10,COUNT:13, FLAG: 20 },   // 10
                BusRecord {CB:11,UMI:11,EC:10,COUNT:13, FLAG: 20 },   // 0
                BusRecord {CB:22,UMI:10,EC:10,COUNT:1, FLAG: 0 },   // 0
                BusRecord {CB:22,UMI:11,EC:10,COUNT:1, FLAG: 0 },    // 1
            ];

            // write plain bus
            // use tempfile::tempfile;
            let dir = tempdir().unwrap();
            let file_path= dir.path().join("buscompress.bus");
            let input_plain = file_path.to_str().unwrap();

            // let x = tempfile::tempfile().unwrap();
            // x.
            let mut  writer = BusWriter::new(
                input_plain, 
                BusHeader::new(16, 12, 1)
            );
            writer.write_records(&v);
            drop(writer);

            // copmress it
            let file_path= dir.path().join("lalalala.busz");
            let copmressed_output = file_path.to_str().unwrap();
            compress_busfile(
                input_plain,
                copmressed_output,
                100
            );

            // // decode it
            let reader = BuszReader::new(copmressed_output);
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
            let recs: Vec<_> = reader.collect();

            let x = BusReader::new(input_compressed_true);
            assert_eq!(x.collect::<Vec<_>>(), recs);

        }

        #[test]
        fn test_compress1() {
            // let input_compressed = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.busz"; 
            let input_plain = "/home/michi/bus_testing/bus_output_shorter/output.corrected.sort.bus";
            let copmressed_output = "/tmp/buscompress_testing.busz";
            compress_busfile(
                input_plain,
                copmressed_output,
                10000
            );
        }

        // #[test]
        fn test_compress_full() {
            // let input_compressed = "/home/michi/bus_testing/bus_output/output.corrected.sort.busz"; 
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
            reader.next();
            println!("==========================================");
            println!("==========================================");
            let mut reader = BuszReader::new("/tmp/lalalala.busz");

            reader.next();
            // let records:Vec<_> = reader.collect();
            // println!("{}", records.len())
        }
        #[test]
        fn test_decompress(){
            // decompress a busfile, check that the contents match the true (uncompressed version)
            let input_compressed = "/home/michi/bus_testing/bus_output/output.corrected.sort.busz"; 
            let input_plain = "/home/michi/bus_testing/bus_output/output.corrected.sort.bus";

            let dir = tempdir().unwrap();
            let file_path= dir.path().join("buscompress_lala.bus");
            let output = file_path.to_str().unwrap();

            let start = std::time::Instant::now();
            decompress_busfile(
                input_compressed,
                output);

            let elapsed = start.elapsed().as_millis();
            println!("decoding: {elapsed} ms");


            let r = BusReader::new(output);
            let records:Vec<_> = r.collect();
            let r_original = BusReader::new(input_plain);
            let records_original:Vec<_> = r_original.collect();

            assert_eq!(records.len(), records_original.len());
            assert_eq!(records, records_original);
        }

        #[test]
        fn test_iterator(){
            let reader = BuszReader::new("/home/michi/bus_testing/bus_output_shortest/output.corrected.sort.busz");
            let records:Vec<BusRecord> = reader.collect();

            let r_original = BusReader::new("/home/michi/bus_testing/bus_output_shortest/output.corrected.sort.bus");
            let records_original:Vec<_> = r_original.collect();

            assert_eq!(records.len(), records_original.len());
            assert_eq!(records, records_original);
        }   
    }

}