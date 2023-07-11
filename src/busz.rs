//! Dealing with the [busz compression format](https://github.com/BUStools/BUSZ-format)

use itertools::Itertools;
use serde::{Serialize, Deserialize};

use crate::io::BusRecord;

#[derive(Serialize, Deserialize, Debug, PartialEq, Eq, Clone)]
pub struct BuszSpecificHeader {
    block_size: u32,
    pfd_block_size: u32,
    lossy_umi: u32
}

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

use fibonacci_codec::Encode;

fn compress_barcodes2(records: &[BusRecord]) -> bit_vec::BitVec {
    let runlength_codec = RunlengthCodec {RLE_VAL: 0};

    let mut cb_iter = records.iter().map(|r| r.CB);
    let mut last_el = cb_iter.next().unwrap();
    let mut delta_encoded = Vec::new();
    delta_encoded.push(last_el);
    for el in cb_iter{
        delta_encoded.push(el-last_el);
        last_el=el;
    }
    println!("Delta enc: {:?}",delta_encoded);
    let runlen_encoded = runlength_codec.encode(delta_encoded.into_iter());
    println!("run enc: {:?}",runlen_encoded);

    //fibbonaci encoding
    dbg!(&runlen_encoded);
    let enc = dbg!(runlen_encoded.fib_encode().unwrap());
    enc
}
// TODO: the first value encoded is a little funny, it gets incremented 2x!!
// periodic runlength encoding since the UMI resets to smaller values (when changing CB)
// and we cant handle negative differences!
fn compress_umis(records: &[BusRecord]) -> bit_vec::BitVec {
    let runlength_codec = RunlengthCodec {RLE_VAL: 0};
    let mut periodic_delta_encoded = Vec::new();
    let iii = records.iter();
    let mut last_record = &records[0];
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
    println!("Delta enc: {:?}",periodic_delta_encoded);
    let runlen_encoded = runlength_codec.encode(periodic_delta_encoded.into_iter());
    println!("run enc: {:?}",runlen_encoded);
    //fibbonaci encoding
    dbg!(&runlen_encoded);
    let enc = dbg!(runlen_encoded.fib_encode().unwrap());
    enc
}

// enum RLEItem{
//     RLE(u64),  // a run of the specified element, its length encoded
//     Other(u64)  // some other item (type u64) only a single occurance
// }

/// Compress barcodes of rows using delta-runlen(0)-fibonacci encoding
fn compress_barcodes(records: &[BusRecord]) -> bit_vec::BitVec {
    let RLE_VAL = 0;
    let mut last_bc = 0;
    let mut runlen = 0;

    let mut to_encode: Vec<u64> = Vec::new();
    for r in records {
        println!("----------------",);
        let barcode = dbg!(r.CB) - dbg!(last_bc); // delta encoding
        println!("{:?}, delta {}", r, barcode);


        if barcode == 0 {
            println!("Extending {} run of {}", runlen, r.CB);
            runlen+= 1;
        }
        else {
            println!("Finished {} run of {}",runlen, r.CB);
            // finish the previous run and encode it
            // we're encoding a runlength (char 0 for x times)
            if runlen > 0 {
                // encode the value
                to_encode.push(RLE_VAL+1); // since we cant encode zero 

                // encode the runlength
                to_encode.push(runlen);

                runlen = 0
            }
            // we're just encoding a single value
            println!("adding single item {}",barcode);

            to_encode.push(barcode+1); // since we cant encode zero 
        }
        last_bc = r.CB
    }
    println!("Aftermath");
    if runlen >0 {
        // encode the value
        to_encode.push(RLE_VAL+1); // since we cant encode zero 

        // encode the runlength
        to_encode.push(runlen);    
    }

    //actual encoding
    dbg!(&to_encode);
    let enc = dbg!(to_encode.fib_encode().unwrap());
    enc
}

/// Run length encoder/decoder
/// It only encodes the runlength of a special item (RLE_VAL), all other values are encoded as is
/// Encoding will yield a stream of u64s, encoding either (RLE_item, runlength) or (items)
struct RunlengthCodec {
    RLE_VAL : u64 // 
}
impl RunlengthCodec {
    // this currenly shifts all values by +1 (since the subsequent fibonacci encoding cant handle 0)
    pub fn encode(&self, input: impl Iterator<Item=u64>) -> Vec<u64>{
        let mut encoded: Vec<u64> = Vec::new();
        let mut runlen = 0;

        for x in input {
            if x == 0 {
                runlen+= 1;
            }
            else {
                // finish the previous run and encode it
                // we're encoding a runlength (char 0 for x times)
                if runlen > 0 {
                    // encode the value
                    encoded.push(self.RLE_VAL+1); // since we cant encode zero 
    
                    // encode the runlength
                    encoded.push(runlen);
    
                    runlen = 0
                }
                // we're just encoding a single value
                encoded.push(x+1); // since we cant encode zero 
            }
        }
        println!("Aftermath");
        if runlen >0 {
            // encode the value
            encoded.push(self.RLE_VAL+1); // since we cant encode zero 
    
            // encode the runlength
            encoded.push(runlen);    
        }
        encoded
    }

    // currently shifts all values by -1 
    pub fn decode(&self, input: Vec<u64>) -> Vec<u64>{
        let mut decoded = Vec::new();
        let mut iii = input.iter();
        loop {
            // if there's still some item in the stream
            if let Some(item) = iii.next() {
                println!("{}", item);
                if item == &(self.RLE_VAL+ 1) { // everything is shifted by 1 in the encoded stream
                    let runlen = *iii.next().unwrap(); // this shouldnt fail, each RLE is followed by runlen
                    for _ in 0..runlen {
                        decoded.push(self.RLE_VAL);  // shifting the value by -1, or equivalently use RLE
                    }
                } else {
                    decoded.push(*item-1)
                }
            }
            else{
                break
            }
        }
        decoded
    }
}




#[cfg(test)]
mod test {
    use crate::{busz::{CompressedBlockHeader, compress_barcodes2, compress_umis}, io::BusRecord};

    use super::compress_barcodes;

    mod RLE_Codec {
        use crate::busz::RunlengthCodec;

        #[test]
        fn test_encode_decode_single_run(){
            // a CODEC which compresses runs of thevalue 0
            let c = RunlengthCodec { RLE_VAL: 0 };

            let plain = vec![0,0,0];
            let enc = c.encode(plain.clone().into_iter());
            assert!(plain.len()> enc.len());

            let dec = c.decode(enc);
            assert_eq!(plain, dec)
        }
        #[test]
        fn test_encode_decode_no_run(){
            // a CODEC which compresses runs of thevalue 0
            let c = RunlengthCodec { RLE_VAL: 0 };
            let plain = vec![1,2,1];
            let enc = c.encode(plain.clone().into_iter());
            assert!(plain.len()== enc.len()); // cant be compressed

            let dec = c.decode(enc);
            assert_eq!(plain, dec)
        }
        #[test]
        fn test_encode_decode_minxed_run(){
            // a CODEC which compresses runs of thevalue 0
            let c = RunlengthCodec { RLE_VAL: 0 };
            let plain = vec![1,0,0, 2,0,1,0];
            let enc = c.encode(plain.clone().into_iter());
            assert!(plain.len()< enc.len()); // here the encoding is acually longer! dueto the freqeunt 1-runs of zeros

            let dec = c.decode(enc);
            assert_eq!(plain, dec)
        }
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
        let enc = compress_barcodes(&v);
        let decoded: Vec<_> = fibonacci_codec::fib_decode_u64(enc).map(|x|x.unwrap()).collect();
        
        assert_eq!(decoded, vec![
            1,2,  // two zers
            2,    // a single 1
            1,3   // three zero
            ]);

        let enc2 = compress_barcodes2(&v);
        assert_eq!(
            compress_barcodes(&v), 
            enc2
        );
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
}