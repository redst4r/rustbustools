use core::panic;

use bit_vec::BitVec;
use fibonacci_codec::{Encode};
use itertools::{Itertools, izip};

use crate::{busz::display_u64_in_bits};

/// # TODO
/// - the mess with the big-endian/little endian and MSB/LSB stuff: it works currently, but hard to read
/// -  encoding the primary buffer (where everything fits within b_bits) using `bitpacker` crate!

#[derive(Debug)]
struct NewPFDParams {
    b_bits: usize,
    min_element: u64
}

pub struct NewPFDCodec {
    blocksize: usize,
}

impl NewPFDCodec {

    pub fn new(blocksize: usize) -> Self {
        NewPFDCodec { blocksize }
    }

    pub fn encode(&self, input_stream: impl Iterator<Item=u64>) -> BitVec {

        let mut encoded_blocks : Vec<BitVec> = Vec::new();
        let mut block_counter = 0;
        // chunk the input into blocks
        for b in &input_stream.chunks(self.blocksize) {

            let block_elements: Vec<_> = b.collect();

            println!("Encoding block {block_counter} with {} elements: {:?}", block_elements.len(), block_elements);
            block_counter+=1;

            let params = NewPFDBlock::get_encoding_params(&block_elements, 0.9);
            let mut enc = NewPFDBlock::new(params.b_bits, self.blocksize);
            let enc_block = enc.encode(&block_elements, params.min_element as u64);
            encoded_blocks.push(enc_block);
        }
        println!("Concat blocks: {:?}", encoded_blocks);
        // concat
        let mut iii = encoded_blocks.into_iter();
        let mut merged_blocks = iii.next().unwrap();
        for mut b in iii {
            merged_blocks.append(&mut b);
        }
        merged_blocks
    }

    pub fn decode(&self, mut encoded: BitVec) -> Vec<u64> {
        let mut elements: Vec<u64> = Vec::new();
        let mut i = 0;
        println!("Decoding Blocks: {:?}, len={}", encoded, encoded.len());
        while encoded.len() > 0 {
            println!("Decoding Block {}: {:?}, len={}", i, encoded, encoded.len());
            i+=1;
            // each call shortens wth encoded BitVec
            let (els, remaining_buffer) = decode_newpfdblock(&mut encoded, self.blocksize);
            println!("----remaining size {}, {:?}", encoded.len(), encoded);

            for el in els {
                elements.push(el);
            }
            println!("Decoded {} elements", elements.len());
            encoded = remaining_buffer;

        }
        elements
    }
}

/// A Block of NewPFD compressed data
struct NewPFDBlock {
    // blocksize: usize,
    b_bits: usize,  // The number of bits each num in `pfd_block` is represented with.
    primary_buffer: Vec<bit_vec::BitVec>,
    exceptions: Vec<u64>,
    index_gaps : Vec<u64>
}

impl NewPFDBlock {

    pub fn new(b_bits: usize, blocksize: usize) -> Self {
        let pb = Vec::with_capacity(blocksize);
        let exceptions: Vec<u64> = Vec::new();
        let index_gaps : Vec<u64> = Vec::new();
        NewPFDBlock { 
            // blocksize: blocksize, 
            b_bits: b_bits, 
            primary_buffer: pb, 
            exceptions: exceptions, 
            index_gaps: index_gaps 
        }
    }

    /// turns input data into the NewPFD storage format:
    /// represent items by b_bits, store any expections separately 
    // pub fn from_data(input: &[u64], b_bits: usize) -> Self{

    // }

    /// determine the bitwidth used to store the elements of the block
    /// and the minimum element of the block 
    pub fn get_encoding_params(input: &[u64], percent_threshold: f64) -> NewPFDParams {

        let mut percent_ix = (input.len() as f64 ) * percent_threshold;
        percent_ix -= 1.0; // e..g if len = 4, percent=0.75  percent_ix should be pointing to the 3rd element, ix=2
        let ix = percent_ix.round() as usize;
        
        let (minimum, mut n_bits) = if ix == 0 {
            let the_element = input[0];
            (the_element, u64::BITS - the_element.leading_zeros())
        } else {
            let mut sorted = input.iter().sorted();
            let minimum = *sorted.next().unwrap();
            let the_element = sorted.nth(ix-1).unwrap();   //-1 since we took the first el out 
            let n_bits = u64::BITS - the_element.leading_zeros(); 
            (minimum, n_bits)
        };

        // weird exception: IF there's only a single bit to encode (blocksize==1)
        // and that bit happens to be zero -> n_bits==0 which will cause problems down the road
        // hence set n_bits>=1
        if n_bits == 0 {
            n_bits = 1;
        } 
        NewPFDParams { b_bits: n_bits as usize, min_element: minimum}
    }

    pub fn encode(&mut self, input: &[u64], min_element: u64) -> BitVec {

        // assert!(input.len() <= self.blocksize);

        // this is the maximum element we can store using b_bits
        let max_elem_bit_mask = (1_u64 << self.b_bits) - 1;

        // index of the last exception we came across
        let mut last_ex_idx = 0;

        for (i,x) in input.iter().enumerate() {
            let mut diff = x- min_element; // all elements are stored relative to the min

            // if its an exception, ie to big to be stored in b-bits
            // we store the overflowing bits seperately
            // and remember where this exception is
            if diff > max_elem_bit_mask {
                println!("Exception {}", diff);
                self.exceptions.push(diff >> self.b_bits);
                self.index_gaps.push((i as u64)- last_ex_idx);
                last_ex_idx = i as u64;

                // write the stuff that fits into the field
                diff &= max_elem_bit_mask;  // bitwise and with 111111...1111 to get the `b_bits` least significant bits
            }

            // put the rest into the primary buffer
            // turn the element to be stored into a bitvec (it'll be bigger than intended, but we'll chop it)
            // BigEndian: least significant BYTES are on the right!
            let mut bvec = bit_vec::BitVec::from_bytes(&diff.to_be_bytes());
            println!("Diff in bytes: {:?}", diff.to_be_bytes());

            let cvec = bvec.split_off(bvec.len()-self.b_bits);  // splits off the higest bits (which should be zero anyway!)
            assert_eq!(cvec.len(), self.b_bits);
            assert!(!bvec.any()); // the chopped off part should be 0
            println!("{:?}", cvec);
            self.primary_buffer.push(cvec);
        }
        println!("encode: Primary {:?}", self.primary_buffer);
        println!("encode: b_bits {}", self.b_bits);
        println!("encode: min_el {}", min_element);
        println!("encode: n_ex {}", self.exceptions.len());
        println!("encode: Exceptions {:?}", self.exceptions);
        println!("encode: Gaps: {:?}", self.index_gaps);

        // merge the primary buffer into a single bitvec, the body of the block
        let mut body = BitVec::with_capacity(self.b_bits * self.primary_buffer.len());
        for mut b in self.primary_buffer.iter_mut() { 
            body.append(&mut b);
        }
        // println!("{:?}", body);

        // now, as the BlockHeader in front of the merged_primary buffer we store via fibonacci encoding:
        // 1. b_bits
        // 2. min_element
        // 3.n_exceptions
        // 4. index_gaps
        // 5. exceptions
        let mut to_encode: Vec<u64> = vec![
            1 + self.b_bits as u64, 
            1+ min_element, 
            1+ self.exceptions.len() as u64
        ];
        // adding the exceptions
        to_encode.extend(self.index_gaps.iter().map(|x| x+1)); //shifting the gaps +1
        to_encode.extend(self.exceptions.iter());

        // merge header + body
        let mut header = to_encode.fib_encode().unwrap();
        header.extend(body);
        header
    }

 
}

// // turn a bitvector (64 elements) into a u64
// fn bits64_to_u64(x: BitVec) -> u64{
//     assert_eq!(x.len(), 64);
//     u64::from_le_bytes(x.to_bytes()[..8].try_into().unwrap())
// }

/// Decoding a block of NewPFD from a BitVec containing a series of blocks
/// 
/// This pops off the front of the BitVec, removing the block from the stream
/// 
/// 1. decode (Fibbonacci) metadata+expeptions+gaps
/// 2. Decode `blocksize` elements (each of sizeb_bits)
/// 3. aseemble the whole thing, filling in expecrionts etc
/// 
/// We need to know the blocksize, otherwise we'd start decoding the header of 
/// the next block
/// 
/// # NOte:
/// The bitvector x gets changed in this function. Oddly that has weird effects on this 
/// variable outside the function (before return the x.len()=9, outside the function it is suddenly 3)
/// Hence we return the remainder explicitly
pub fn decode_newpfdblock(x: &mut BitVec, blocksize: usize) -> (Vec<u64>, BitVec) {

    println!("********Start of Block************");
    println!("decode_newpfdblock {:?}", x);
    // The header, piece by piece
    let (mut _b_bits, mut x) = bitvec_single_fibbonacci_decode(x);
    let mut b_bits = _b_bits as usize;
    let (mut min_el, mut x) = bitvec_single_fibbonacci_decode(&mut x);
    let (mut n_exceptions, mut x) = bitvec_single_fibbonacci_decode(&mut x);
    b_bits -= 1;
    min_el -= 1;
    n_exceptions -= 1;
    
    println!("Decoded Header {b_bits} {min_el} {n_exceptions}");

    println!("Decoding gaps");
    let mut index_gaps = Vec::with_capacity(n_exceptions as usize);
    for _ in 0..n_exceptions { 
        let (mut ix, x2) = bitvec_single_fibbonacci_decode(&mut x);
        ix -= 1;  // shift in encode
        x = x2;
        index_gaps.push(ix);
    }

    println!("Decoding exceptions");
    let mut exceptions = Vec::with_capacity(n_exceptions as usize);
    for _ in 0..n_exceptions { 
        let (ex, x2) = bitvec_single_fibbonacci_decode(&mut x);
        exceptions.push(ex);
        x = x2;
    }
    let index: Vec<u64> = index_gaps
        .into_iter()
        .scan(0, |acc, x| {
            *acc += x;
            Some(*acc)
        })
        .collect();

    // the body of the block
    println!("decoding body {x:?}");

    // note that a block can be shorter than sepcifiec if we ran out of elements
    let n_elements = if (x.len() /b_bits) < blocksize {
        assert_eq!(x.len() % b_bits, 0); 
        let n_elem_truncate = x.len() / b_bits;
        n_elem_truncate
    } else {
        blocksize
    };

    let mut decoded_primary: Vec<u64> = Vec::with_capacity(n_elements);

    for i in 0..n_elements {
        println!("++++++++++++++Element {}/{}++++++++++++*", i+1, n_elements);
        println!("{:?}", x);
        // split off a block into x
        let second_half = x.split_off(b_bits);
        println!("x_split {:?}, len={}", x, x.len());
        // note that .to_Bytes will pad the bits to a multiple of 8! trailing bits will be added as 0
        // if x == 0001 = 1 
        //  pad to 00010000  =16
        // if x == 0010 = 2
        //  pad to 00100000  =32

        // hence lets manually pad the bitvec on the left with zeros
        let n_pad = if x.len() < 8 {
            8-x.len()
        } else
        {
            x.len() % 8
        };
        let mut padded_x = BitVec::from_elem(n_pad, false);
        padded_x.append(&mut x);
        println!("padded_x {:?},len={}", padded_x, padded_x.len());


        let mut bytes = padded_x.to_bytes();
        println!("bytes {:?}, len={}", bytes, bytes.len());

        // fill with zero bytes
        let npad = 8-bytes.len();
        let mut pad: Vec<u8> = Vec::with_capacity(npad);
        for _ in 0..npad {
            pad.push(0)
        }

        if true {
            // little endian
            bytes.extend(pad);

            println!("bytes  padded {:?}", bytes);
            
            assert_eq!(bytes.len(), 8);
            decoded_primary.push(u64::from_le_bytes(bytes[..8].try_into().unwrap()));
        }
        else {
            // // big endian:
            pad.extend(bytes);
            println!("bytes padded {:?}", pad);
            assert_eq!(pad.len(), 8);
            decoded_primary.push(u64::from_le_bytes(pad[..8].try_into().unwrap()));
        }

        x = second_half;
        println!("remaining size {}", x.len())
    }
    println!("********End of Block************\n");
    println!("decode Index: {:?}", index);
    println!("decode Excp: {:?}", exceptions);
    println!("decode prim: {:?}", decoded_primary);
    println!("decode min: {}", min_el);

    // puzzle it together
    for (i, highest_bits) in izip!(index, exceptions) {
        let lowest_bits = decoded_primary[i as usize];
        println!("high: {}", &display_u64_in_bits(highest_bits));
        println!("low:  {}", &display_u64_in_bits(lowest_bits));

        let el = (highest_bits << b_bits) | lowest_bits;
        println!("high_s{}", &display_u64_in_bits(highest_bits << b_bits));
        println!("merge {}", &display_u64_in_bits(el));

        // println!("i:{i}, over: {highest_bits} -> {el}");
        let pos = dbg!(i as usize);
        decoded_primary[pos] = el;
    }

    //shift up againsty min_element
    let decoded_final: Vec<u64> = decoded_primary.iter().map(|x|x+min_el).collect();

    println!("!!!remaining size {}, {:?}", x.len(), x);
    println!("Final decode {:?}", decoded_final);
    (decoded_final, x)
}



/// decode a single number from the Bitvector, truncates the original bitvecotr
fn bitvec_single_fibbonacci_decode(x: &mut BitVec) -> (u64, BitVec){

    // find the location of the first occurance of 11
    let mut lastbit = false;
    let mut ix: Option<usize> = None;
    for (i, b) in x.iter().enumerate() {
        if lastbit & b { // we ran into a 11
            ix = Some(i);
            break;
        }
        lastbit = b
    }
    if let Some(i) = ix{
        let remainder = x.split_off(i + 1); // i marks the last `1`: 00011, but the  

        // println!("{x:?}");
        // println!("{:?}", remainder);

        let dec: Vec<_> = fibonacci_codec::fib_decode_u64(x.iter())
            .map(|a|a.unwrap())
            .collect();
        assert_eq!(dec.len(), 1);

        // println!("single dec {x:?} -> {}", dec[0]);

        (dec[0], remainder)
        // a
    }
    else {
        let error = format!("error in decoding {} {:?}",x.len(), x);
        panic!("{}", &error);
    }

}

#[cfg(test)]
mod test {
    use crate::newpfd::{decode_newpfdblock, NewPFDCodec};

    use super::NewPFDBlock;

    #[test]
    fn test_encode() {
        let mut n = NewPFDBlock::new(4,  10);
        let input = vec![0,1,2,16, 1, 17];
        n.encode(&input, 0);
    }

    #[test]
    fn test_params() {
        let input = vec![0,10, 100, 1000];
        let params = NewPFDBlock::get_encoding_params(&input, 0.75);
        assert_eq!(params.min_element, 0);
        assert_eq!(params.b_bits, 7);
    }

    #[test]
    fn test_params_min_el() {
        let n = NewPFDBlock::new(4,  10);
        let input = vec![1,10, 100, 1000];
        let params = dbg!(NewPFDBlock::get_encoding_params(&input, 0.75));
        assert_eq!(params.min_element, 1);
        assert_eq!(params.b_bits, 7);
    }


    #[test]
    fn test_encode_decode_less_elements_than_blocksize() {
        let blocksize = 10; 
        let mut n = NewPFDBlock::new(4,  blocksize);
        let input = vec![0,1,2,16, 1, 17, 34, 1];
        let mut encoded = n.encode(&input, 0);
        println!("Enc length {}", encoded.len());
        println!("Plain length {}", input.len()*64);


        let (decoded,_) = decode_newpfdblock(&mut encoded, blocksize);
        assert_eq!(decoded, input);
    }
    #[test]
    fn test_encode_no_exceptions() {
        let mut n = NewPFDBlock::new(2,  10);
        let input = vec![0,1,0, 1];
        let _encoded = n.encode(&input, 0);

        assert_eq!(n.exceptions.len(), 0);
        assert_eq!(n.index_gaps.len(), 0);
    }

    #[test]
    fn test_NewPFD_codec_encode_decode() {
        let n = NewPFDCodec::new(2);
        let input = vec![0_u64,1,0, 1];
        let encoded = n.encode(input.into_iter());
        println!("=============================");
        println!("=============================");
        println!("=============================");
        let decoded = n.decode(encoded);
        assert_eq!(decoded, vec![0_u64,1,0, 1]);
    }


    #[test]
    fn test_NewPFD_codec_encode_decode_blocksize1() {
        // blocksize==1 exposes some edge cases, like a single 0bit in the block
        let n = NewPFDCodec::new(1);
        let input = vec![0_u64,1,0, 1];
        let encoded = n.encode(input.into_iter());

        println!("=============================");
        println!("=============================");
        println!("=============================");
        let decoded = n.decode(encoded);
        assert_eq!(decoded, vec![0_u64,1,0, 1]);
    }

    #[test]
    fn test_NewPFD_codec_encode_decode_nonzero_min_el() {
        let n = NewPFDCodec::new(4);
        let input = vec![1_u64,2,2, 1];
        let encoded = n.encode(input.into_iter());
        let decoded = n.decode(encoded);
        assert_eq!(decoded, vec![1_u64,2,2, 1]);
    }


}