//! NewPFD with bitvec backend (instead of the old, bit_vec crate)
//! 
//! 
//! 

use bitvec::{prelude as bv, field::BitField};
use itertools::{izip, Itertools};

use crate::{myfibonacci::{self, bitslice_to_fibonacci}, busz::round_to_multiple, newpfd::NewPFDCodec};

/// elements in the primary buffer are stored using b_bits
/// decode these as u64s
fn decode_primary_buf_element(mut x: &bv::BitSlice<u8, bv::Msb0>) -> u64 {
    let a:u64 = x.load_be();
    a
}

 
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
fn decode_newpfdblock(buf: &bv::BitSlice<u8, bv::Msb0>, blocksize: usize) -> (Vec<u64>, usize) {

    // println!("********Start of NewPFD Block************");
    // println!("decode_newpfdblock \n{}", bitstream_to_string(buf));
    // The header, piece by piece

    let mut buf_position = 0;

    let mut fibdec = myfibonacci::MyFibDecoder::new(buf);
    let _b_bits = bitslice_to_fibonacci( fibdec.next().unwrap());
    let mut b_bits = _b_bits as usize;

    let mut min_el = bitslice_to_fibonacci( fibdec.next().unwrap());
    let mut n_exceptions = bitslice_to_fibonacci( fibdec.next().unwrap());

    b_bits -= 1;
    min_el -= 1;
    n_exceptions -= 1;
    
    // println!("Decoded Header b_bits: {b_bits} min_el: {min_el} n_exceptions: {n_exceptions}");

    // println!("Decoding gaps");
    let mut index_gaps = Vec::with_capacity(n_exceptions as usize);
    for _ in 0..n_exceptions { 
        let ix = bitslice_to_fibonacci( fibdec.next().unwrap()) - 1; // shift in encode
        index_gaps.push(ix);
    }

    // println!("Decoding exceptions");
    let mut exceptions = Vec::with_capacity(n_exceptions as usize);
    for _ in 0..n_exceptions { 
        let ex = bitslice_to_fibonacci( fibdec.next().unwrap());
        exceptions.push(ex);
    }

    let index: Vec<u64> = index_gaps
        .into_iter()
        .scan(0, |acc, i| {
            *acc += i;
            Some(*acc)
        })
        .collect();

    // println!("remain {:?}, len {}", x, x.len(), );
    // need to remove trailing 0s which where used to pad to a muitple of 32
    // let bits_after = x.len();
    // let delta_bits  = bits_before-bits_after;
    let delta_bits = fibdec.get_bits_processed();
    let padded_bits =  round_to_multiple(delta_bits, 32) - delta_bits;

    assert!(!buf[buf_position+delta_bits..buf_position + delta_bits + padded_bits].any());
    buf_position = buf_position + delta_bits + padded_bits;


    // the body of the block
    let buf_body = &buf[buf_position..];
    let mut body_pos = 0;
    // println!("buf_body \n{}", bitstream_to_string(buf_body));
    
    // note that a block can be shorter than sepcifiec if we ran out of elements
    // however, as currently implements (by bustools)
    // the primary block always has blocksize*b_bits bitsize!
    // i.e. for a partial block, we cant really know how many elements are in there
    // (they might all be 0)
    /*
    let n_elements = if (x.len() /b_bits) < blocksize {
        assert_eq!(x.len() % b_bits, 0); 
        let n_elem_truncate = x.len() / b_bits;
        n_elem_truncate
    } else {
        blocksize
    };
    */
    let n_elements = blocksize;
    let mut decoded_primary: Vec<u64> = Vec::with_capacity(n_elements);

    for i in 0..n_elements {
        // println!("++++++++++++++Element {}/{}++++++++++++*", i+1, n_elements);
        // println!("{:?}", x);
        // split off an element into x
        let bits = &buf_body[body_pos..body_pos+b_bits];
        body_pos+=b_bits;
        // println!("buf_pos {}",         body_pos);
        // println!("n_El {}",         i);
        // println!("prim bits {}", bitstream_to_string(bits));

        decoded_primary.push(decode_primary_buf_element(bits));
        // println!("remaining size {}", x.len())
    }
    // println!("********End of NEWPFD  Block************\n");
    // println!("decode Index: {:?}", index);
    // println!("decode Excp: {:?}", exceptions);
    // println!("decode prim: {:?}", decoded_primary);
    // println!("decode min: {}", min_el);

    // puzzle it together
    for (i, highest_bits) in izip!(index, exceptions) {
        let lowest_bits = decoded_primary[i as usize];
        // println!("high: {}", &display_u64_in_bits(highest_bits));
        // println!("low:  {}", &display_u64_in_bits(lowest_bits));

        let el = (highest_bits << b_bits) | lowest_bits;
        // println!("high_s{}", &display_u64_in_bits(highest_bits << b_bits));
        // println!("merge {}", &display_u64_in_bits(el));

        // println!("i:{i}, over: {highest_bits} -> {el}");
        let pos = dbg!(i as usize);
        decoded_primary[pos] = el;
    }

    //shift up againsty min_element
    let decoded_final: Vec<u64> = decoded_primary.iter().map(|x|x+min_el).collect();

    // println!("!!!remaining size {}, {:?}", x.len(), x);
    // println!("Final decode {:?}", decoded_final);
    (decoded_final, buf_position+body_pos)
}


pub fn decode(newpfd_buf: &bv::BitSlice<u8, bv::Msb0>, n_elements: usize, blocksize: usize) -> (Vec<u64>, usize){

    let mut pos = 0;
    let mut elements: Vec<u64> = Vec::with_capacity(n_elements);
    while elements.len() < n_elements {
        // println!("Decoding newPFD Block {}: {:?}, len={}", i, encoded, encoded.len());
        // each call shortens wth encoded BitVec

        let current_block = &newpfd_buf[pos..];
        let (els, bits_consumed) = decode_newpfdblock(current_block, blocksize);

        pos+= bits_consumed;
        // println!("----remaining size {}, {:?}", encoded.len(), encoded);
        // println!("Decoded {} elements", els.len());

        for el in els {
            elements.push(el);
        }
    }

    // trucate, as we retrieved a bunch of zeros from the last block
    elements.truncate(n_elements);

    (elements, pos)
}

fn bitstream_to_string(buffer: &bv::BitSlice<u8, bv::Msb0>) -> String{
    let s = buffer.iter().map(|x| if *x==true{"1"} else {"0"}).join("");
    s
}

#[test]
fn test_larger_ecs_22() {
    let input = vec![264597, 760881, 216982, 203942, 218976];
    // let input = vec![10, 22, 12, 13, 14];
    // let input = vec![1,2,3,4,5];
    let n = NewPFDCodec::new(32);

    // println!("{:?}", NewPFDBlock::get_encoding_params(&input, 0.9));
    let encoded = n.encode(input.iter().cloned());
    // println!("Encoded:\n{:?}", encoded);
    let encoded_bv: bv::BitVec<u8, bv::Msb0> = bv::BitVec::from_iter(encoded.iter());
    // println!("Encoded bv:\n{}", bitstream_to_string(&encoded_bv));


    let (decoded, _) = decode(&encoded_bv.as_bitslice(), 5, 32);
    
    assert_eq!(decoded, input);
}