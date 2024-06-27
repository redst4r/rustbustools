use bitvec::field::BitField;
use itertools::Itertools;

use super::BuszBitSlice;

/// turn a bitslice into an array of bytes
/// the first 8 bits (bits[..8]) will become the first byte in the result
/// i.e. a sort of BigEndian encoding
pub (crate) fn bitslice_to_bytes(bits: &BuszBitSlice) -> Vec<u8>{

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

/// for debugging: display a BitVec as a string of bits
pub (crate) fn bitstream_to_string(buffer: &BuszBitSlice) -> String{
    let mut s = String::new();
    let x = buffer.iter().map(|x| if *x{"1"} else {"0"});
    for bit64 in &x.into_iter().chunks(64){
        let concat = bit64.collect::<Vec<_>>().join("");
        s.push_str(&concat);
        s.push('\n');
    }
    s
}

// all the encoded parts (cb,umi,ec...) are padded with zeros until a mutiple of 64
// which we need to remove
pub(crate) fn calc_n_trailing_bits(bits_processed: usize) -> usize {
    let padded_size = round_to_multiple(bits_processed, 64);
    let zeros_toremoved = padded_size - bits_processed;
    zeros_toremoved
}

// // // turn a bitvector (64 elements) into a u64
// pub fn bits64_to_u64(x: bit_vec::BitVec) -> u64{
//     assert_eq!(x.len(), 64);
//     u64::from_le_bytes(x.to_bytes()[..8].try_into().unwrap())
// }

/// swaps endianness of the byte-vector
/// assuming 8byte (u64) words
/// simply grabs each 8byte section and reverses the order within
/// # params
/// * wordsize: 4 (bytes) for u32, 8(bytes) for u64
pub (crate) fn swap_endian(bytes: &[u8], wordsize: usize) -> Vec<u8>{
    let mut swapped_endian: Vec<u8> = Vec::with_capacity(bytes.len());
    for bytes in bytes.chunks(wordsize){
        swapped_endian.extend(bytes.iter().rev());
    }
    swapped_endian
}

/// shortcut for doing swap_endian(swap_endian(x, 8), 4)
/// essentially does tihs:
/// ...ABCD|EFGH...: ....EFGH|ABCD...
pub (crate) fn swap_endian8_swap_endian4(bytes: &[u8], ) -> Vec<u8>{
    let mut swapped_endian: Vec<u8> = Vec::with_capacity(bytes.len());
    for bytes in bytes.chunks(8){
        swapped_endian.extend(&bytes[4..]);
        swapped_endian.extend(&bytes[..4]);
    }
    swapped_endian
}

/// does the endian8-endian4 swap, but without allocating anything
/// just does it in the vec itself
pub (crate) fn swap_endian8_swap_endian4_inplace(x: &mut [u8]) {
    // let n_chunks: usize = (x.len() / 8).try_into().expect("must be multiple of 8 (u64=8 bytes!)");
    assert_eq!(x.len() % 8, 0, "must be multiple of 8 (u64=8 bytes!");
    let n_chunks: usize = x.len() / 8;
    for i in 0..n_chunks {
        let pos = i * 8;

        x.swap(pos, pos+4);
        x.swap(pos+1, pos+5);
        x.swap(pos+2, pos+6);
        x.swap(pos+3, pos+7);
    }
}




// pub fn display_u64_in_bits(x: u64) -> String{
//     let s: Vec<u64> = (0..64).rev().map (|n| (x >> n) & 1).collect();
//     s.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join("")
// }

// pub fn display_u32_in_bits(x: u32) -> String{
//     let s: Vec<u32> = (0..32).rev().map (|n| (x >> n) & 1).collect();
//     s.into_iter().map(|x| x.to_string()).collect::<Vec<_>>().join("")
// }

/// round an integer to the next bigger multiple
/// ```bash, no_run  // TODO STUPID, wont run on private modules
///  use bustools::busz::utils::round_to_multiple;
///  assert_eq!(round_to_multiple(10,10), 10);
///  assert_eq!(round_to_multiple(11,10), 20);
///  assert_eq!(round_to_multiple(6,5), 10);
/// ```
pub (crate) fn round_to_multiple(i: usize, multiple: usize) -> usize {
    // ((i+multiple-1)/multiple)*multiple
    i.next_multiple_of(multiple)  // rust 1.73
}

/// set the lowest x bits in a 32bit vecotr (represented as u32)
pub (crate) fn setbits_u32(x: u8) -> u32 {
    u32::MAX >> (32 - x)
}

/// set the lowest x bits in a 64bit vecotr (represented as u32)
pub (crate) fn setbits_u64(x: u8) -> u64 {
    u64::MAX >> (64 - x)
}


#[cfg(test)]
mod test {
    use std::io::Read;
    use bitvec::{bits, order::Msb0};
    use fastfibonacci::utils::create_bitvector;
    use super::*;
    use crate::busz::utils::{swap_endian, setbits_u32, setbits_u64};
    #[test]
    fn test_bitslice_to_bytes() {
        let b = create_bitvector(vec![
            0,0,0,0, 0, 0, 1, 1,
            1,1,1,1, 1, 1, 1, 1,
        ]);
        
        assert_eq!(
            bitslice_to_bytes(&b),
            vec![3, 255]
        );

        assert_eq!(
            b.bytes().map(|x| x.unwrap()).collect::<Vec<_>>(),
            vec![3, 255]
        )
    }

    #[test]
    fn test_swap_inmem() {
        let b = bits![u64, Msb0; 
            0,0,0,0, 0, 0, 0, 1,
            0,0,0,0, 0, 0, 1, 1,
            0,0,0,0, 0, 1, 1, 1,
            0,0,0,0, 1, 1, 1, 1,
            0,0,0,1, 1, 1, 1, 1,
            0,0,1,1, 1, 1, 1, 1,
            0,1,1,1, 1, 1, 1, 1,
            1,1,1,1, 1, 1, 1, 1,
        ];
        let before = b.bytes().map(|x| x.unwrap()).collect::<Vec<_>>();

        println!("{:?}", before);


        let bytes_it = b.bytes().map(|x| x.unwrap());

        let mut swapped_endian: Vec<u8> = Vec::new();
        for bytes in bytes_it.chunks(8).into_iter(){
            let b1: Vec<_> = bytes.collect();
            swapped_endian.extend(&b1[4..]);
            swapped_endian.extend(&b1[..4]);
        }

        println!("{:?}", swapped_endian)
    }
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
    fn test_setbits_u32() {
        assert_eq!(setbits_u32(3), 7);
        assert_eq!(setbits_u32(2), 3);
        assert_eq!(setbits_u32(1), 1);
    }
    #[test]
    fn test_setbits_u64() {
        assert_eq!(setbits_u64(3), 7);
        assert_eq!(setbits_u64(2), 3);
        assert_eq!(setbits_u64(1), 1);

    }
    #[test]
    fn test_swap_endian8_swap_endian4() {
        let x = &[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16, 17,18,19,20,21,22,23,24];
        let y1 = swap_endian8_swap_endian4(x);
    
        let y2 = swap_endian(
            &swap_endian(x, 8),
            4
        );
        assert_eq!(y1, y2);
    
        let y3 = swap_endian(
            &swap_endian(x, 4),
            8
        );
        assert_eq!(y1, y3);
    
        // inplace swap
        let mut a: Vec<u8> = x.to_vec();
        // println!("before {:?}", a);
        swap_endian8_swap_endian4_inplace(&mut a);
        // println!("after {:?}", a);
        assert_eq!(y1, a);
    }

}