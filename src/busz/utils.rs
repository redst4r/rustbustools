use bitvec::{prelude as bv, field::BitField};
use itertools::Itertools;

/// turn a bitslice into an array of bytes
/// the first 8 bits (bits[..8]) will become the first byte in the result
/// i.e. a sort of BigEndian encoding
pub (crate)  fn bitslice_to_bytes(bits: &bv::BitSlice<u8, bv::Msb0>) -> Vec<u8>{

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


pub (crate) fn bitstream_to_string(buffer: &bv::BitSlice<u8, bv::Msb0>) -> String{
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

// // turn a bitvector (64 elements) into a u64
pub fn bits64_to_u64(x: bit_vec::BitVec) -> u64{
    assert_eq!(x.len(), 64);
    u64::from_le_bytes(x.to_bytes()[..8].try_into().unwrap())
}

/// swaps endianness of the byte-vector
/// assuming 8byte (u64) words
/// simply grabs each 8byte section and reverses the order within
/// # params
/// * wordsize: 4 (bytes) for u32, 8(bytes) for u64
pub fn swap_endian(bytes: &[u8], wordsize: usize) -> Vec<u8>{
    let mut swapped_endian: Vec<u8> = Vec::with_capacity(bytes.len());
    for bytes in bytes.chunks(wordsize){
        swapped_endian.extend(bytes.iter().rev());
    }
    swapped_endian
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
/// ```bash, no_run  // TODO STUPID, wont run on private modules
///  use bustools::busz::utils::round_to_multiple;
///  assert_eq!(round_to_multiple(10,10), 10);
///  assert_eq!(round_to_multiple(11,10), 20);
///  assert_eq!(round_to_multiple(6,5), 10);
/// ```
pub fn round_to_multiple(i: usize, multiple: usize) -> usize {
    ((i+multiple-1)/multiple)*multiple
}

pub fn setbits_u32(x: u8) -> u32 {
    u32::MAX >> (32 - x)
}

pub fn setbits_u64(x: u8) -> u64 {
    u64::MAX >> (64 - x)
}


#[cfg(test)]
mod test {
    use crate::busz::utils::swap_endian;

    #[test]
    fn endian_swapping() {
        let v = vec![0_u8,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15];
        let a = swap_endian(&v, 8);
        let b = swap_endian(&a, 4);
        let c = swap_endian(&b, 4);
        let d = swap_endian(&c, 8);
        assert_eq!(v,d);
    }

}