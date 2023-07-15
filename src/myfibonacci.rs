use bitvec::prelude::*;
use itertools::izip;

/// Iterative fibonacci.
///
/// <https://github.com/rust-lang/rust-by-example>
struct Fibonacci {
    curr: u64,
    next: u64,
}

impl Iterator for Fibonacci {
    type Item = u64;
    fn next(&mut self) -> Option<u64> {
        let new_next = self.curr + self.next;

        self.curr = self.next;
        self.next = new_next;

        Some(self.curr)
    }
}
/// A "constructor" for Iterative fibonacci.
fn iterative_fibonacci() -> Fibonacci {
    Fibonacci { curr: 1, next: 1 }
}


// not sure what the significance of those settings is
// in busz, converting byte buffers to BitSlices seems to require u8;Msb01
type MyBitSlice = BitSlice<u8, Msb0>;
type MyBitVector = BitVec<u8, Msb0>;

// let v: Vec<_> = iterative_fibonacci().take(65 - 1).collect();
// println!("{:?}", v);
const FIB64: &[u64]= &[1, 2, 3, 5, 8, 13, 21, 34, 55, 89, 144, 233, 377, 610, 987, 1597, 2584, 4181, 6765, 10946, 17711, 28657, 46368, 75025, 121393, 196418, 317811, 514229, 832040, 1346269, 2178309, 3524578, 5702887, 9227465, 14930352, 24157817, 39088169, 63245986, 102334155, 165580141, 267914296, 433494437, 701408733, 1134903170, 1836311903, 2971215073, 4807526976, 7778742049, 12586269025, 20365011074, 32951280099, 53316291173, 86267571272, 139583862445, 225851433717, 365435296162, 591286729879, 956722026041, 1548008755920, 2504730781961, 4052739537881, 6557470319842, 10610209857723, 17167680177565];

/// convert a bitslice holding a fibbonacci encoding into the numerical representation
pub fn bitslice_to_fibonacci(b: &MyBitSlice) -> u64{
    // omits the initial 1, i.e.
    // fib = [1,2,3,5,...]
    // let fib: Vec<_> = iterative_fibonacci().take(b.len() - 1).collect(); // shorten by one as we omit the final bit
    // println!("{:?}", fib);
    // b.ends_with(&[true, true].into());
    
    // TODO make sure its a proper fib-encoding (no 11 except the end)
    let mut sum = 0;
    for (bit, f) in izip!(&b[..b.len()-1], FIB64) {
        if *bit == true {
            sum+=f;
        }
    }
    sum
}

#[derive(Debug)]
pub struct MyFibDecoder <'a> {
    buffer: &'a MyBitSlice,
    current_pos: usize, // where we are at in the buffer (the last split), i.e the unprocessed part is buffer[current_pos..]
}

impl <'a> MyFibDecoder<'a> {
    pub fn new(buffer: &'a MyBitSlice) -> Self {
        MyFibDecoder { buffer, current_pos:0}
    }

    pub fn get_remaining_buffer(&self) -> &'a MyBitSlice{
        &self.buffer[self.current_pos..]
    }

    /// how far did we process into the buffer (pretty much the first bit after a 11)
    pub fn get_bits_processed(&self) -> usize{
        self.current_pos
    }
}

impl <'a> Iterator for MyFibDecoder<'a> {
    type Item=&'a MyBitSlice;

    fn next(&mut self) -> Option<Self::Item> {
        // let pos = self.current_pos;
        let mut lastbit = false;

        let current_slice = &self.buffer[self.current_pos..];
        // println!("currentslice {:?}", current_slice);
        for (idx, b) in current_slice.iter().enumerate() {
            if *b & lastbit {
                // found 11
                // let the_hit = Some(&self.buffer[self.current_pos..self.current_pos+idx+1]);
                let the_hit = &current_slice[..idx+1];
                self.current_pos += idx; 
                self.current_pos += 1;

                // let decoded = bitslice_to_fibonacci(the_hit);
                return Some(the_hit);
            }
            lastbit = *b
        }
        None
    }
}

#[cfg(test)]
mod test {
    use crate::myfibonacci::{bitslice_to_fibonacci, MyFibDecoder, MyBitVector};
    use bitvec::prelude::*;


    #[test]
    fn test_fib_(){
        let v: Vec<bool> = vec![1,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());
        assert_eq!(
            bitslice_to_fibonacci(&b),
            1
        );

        let v: Vec<bool> = vec![0,1,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());
        assert_eq!(
            bitslice_to_fibonacci(&b),
            2
        );
        let v: Vec<bool> = vec![0,0,1,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());
        assert_eq!(
            bitslice_to_fibonacci(&b),
            3
        );

        let v: Vec<bool> = vec![1,0,1,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());
        assert_eq!(
            bitslice_to_fibonacci(&b),
            4
        );

        let v: Vec<bool> = vec![1,0,0,0,0,1,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());
        assert_eq!(
            bitslice_to_fibonacci(&b),
            14
        );
        let v: Vec<bool> = vec![1,0,1,0,0,1,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());
        assert_eq!(
            bitslice_to_fibonacci(&b),
            17
        );
    }
    #[test]
    fn test_myfib_decoder() {
        let v: Vec<bool> = vec![0,0,1,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());

        // println!("full : {:?}", b);
        let mut my = MyFibDecoder {buffer: b.as_bitslice(), current_pos:0};

        assert_eq!(
            my.next(), 
            Some(BitVec::from_iter(vec![false, false, true, true].into_iter()).as_bitslice())
        );
        assert_eq!(
            my.next(), 
            None
        );
    }

    #[test]
    fn test_myfib_decoder_consecutive_ones() {
        let v: Vec<bool> = vec![0,0,1,1,1,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());

        println!("full : {:?}", b);
        let mut my = MyFibDecoder {buffer: b.as_bitslice(), current_pos:0};

        assert_eq!(
            my.next(), 
            Some(BitVec::from_iter(vec![false, false, true, true].into_iter()).as_bitslice())
        );
        assert_eq!(
            my.next(), 
            Some(BitVec::from_iter(vec![true, true].into_iter()).as_bitslice())

        );
    }

    #[test]
    fn test_myfib_decoder_nothing() {
        let v: Vec<bool> = vec![0,0,1,0,1,0,1].iter().map(|x|*x==1).collect();
        let b: MyBitVector  = BitVec::from_iter(v.into_iter());

        println!("full : {:?}", b);
        let mut my = MyFibDecoder {buffer: b.as_bitslice(), current_pos:0};

        assert_eq!(
            my.next(), 
            None
        );
    }
}