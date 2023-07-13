
// enum RLEItem{
//     RLE(u64),  // a run of the specified element, its length encoded
//     Other(u64)  // some other item (type u64) only a single occurance
// }

/// Run length encoder/decoder
/// It only encodes the runlength of a special item (RLE_VAL), all other values are encoded as is
/// Encoding will yield a stream of u64s, encoding either (RLE_item, runlength) or (items)
pub struct RunlengthCodec {
    #[allow(non_snake_case)]
    pub RLE_VAL : u64 // 
}
impl RunlengthCodec {
    // this currenly shifts all values by +1 (since the subsequent fibonacci encoding cant handle 0)
    pub fn encode(&self, input: impl Iterator<Item=u64>) -> Vec<u64>{
        let mut encoded: Vec<u64> = Vec::new();
        let mut runlen = 0;

        for x in input {
            if x == self.RLE_VAL {
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
        // println!("Aftermath");
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
                // println!("{}", item);
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
    use crate::runlength_codec::RunlengthCodec;

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

        let dec = c.decode(dbg!(enc));
        assert_eq!(plain, dec)
    }

    #[test]
    fn test_encode_decode_minxed_run_rle1(){
        // a CODEC which compresses runs of thevalue 0
        let c = RunlengthCodec { RLE_VAL: 1 };
        let plain = vec![0,1,1, 1,1, 0];
        let enc = c.encode(plain.clone().into_iter());
        let dec = c.decode(dbg!(enc));
        assert_eq!(plain, dec)
    }
}