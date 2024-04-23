use std::{collections::VecDeque, fs::File, io::{BufReader, Read}};
use crate::{io::{BusHeader, BusRecord, BusWriterPlain, DEFAULT_BUF_SIZE, BUS_HEADER_SIZE, CUGIterator, BusParams}, busz::{BUSZ_HEADER_SIZE, utils::{swap_endian, calc_n_trailing_bits, bitstream_to_string}}};
use bitvec::prelude as bv;
use itertools::izip;
use fastfibonacci::FbDec;
use super::{utils::{bitslice_to_bytes, swap_endian8_swap_endian4, swap_endian8_swap_endian4_inplace}, BuszHeader, CompressedBlockHeader, PFD_BLOCKSIZE};

/// Reading a compressed busfile
/// 
/// The analog of [`crate::io::BusReader`] for compressed bus files
/// 
/// # Example
/// ```rust, no_run
/// use bustools::busz::BuszReader;
/// let reader = BuszReader::new("/some/file.busz");
/// for record in reader {
///     // ...
/// }
/// ```
/// ## Iterating over cells
/// ```rust, no_run
/// use bustools::busz::BuszReader;
/// use bustools::iterators::CellGroupIterator; //need to bring that trait into scope
/// let reader = BuszReader::new("/some/file.busz");
/// for (cb,record) in reader.groupby_cb() {
///     // ...
/// }
/// ```
/// 
/// 
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

        
        // conversion to big-endian to make the reading work right
        let bigendian_buf = swap_endian(&block_buffer, 8);

        let theblock: &bv::BitSlice<u8, bv::Msb0> = bv::BitSlice::from_slice(&bigendian_buf);
        let mut block = BuszBlock::new(theblock,nrecords as usize);
        
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
    buffer: &'a bv::BitSlice<u8, bv::Msb0>,  // TODO chould probably own it, we're doing some swaping and shifting, and it probablyh wont be usable after
    pos: usize, // where we are currently in the buffer
    n_elements: usize ,// how many busrecords are stored in the block
    state: BuszBlockState,
    debug: bool
}

impl <'a> BuszBlock <'a> {
    fn new(buffer: &'a bv::BitSlice<u8, bv::Msb0>, n_elements: usize) -> Self {
        // TODO: warning, buffer must be conveted in a special way if the bytes come out of a file
        // see BuszReader::load_busz_block_faster
        // cant do it in herer due to lifetime issues
        BuszBlock { 
            buffer, 
            pos: 0, 
            n_elements, 
            state: BuszBlockState::Cb,
            debug: false
        }
    }

    // just a streamlined way to generate a fibonacci decoder
    // which is used in parse_xxx()
    // fn fibonacci_factory(stream: &bv::BitSlice<u8, bv::Msb0>) -> fastfibonacci::fibonacci::FibonacciDecoder { 
    //     fastfibonacci::fibonacci::FibonacciDecoder::new(stream, false)
    // }

    fn fibonacci_factory(stream: &bv::BitSlice<u8, bv::Msb0>) -> fastfibonacci::fast::FastFibonacciDecoder<u8>{
        fastfibonacci::fast::get_u8_decoder(stream, false)
    }

    /// The whole issue with decoding is the way bustools stored the different parts (CB/UMI...)
    /// heres the layout
    /// 
    /// | CB     |   UMI   |    EC    |    COUNT   |   FLAG    |
    /// |  a*64  |   b*64  |c*32|c*32 |    e*64    |   f *64   |  size in bits
    /// 
    /// In other words the CB/UMI/COUNT/FLAG as stored in multiples of 64 (adding padding if needed)
    /// 
    /// When bustools saves to disk it uses little endian for all u64,u32
    /// The individual bytes look like this 
    /// | CB         |   UMI       |    EC     |    COUNT   |   FLAG    |
    /// |HGFEDCBA... |HGFEDCBA...  |4321 |8765 | HGFEDCBA   | HGFEDCBA ...
    /// 
    /// i.e. the stream of u8 we get from reading the busfile is DCBA...HGFE...
    /// to convert to a bitstream we swap *the entire stream* to big endian (since we cant now how long de segments are)
    /// | CB         |   UMI       |    EC     |    COUNT   |   FLAG    |
    /// |ABCDEFGH... | ABCDEFGH... | 5678 1234 | ABCDEFGH...
    /// Note how this screws up the order of the two u32s in EC (should really be 1234|5678 for proper iterating)
    ///
    /// Once we get to the start of the EC segment we have to undo this mess
    /// Initially: 
    /// - convert the bitvector (at this stage) to byte-vector, undo the 8byte endian swap and do a 4bye endian swap,convert Vec<u8> to bitvec
    /// - after the EC block,undo 4byte swap, redo 8byte swap, convert Vec<u8> to bitvec
    /// That is ALOT of unneeded swapping and conversion
    /// 
    /// New plan: 
    /// - swaping until we get to the EC part
    /// - subset the bitstream to bits[ec_start...]
    /// - just switch the postitions of every pair of 4bytes: 5678 1234 - > 1234 5678  (this is exactly what swap_endian8_swap_endian4() is doing)
    /// - parse the EC
    /// - undo the swtiching for the remainder bits[countstart...] by running swap_endian8_swap_endian4 again
    ///
    /// BIG ISSUE:  The EC might have an uneven number of u32, hence it does NOT align with a u64 border
    ///       ----------------------------u64-----------------------------|--------------------------u64------------
    /// `aaaaaaaabbbbbbbbccccccccdddddddd|eeeeeeeeffffffffgggggggghhhhhhhh|iiiiiiiijjjjjjjjkkkkkkkkllllllllmmmmmmmmnnnnnnnnoooooooopppppppp`
    /// `EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE|CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC|CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC`  // E= EC, C=COUNT
    /// 
    /// Solution: we just check how many bits processed after EC:
    /// - if n*64 : simple undo the swap_endian8_swap_endian4 on the remaining buffer
    /// - if n*32 : the first 4 bytes are already in the correct place and order; just swap anything after the first 4 bytes via swap_endian8_swap_endian4
    /// 
    /// 
    /// last thought: to avoid conversion we could run this endianess-swap in the bistream itself
    ///
    /// =======================
    /// parse the flags
    /// since we already have the correct buffer alignment in count_buffer, we could just do decodeing of flags right there
    /// (cant really store the modified buffer in the struct to pass it to parse_flag)
    /// ACTUALLY, since flag is usually all 0, the entire flags are compressed into a single 0x runlength
    /// and will only occupy a single u64! Just doing the swap thing will be very quick
    // ======================= 

    fn parse_cb(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Cb);
        assert_eq!(self.pos, 0, "must be at beginning of buffer to parse CBs");

        //decode the CBs
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb

        let cb_buffer = &self.buffer[self.pos..];
        let mut fibdec = Self::fibonacci_factory(cb_buffer);

        const CB_RLE_VAL:u64 = 0 + 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut cb_delta_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        while cb_delta_encoded.len() < self.n_elements {
            let item = fibdec.next().unwrap();
            if item == CB_RLE_VAL {
                let runlength = fibdec.next().unwrap() as usize;
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

        // need to remove the trailing 0 which come from padding to the next mutiple of u64 in th CB block
        // we know that the CB block has a length of n x 64
        let bits_processed = fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;

        assert!(
            !self.buffer[self.pos + bits_processed..self.pos + bits_processed + zeros_toremoved].any()
        );

        // adcance the pointer to the start of the next section, i.e. the umis
        self.pos = self.pos + bits_processed + zeros_toremoved;
        self.state = BuszBlockState::Umi;

        cb_delta_encoded
    }

    fn parse_umi(&mut self, cbs: &[u64]) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Umi);

        let umi_buffer = &self.buffer[self.pos..];
        let mut fibdec = Self::fibonacci_factory(umi_buffer);

        const UMI_RLE_VAL: u64 = 0 + 1 ; // since shifted

        // guarantee that last_barcode != rows[0].barcode
        let mut last_cb = cbs[0] + 1;
        let mut umi =0_u64;
        let mut umis: Vec<u64> = Vec::with_capacity(self.n_elements);
        while umis.len() < self.n_elements {
            let diff = fibdec.next().unwrap() - 1;

            let current_cb =  cbs[umis.len()];
            if last_cb !=current_cb {
                umi=0;
            }

            if diff == UMI_RLE_VAL - 1 {
                let runlength = fibdec.next().unwrap() as usize; 
                // adding the same element repeatedly
                umis.resize(umis.len() + runlength, umi - 1);
            } else {
                umi+= diff;
                // decoding single element
                umis.push(umi - 1);//due to shift+1
            }
            last_cb = current_cb;
        }

        // need to remove the trailing 0 which come from padding to the next mutiple of u64
        // we know that the CB block has a length of n x 64
        let bits_processed =  fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;

        assert!(
            !self.buffer[self.pos + bits_processed..self.pos + bits_processed + zeros_toremoved].any()
        );

        // adcance the pointer to the start of the next section, i.e. the umis
        self.pos = self.pos + bits_processed + zeros_toremoved;
        self.state = BuszBlockState::Ec;

        umis

    }

    fn parse_ec(&mut self) -> Vec<u64> {
        // here it gets tricky: previously we swapped the entire stream
        // to u64-little endian
        //
        // however, ECs uses u32 as a basis of operatiors
        // and in order to get the bits appear in the right order we need u32-little endian
        // 1. revert the u64 endian swap
        // 2. swap to u32 little endian
        // 3. do decoding
        // 4 swap back to u64-endian

        let ec_buffer = &self.buffer[self.pos..];

        // println!("converting EC bitstream");

        let bytes = bitslice_to_bytes(ec_buffer);       
        let little_endian_32_bytes = swap_endian8_swap_endian4(&bytes);  //TODO performance: could be done in-place
        let remainder_little_endian_32: &bv::BitSlice<u8, bv::Msb0> =  bv::BitSlice::from_slice(&little_endian_32_bytes);

        // IMPORTANT NOTE: the next line copies self.buffer
        // which is subsequently swapped!
        // if WE DONT COPY, the entire self.buffer would be swapped
        // let mut bytes = bitslice_to_bytes(ec_buffer);       
        // swap_endian8_swap_endian4_inplace(&mut bytes);
        // let remainder_little_endian_32: &bv::BitSlice<u8, bv::Msb0> =  bv::BitSlice::from_slice(&bytes);


        let (ecs, bits_consumed) = newpfd::newpfd_bitvec::decode(remainder_little_endian_32, self.n_elements, PFD_BLOCKSIZE);

        self.pos += bits_consumed;
        self.state = BuszBlockState::Count;

        ecs
    }

    fn parse_counts(&mut self) -> Vec<u64> {
        // ===============================
        // count decoding
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
        // println!("Decoding Counts {:?}", remainder_little_endian_64);
        assert_eq!(self.state, BuszBlockState::Count);

        // things get complicated here: in EC-we had to swap_endian and accomodate u32 encoding
        // this was done within `parse_ec` and did not affect the self.buffer
        // Now we're back to u64. However somethings off. In particular, trying to same procedure 
        // as with the old `decompress_busz_block` (undoing the u32 swap from ec) doesnt work
        // Its an alignment issue!! The EC block might end in the middle of a u64:
        //       ----------------------------u64-----------------------------|--------------------------u64------------
        // `aaaaaaaabbbbbbbbccccccccdddddddd|eeeeeeeeffffffffgggggggghhhhhhhh|iiiiiiiijjjjjjjjkkkkkkkkllllllllmmmmmmmmnnnnnnnnoooooooopppppppp`
        // `EEEEEEEEEEEEEEEEEEEEEEEEEEEEEEEE|CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC|CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC`  // E= EC, C=COUNT
        // 
        //initially stuf has been swapped around using the whole buffer
        // now, trying to undo the swapping on only the COUNT-part wont work
        // The original byte order was `hgfedcba|ponmlkji` -> endian-swap-64 -> abcdefgh|ijklmnop
        // undoing that swap on the count part only: efgh|ijklmnop (without taking into accont the alignment) will lead to lkjihgfe|....ponm|
        //
        // Solution: operate swaping on the entire buffer, THEN select the COUNTs part


        // this does the trick, but copies things (to_bitvec)
        // kind of tricky to not copy in the if, only copy in the else (due to lifetimes)
        
        // using CoW for lazy eval: if aligned with u64, just use the existing buffer
        // if not, do the above shifting, whihc requires allocation
        let _count_buffer = if self.pos % 64 == 0 {
            self.buffer[self.pos..].to_bitvec()
            // Cow::Borrowed(&self.buffer[self.pos..])
        } else {

            // note that we acutally have to backtrack a bit in the original buffer
            // self.buffer looks like this   (ABCDEF... being the actual bytes of the count buffer)
            // EFGH1234| QRSTABCD|YZABMNOP|....UVWX
            //     ^ pos
            // due to the swapping inside EC, the correct EC buffer was processed, and the position moved correctly
            // in the original buffer however, the some of count-buffer bytes ARE BEFORE the current cursor
            //
            // let move it back 4 bytes
            // EFGH1234| QRSTABCD|YZABMNOP|....UVWX
            // ^ pos
            // swap:
            // 1234EFGH|ABCDQRST|MNOPYZAB|UVWX....
            // ^ pos
            //
            // remove the first 4
            // 1234EFGH|ABCDQRST|MNOPYZAB|UVWX....
            //     ^ pos
            // and swap again
            // 1234|ABCDEFGH|MNOPQRST|UVWXYZAB....
            //     ^ pos
            let tmp = &self.buffer[self.pos-32..];
            let bytes_behind_ec = bitslice_to_bytes(tmp);
            let swapped = swap_endian8_swap_endian4(&bytes_behind_ec);
            let swapped_trunc = &swapped[4..];
            let swapped_final = swap_endian8_swap_endian4(swapped_trunc);
            // bv::BitVec::from_slice(&swapped_final)

            // println!("triggered count");
            bv::BitVec::from_vec(swapped_final)
            // Cow::Owned(bv::BitVec::from_vec(swapped_final))
        };
        // let count_buffer = _count_buffer.as_bitslice();
        let count_buffer = _count_buffer.as_ref();
        
        // finally decoding the COUNT bytes
        let mut fibdec = Self::fibonacci_factory(count_buffer);

        const COUNT_RLE_VAL: u64 = 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut counts_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        while counts_encoded.len() < self.n_elements {
            let item = fibdec.next().unwrap();

            if item == COUNT_RLE_VAL {
                let runlength = fibdec.next().unwrap() as usize;
                // just append the RLE repeatedly:
                counts_encoded.resize(counts_encoded.len() + runlength, COUNT_RLE_VAL);
            } else {
                //decode single element
                counts_encoded.push(item);//due to shift+1
            }
        }

        let bits_processed =  fibdec.get_bits_processed();       
        let zeros_toremoved = calc_n_trailing_bits(bits_processed);

        assert!(
            !count_buffer[bits_processed..bits_processed + zeros_toremoved].any(),
            "expected zeros only in this buffer"
        );

        // adcance the pointer to the start of the next section, i.e. the umis
        self.pos = self.pos + bits_processed + zeros_toremoved;
        self.state = BuszBlockState::Flag;

        counts_encoded
    }

    fn parse_flags(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Flag);

        // same issue as in count: some misalginemtn and endianess issues. 
        // do the same thing as in counts
        // if the alignment is off, fix and swap

 
        let _flag_buffer = if self.pos % 64 == 0 {
            self.buffer[self.pos..].to_bitvec()
        } else {

            // note that we acutally have to backtrack a bit in the original buffer
            // self.buffer looks like this   (ABCDEF... being the actual bytes of the count buffer)
            // EFGH1234| QRSTABCD|YZABMNOP|....UVWX
            //     ^ pos
            // due to the swapping inside EC, the correct EC buffer was processed, and the position moved correctly
            // in the original buffer however, the some of count-buffer bytes ARE BEFORE the current cursor
            //
            // let move it back 4 bytes
            // EFGH1234| QRSTABCD|YZABMNOP|....UVWX
            // ^ pos
            // swap:
            // 1234EFGH|ABCDQRST|MNOPYZAB|UVWX....
            // ^ pos
            //
            // remove the first 4
            // 1234EFGH|ABCDQRST|MNOPYZAB|UVWX....
            //     ^ pos
            // and swap again
            // 1234|ABCDEFGH|MNOPQRST|UVWXYZAB....
            //     ^ pos
            let tmp = &self.buffer[self.pos-32..];
            let bytes_behind_ec = bitslice_to_bytes(tmp);
            // tmp.load()

            let swapped = swap_endian8_swap_endian4(&bytes_behind_ec);
            let swapped_trunc = &swapped[4..];
            let swapped_final = swap_endian8_swap_endian4(swapped_trunc);
            // bv::BitVec::from_slice(&swapped_final)

            // println!("triggered flag");
            bv::BitVec::from_vec(swapped_final)
        };
        let flag_buffer = _flag_buffer.as_bitslice();

        // decode flags
        let mut fibdec = Self::fibonacci_factory(flag_buffer);

        const FLAG_RLE_VAL: u64 = 0+1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut flag_decoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        while flag_decoded.len() < self.n_elements {

            let item = fibdec.next().unwrap();
            if item == FLAG_RLE_VAL {
                let runlength = fibdec.next().unwrap() as usize;
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

        let bits_processed = fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed);
        
        assert!(
            !flag_buffer[bits_processed..bits_processed + zeros_toremoved].any()
        );
        // adcance the pointer to the start of the next section, i.e. the umis
        self.pos = self.pos + bits_processed + zeros_toremoved;
        self.state = BuszBlockState::Finished;
        
        assert_eq!(self.pos, self.buffer.len(), "still leftover bits in the buffer!"); 

        flag_decoded
    }


    /// parses the raw bytes in this block
    /// into a list of busrecords
    fn parse_block(&mut self) -> Vec<BusRecord>{
        if self.debug {
            println!("Block-bits:\n{}", bitstream_to_string(&self.buffer[self.pos..]));
            println!("CB pos {}", self.pos);
            }
        let cbs = self.parse_cb();

        if self.debug {
            println!("CBs: {:?}", cbs);

            println!("UMI pos {}", self.pos);
            }
        let umis = self.parse_umi(&cbs);

        if self.debug {
            println!("UMIs: {:?}", umis);

            println!("EC pos {}", self.pos);
            }
        let ec = self.parse_ec();

        if self.debug {
            println!("ECs: {:?}", ec);

            println!("COUNT pos {}", self.pos);
            println!("COUNT buf\n{}", bitstream_to_string(&self.buffer[self.pos..]));
            }
        let count = self.parse_counts();

        if self.debug {
            println!("count: {:?}", count);

            println!("FLAG pos {}", self.pos);
            }
        let flag = self.parse_flags();
        if self.debug {
            println!("flag: {:?}", flag);
        }

        let mut decoded_records = Vec::with_capacity(self.n_elements);
        for (cb,umi,ec,count,flag) in izip!(cbs, umis, ec, count, flag) {
            decoded_records.push(
                BusRecord {
                    CB: cb, 
                    UMI: umi, 
                    EC:ec.try_into().unwrap(), 
                    COUNT:count.try_into().unwrap(), 
                    FLAG:flag.try_into().unwrap()
                }
            );
        }
        decoded_records
    }
}

/// Decompress the `input` busz file into a plain busfile, `output`
pub fn decompress_busfile(input: &str, output: &str) {
    let reader = BuszReader::new(input);
    let mut writer = BusWriterPlain::new(
        output,
        reader.params.clone()
    );

    for r in reader {
        writer.write_record(&r);
    }
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

}