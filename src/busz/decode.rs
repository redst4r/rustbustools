use std::{io::{BufReader, SeekFrom, Read, Seek}, collections::VecDeque, fs::File};
use crate::{io::{BusHeader, BusRecord, BusWriter, DEFAULT_BUF_SIZE, BUS_HEADER_SIZE, CUGIterator}, busz::{BUSZ_HEADER_SIZE, utils::{swap_endian, calc_n_trailing_bits, bitstream_to_string}, runlength_codec::RunlengthCodec}};
use bitvec::prelude as bv;
use itertools::izip;
use super::{BuszHeader, CompressedBlockHeader, utils::{bitslice_to_bytes, swap_endian8_swap_endian4}, PFD_BLOCKSIZE};

/// Reading a compressed busfile
/// 
/// The analog of [crate::io::BusReader] for compressed bus files
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
#[derive(Debug)]
pub struct BuszReader {
    bus_header: BusHeader,
    busz_header: BuszHeader,
    reader: BufReader<File>,
    buffer: VecDeque<BusRecord>
}

impl BuszReader {
    /// main constructor for busreader, buffersize is set to best performance
    pub fn new(filename: &str) -> Self {
        BuszReader::new_with_capacity(filename, DEFAULT_BUF_SIZE)
    }

    /// Creates a buffered reader over busfiles, with specific buffersize
    pub fn new_with_capacity(filename: &str, bufsize: usize) -> Self {

        // read the regular busheader
        let bus_header = BusHeader::from_file(filename);
        let mut file_handle = File::open(filename).expect("FAIL");

        // advance the file point 20 bytes (pure header) + header.tlen (the variable part)
        let to_seek = BUS_HEADER_SIZE as u64 + bus_header.get_tlen() as u64;
        let _x = file_handle.seek(SeekFrom::Start(to_seek)).unwrap();

        // now we're at the start of the BusZHeader
        let mut buszheader_bytes = [0_u8; BUSZ_HEADER_SIZE];
        file_handle.read_exact(&mut buszheader_bytes).unwrap();

        let bzheader = BuszHeader::from_bytes(&buszheader_bytes);
        //FH is now at the start of the actual contents
        let buf = BufReader::with_capacity(bufsize, file_handle);

        let buffer = VecDeque::with_capacity(bzheader.block_size as usize);
        BuszReader { bus_header, busz_header:bzheader, reader: buf, buffer }
    }

    /// takes the next 8 bytes (u64) out of the stream and interprets it as a buszheader
    /// if we encounter the zero byte [00000000], this indicates the EOF
    fn load_busz_header(&mut self) -> Option<CompressedBlockHeader>{
        // get block header
        let mut blockheader_bytes = [0_u8;8];
        self.reader.read_exact(&mut blockheader_bytes).unwrap();
        if blockheader_bytes == [0,0,0,0,0,0,0,0] {  // EOF
            return None
        }
        let h = CompressedBlockHeader { header_bytes: u64::from_le_bytes(blockheader_bytes)};
        // println!("H bytes {}, H records {}", H.get_blocksize_and_nrecords().0, H.get_blocksize_and_nrecords().1);
        Some(h)
    }

    /// loads a BusZ-block from the stream
    /// we parse the header, which tells us how many bytes and records the block has
    /// with that we pull that amount of bytes out of the stream and parse them into Records
    fn load_busz_block_faster(&mut self) -> Option<Vec<BusRecord>>{

        let h: Option<CompressedBlockHeader> = self.load_busz_header();
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
        let theblock: &bv::BitVec<u8, bv::Msb0> = &bv::BitVec::from_slice(&bigendian_buf);
        let mut block = BuszBlock::new(theblock.as_bitslice(),nrecords as usize);
        
        let records =block.parse_block();
    
        Some(records)
    }
}

impl Iterator for BuszReader {
    type Item=BusRecord;

    fn next(&mut self) -> Option<Self::Item> {
        if self.buffer.is_empty(){
            // println!("buffer empty, loading new block");
            // pull in another block
            let block =self.load_busz_block_faster(); 
            match block {
                Some(records) => {self.buffer = records.into()},
                None => {return None}
            };
        }
        match self.buffer.pop_front() {
            Some(record) => Some(record),
            None => panic!("cant happen")
        }
    }
}

// To make `grouby_cb()` etc work
impl CUGIterator for BuszReader {}

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
    buffer: &'a bv::BitSlice<u8, bv::Msb0>,
    pos: usize, // where we are currently in the buffer
    n_elements: usize ,// how many busrecords are stored in the block
    state: BuszBlockState,
    debug: bool
}

impl <'a> BuszBlock <'a> {
    // pub fn new(buffer: &'a BitSlice<u8, Msb0>, n_elements: usize) -> Self {
    pub fn new(buffer: &'a bv::BitSlice<u8, bv::Msb0>, n_elements: usize) -> Self {
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
    fn fibonacci_factory(stream: &bv::BitSlice<u8, bv::Msb0>) -> newpfd::fibonacci::FibonacciDecoder{
        newpfd::fibonacci::FibonacciDecoder::new(stream)
    }

    fn parse_cb(&mut self) -> Vec<u64> {
        assert_eq!(self.state, BuszBlockState::Cb);
        assert_eq!(self.pos, 0, "must be at beginning of buffer to parse CBs");

        //decode the CBs
        // note that its RLE encoded, i.e. each fibbonacci number is not really a cb
        // println!("Decoding CBs {:?}", block_body);

        let cb_buffer = &self.buffer[self.pos..];
        // let bits_before = self.buffer.len();

        let mut fibdec = Self::fibonacci_factory(cb_buffer);
        // let mut fibdec = FibbonacciDecoder { bitstream: block_body};

        const CB_RLE_VAL:u64 = 0 + 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut cb_delta_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            // println!("Counter {counter}");
            // println!("Bits {:?}", fibdec.bitstream);
            // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
            let item = fibdec.next().unwrap();
            // println!("Item {item}");
            if item == CB_RLE_VAL {
                let runlength = fibdec.next().unwrap();
                // let (runlength, remainder) = bitvec_single_fibbonacci_decode(&mut block_body);
                // block_body = remainder;  //weird scope thing: if we call the above with return named block_body, it'll be its own var and not progagate
                // println!("Decoding run of  {runlength}");
                // println!("Body after   {:?}", fibdec.bitstream);

                for _ in 0..runlength{
                    cb_delta_encoded.push(CB_RLE_VAL -1 );  //due to shift+1
                    counter+= 1;
                }
            } else {
                // println!("Decoding sigle element");
                cb_delta_encoded.push(item - 1);//due to shift+1
                counter+= 1;
            }
        }
        // undo delta encoding, i,e. sumsum
        cb_delta_encoded.iter_mut().fold(0, |acc, x| {
            *x += acc;
            *x
        });
        // println!("cb_delta_encoded: {:?}", cb_delta_encoded);
        // println!("decoded {:?} CBs", cb_delta_encoded.len());
        // println!("Bitstream after CB {:?}", fibdec.bitstream);

        // need to remove the trailing 0 which come from padding to the next mutiple of u64 in th CB block
        // we know that the CB block has a length of n x 64
        let bits_processed = fibdec.get_bits_processed();

        // let padded_size = round_to_multiple(bits_processed, 64);
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
        // println!("before {bits_before}, after {bits_after}; zeros to remove: {zeros_toremoved}");
        // println!("CBs: {}", cb_delta_encoded.len());
        // println!("Cbs: bits proc {}, bits padded {}", bits_processed, zeros_toremoved);
        // println!("CBs: buffer at {}", self.pos + bits_processed + zeros_toremoved);

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
        let mut counter = 0;

        // guarantee that last_barcode != rows[0].barcode
        let mut last_cb = cbs[0] + 1;
        let mut umi =0_u64;
        let mut umis: Vec<u64> = Vec::with_capacity(self.n_elements);
        while counter < self.n_elements {
            let diff = fibdec.next().unwrap() - 1;
            // println!("Decoded item {diff}");
            let current_cb =  cbs[counter];
            if last_cb !=current_cb {
                umi=0;
            }
            // println!("Item {item}");
            if diff == UMI_RLE_VAL - 1 {
                let runlength = fibdec.next().unwrap();
                // let (runlength, remainder) = bitvec_single_fibbonacci_decode(&mut block_body);
                // block_body = remainder;  //weird scope thing: if we call the above with return named block_body, it'll be its own var and not progagate
                // println!("Decoding run of  {runlength}");
                // println!("Body after   {:?}", fibdec.bitstream);

                for _ in 0..runlength{
                    umis.push(umi -1 );  //due to shift+1
                    counter+= 1;
                }
            } else {
                umi+= diff;
                // println!("Decoding sigle element");
                umis.push(umi - 1);//due to shift+1
                counter+= 1;
            }
            last_cb = current_cb;
        }
        // println!("UMIs {:?}", umis);
        // println!("decoded {:?} UMIs", umis.len());

        // need to remove the trailing 0 which come from padding to the next mutiple of u64
        // we know that the CB block has a length of n x 64
        let bits_processed =  fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
        // println!("UMIs: {}", umis.len());
        // println!("UMIs: bits proc {}, bits padded {}", bits_processed, zeros_toremoved);
        // println!("UMIs: buffer at {}", self.pos + bits_processed + zeros_toremoved);

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
        let little_endian_32_bytes = swap_endian8_swap_endian4(&bytes);

        let remainder_little_endian_32: &bv::BitSlice<u8, bv::Msb0> =  bv::BitSlice::from_slice(&little_endian_32_bytes);

        // println!("len remainder_little_endian_32 {}", remainder_little_endian_32.len());
        // println!("remainder_little_endian_32\n{}", bitstream_to_string(remainder_little_endian_32));
        // println!("main buf: remainder_little_endian_32\n{}", bitstream_to_string(ec_buffer));

        let (ecs, bits_consumed) = newpfd::newpfd_bitvec::decode(remainder_little_endian_32, self.n_elements, PFD_BLOCKSIZE);

        // println!("EC buffer used:\n{}", bitstream_to_string(&remainder_little_endian_32[..bits_consumed]));

        self.pos += bits_consumed;
        // println!("after EC main buffer: {}", bitstream_to_string(&self.buffer[self.pos..]));
        // println!("after EC main buffer: {}", bitstream_to_string(&ec_buffer[bits_consumed..]));

        // println!("ECs: {}", ecs.len());
        // println!("ECs: bits proc {}, bits padded 0", bits_consumed);
        // println!("ECs: buffer at {}", self.pos);


        self.state = BuszBlockState::Count;

        ecs
    }

    fn parse_counts(&mut self) ->Vec<u64> {
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

        // undoing initial flip
        let bytes = bitslice_to_bytes(self.buffer);
        // this is what happens in EC: undo u64-swap, apply u32-swap
        let b = swap_endian8_swap_endian4(&bytes);

        // move the bufferposition after EC
        assert_eq!(self.pos % 8, 0);
        let ec_pos = self.pos / 8;

        let bytes_behind_ec = &b[ec_pos..];
        let d = swap_endian8_swap_endian4(&bytes_behind_ec);

        // println!("bytes\n{:?}", bytes);
        // println!("after bytes\n{:?}", d);
        let count_buffer: &bv::BitSlice<u8, bv::Msb0> =  bv::BitSlice::from_slice(&d);

        // println!("Count buffer, reswapped\n{}", bitstream_to_string(count_buffer));


        // let mut count_buffer = &self.buffer[self.pos..];
        // println!("before transform count_buffer\n{}", bitstream_to_string(count_buffer));

        let mut fibdec = Self::fibonacci_factory(count_buffer);

        const COUNT_RLE_VAL: u64 = 1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut counts_encoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            let item = fibdec.next().unwrap();
            // println!("Item {item}");
            if item == COUNT_RLE_VAL {
                let runlength = fibdec.next().unwrap() as usize;
                // for _ in 0..runlength{
                //     counts_encoded.push(COUNT_RLE_VAL);  //due to shift+1
                //     counter+= 1;
                // }
                counts_encoded.resize(counts_encoded.len() + runlength, COUNT_RLE_VAL);
                counter+=runlength;
            } else {
                counts_encoded.push(item);//due to shift+1
                counter+= 1;
            }
        }

        let bits_processed =  fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;

        // println!("COUNT: {}", counts_encoded.len());
        // println!("COUNT: bits proc {}, bits padded {}", bits_processed, zeros_toremoved);
        // println!("COUNT: buffer at {}", self.pos + bits_processed + zeros_toremoved);

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
        // do the same thing: transformations on the entire buffer, the move position to FLAG section
        // undoing initial flip

        let bytes = bitslice_to_bytes(self.buffer);
        // this is what happens in EC: undo u64-swap, apply u32-swap
        let b = swap_endian8_swap_endian4(&bytes);

        // move the bufferposition after EC
        assert_eq!(self.pos % 8, 0);
        let count_pos = self.pos / 8;

        let bytes_behind_count = &b[count_pos..];
        let d = swap_endian8_swap_endian4(&bytes_behind_count);

        // println!("after bytes\n{:?}", d);
        let flag_buffer: &bv::BitSlice<u8, bv::Msb0> =  bv::BitSlice::from_slice(&d);
        // println!("flag_buffer\n{}", bitstream_to_string(flag_buffer));


        let mut fibdec = Self::fibonacci_factory(flag_buffer);

        const FLAG_RLE_VAL: u64 = 0+1;  //since everhthing is shifted + 1, the RLE element is also +1
        let mut flag_decoded: Vec<u64> = Vec::with_capacity(self.n_elements);
        let mut counter = 0;
        while counter < self.n_elements {
            // println!("Counter {counter}");
            // println!("Bits {:?}", fibdec.bitstream);
            // let (item, mut block_body) = bitvec_single_fibbonacci_decode(block_body);
            let item = fibdec.next().unwrap();
            // println!("Item {item}");
            if item == FLAG_RLE_VAL {
                let runlength = fibdec.next().unwrap() as usize;
                // let (runlength, remainder) = bitvec_single_fibbonacci_decode(&mut block_body);
                // block_body = remainder;  //weird scope thing: if we call the above with return named block_body, it'll be its own var and not progagate
                // println!("Decoding run of  {runlength}");
                // println!("Body after   {:?}", fibdec.bitstream);
                // for _ in 0..runlength{
                //     flag_decoded.push(FLAG_RLE_VAL - 1 );  //due to shift+1
                //     counter+= 1;
                // }
                flag_decoded.resize(
                    flag_decoded.len() + runlength,
                    FLAG_RLE_VAL - 1 //due to shift+1
                );
                counter+=runlength;

            } else {
                // println!("Decoding sigle element {}", item-1 );
                flag_decoded.push(item-1); //due to shift+1
                counter+= 1;
            }
        }
        // println!("Falgs: {:?}", flag_decoded);
        // println!("Reamining: {:?}", fibdec.bitstream);

        let bits_processed = fibdec.get_bits_processed();
        let zeros_toremoved = calc_n_trailing_bits(bits_processed); //padded_size - bits_processed;
        
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
    pub fn parse_block(&mut self) -> Vec<BusRecord>{
        if self.debug {
            println!("Block-bits:\n{}", bitstream_to_string(&self.buffer[self.pos..]));
            println!("CB pos {}", self.pos);
            }
        let cbs = self.parse_cb();

        if self.debug {
            println!("CBs: {:?}", cbs);

            println!("UMI pos {}", self.pos);
            }
        // println!("umi: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let umis = self.parse_umi(&cbs);


        if self.debug {
            println!("UMIs: {:?}", umis);

            println!("EC pos {}", self.pos);
            }
        // println!("ec: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let ec = self.parse_ec();


        if self.debug {
            println!("ECs: {:?}", ec);

            println!("COUNT pos {}", self.pos);
            println!("COUNT buf\n{}", bitstream_to_string(&self.buffer[self.pos..]));
            }
        // println!("count: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let count = self.parse_counts();


        if self.debug {
            println!("count: {:?}", count);

            println!("FLAG pos {}", self.pos);
            }
        // println!("flag: {}", bitstream_to_string(&block.buffer[block.pos..]));
        let flag = self.parse_flags();
        if self.debug {
            println!("flag: {:?}", flag);
        }
        // println!("{:?}", cbs);
        // println!("{:?}", umis);
        // println!("{:?}", ec);
        // println!("{:?}", count);
        // println!("{:?}", flag);

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
    let mut writer = BusWriter::new(
        output,
        reader.bus_header.clone()
    );

    for r in reader {
        writer.write_record(&r);
    }
}