use crate::{
    io::{BusReader, BusRecord, BusWriter},
    iterators::CbUmiGroupIterator,
    merger::MultiIterator,
};
use itertools::Itertools;
use std::collections::{BTreeMap, HashMap};
use tempfile::tempdir;

/// Differences to bustools sort
/// 1. bustools sort DOES MERGE records if they share CB/UMI/EC. This implementation does not!
///    We easily could though: in `sort_into_btree` just aggregate the Vec<BusRecord> values

// sorts/inserts an Iterator over records ito a BTreeMap
fn sort_into_btree<I: Iterator<Item = BusRecord>>(
    iterator: I,
) -> BTreeMap<(u64, u64, u32), Vec<BusRecord>> {
    let mut in_mem_sort: BTreeMap<(u64, u64, u32), Vec<BusRecord>> = BTreeMap::new();

    for record in iterator {
        let rlist = in_mem_sort
            .entry((record.CB, record.UMI, record.EC))
            .or_insert(Vec::new());
        rlist.push(record)
    }
    in_mem_sort
}

///
/// Sort a busfile (CB/UMI/EC) in memory!
/// This gets quite bad for larger files
/// uses BTreeMaps internal sorting
pub fn sort_in_memory(busfile: &str, outfile: &str) {
    let reader = BusReader::new(busfile);
    let header = reader.bus_header.clone();

    let in_mem_sort = sort_into_btree(reader);

    // write out
    let mut writer = BusWriter::new(outfile, header);
    for (_cbumi, recordlist) in in_mem_sort {
        writer.write_records(&recordlist);
    }
}

pub fn sort_on_disk(busfile: &str, outfile: &str, chunksize: usize) {
    let reader = BusReader::new(busfile);
    let header = reader.bus_header.clone();

    let mut current_chunk = 0;
    let mut chunkfiles = Vec::new();

    println!("Sorting chunks");
    let tmpdir = tempdir().unwrap();

    for record_chunk in &reader.chunks(chunksize) {
        // sort the chunk in memory
        let in_mem_sort = sort_into_btree(record_chunk);

        //write current sorted file to disk
        let file_path = tmpdir.path().join(format!("tmp_{}.bus", current_chunk));
        let tmpfilename = file_path.to_str().unwrap().to_string();

        // println!("{}", tmpfilename);

        let mut tmpwriter = BusWriter::new(&tmpfilename, header.clone());

        for (_cbumi, recordlist) in in_mem_sort {
            tmpwriter.write_records(&recordlist);
        }
        chunkfiles.push(tmpfilename);
        current_chunk += 1
    }

    // merge all chunks
    println!("Merging {} chunks", chunkfiles.len());
    let mut writer = BusWriter::new(outfile, header);

    let mut iterator_map = HashMap::new();
    for file in chunkfiles.iter() {
        let iter = BusReader::new(file).groupby_cbumi();
        iterator_map.insert(file.to_string(), iter);
    }

    // each file itself is sorted
    // now we only have to merge them
    // if a single cb/umi is split over multiple records, this will put them back together
    let mi = MultiIterator::new(iterator_map);
    for (_cbumi, record_dict) in mi {
        for (_, rlist) in record_dict {
            writer.write_records(&rlist);
        }
    }
    //tmpfiles get clean up once tmpdir is dropped!
}

#[cfg(test)]
mod test {
    use super::{sort_in_memory, sort_on_disk};
    use crate::{
        io::{setup_busfile, BusHeader, BusReader, BusRecord, BusWriter},
        iterators::CbUmiGroupIterator,
    };

    #[test]
    fn test_sort_in_memory() {
        // this is the correct order here:
        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0, };
        let r2 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0, };
        let r3 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0, };
        let r4 = BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0, };
        let r5 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0, };
        let r6 = BusRecord { CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0, };

        let unsorted_records = vec![
            r6.clone(),
            r4.clone(),
            r1.clone(),
            r2.clone(),
            r5.clone(),
            r3.clone(),
        ];
        let (busname, _dir) = setup_busfile(&unsorted_records);

        let outfile = "/tmp/rustbustools_test_sorted.bus";

        sort_in_memory(&busname, outfile);

        let b = BusReader::new(outfile);
        let v: Vec<BusRecord> = b.collect();

        assert_eq!(v, vec![r1, r2, r3, r4, r5, r6]);
    }

    #[test]
    fn test_sort_on_disk() {
        // lets use chunksize 2 and split records over chunks on purpose

        let r1 = BusRecord { CB: 0, UMI: 1, EC: 0, COUNT: 12, FLAG: 0, };
        let r2 = BusRecord { CB: 0, UMI: 1, EC: 1, COUNT: 2, FLAG: 0, };
        let r3 = BusRecord { CB: 0, UMI: 2, EC: 0, COUNT: 12, FLAG: 0, };
        let r4 = BusRecord { CB: 1, UMI: 1, EC: 1, COUNT: 2, FLAG: 0, };
        let r5 = BusRecord { CB: 1, UMI: 2, EC: 1, COUNT: 2, FLAG: 0, };
        let r6 = BusRecord { CB: 2, UMI: 1, EC: 1, COUNT: 2, FLAG: 0, };
        let r7 = BusRecord { CB: 2, UMI: 1, EC: 0, COUNT: 2, FLAG: 0, };

        let unsorted_records = vec![
            // chunk 1
            r6.clone(),
            r4.clone(),
            // chunk 2
            r1.clone(),
            r7.clone(),
            // chunk 3
            r5.clone(),
            r3.clone(),
            // chunk 4
            r2.clone(),
        ];

        let (busname, _dir) = setup_busfile(&unsorted_records);
        let outfile = "/tmp/rustbustools_test_sorted.bus";

        sort_on_disk(&busname, outfile, 2);

        let b = BusReader::new(outfile);

        // the followug doesnt work: r1 and r2 are both (0,1) and hence their order is arbitray
        // assert_eq!(b.collect(), vec![r1, r2, r3, r4, r5, r6, r7]);

        // instead check the sorting of the file implicitely
        let n: usize = b.groupby_cbumi().map(|(_, rlist)| rlist.len()).sum();
        assert_eq!(n, 7)
    }

    use rand::distributions::{Distribution, Uniform};

    #[test]
    fn test_random_file_sort() {
        let cb_len = 16;
        let umi_len = 12;
        // let n_records = 10_000_000;
        // let chunksize = 1_000_000;

        let n_records = 10_000;
        let chunksize = 1_000;

        let cb_distr = Uniform::from(0..10000);
        let umi_distr = Uniform::from(0..10000);
        let mut rng = rand::thread_rng();

        let outfile = "/tmp/test_bus_sort_random.bus";
        let mut writer = BusWriter::new(outfile, BusHeader::new(cb_len, umi_len, 20));
        for _ in 0..n_records {
            let cb = cb_distr.sample(&mut rng);
            let umi = umi_distr.sample(&mut rng);

            let r = BusRecord {
                CB: cb,
                UMI: umi,
                EC: 0,
                COUNT: 1,
                FLAG: 0,
            };
            writer.write_record(&r);
        }
        drop(writer); //make sure everything is written
                      // sort it
        let sorted_out = "/tmp/test_bus_sort_random_sorted.bus";
        sort_on_disk(&outfile, sorted_out, chunksize);

        // check if sorted

        let b = BusReader::new(sorted_out);
        let n: usize = b.groupby_cbumi().map(|(_, rlist)| rlist.len()).sum();
        assert_eq!(n, n_records)
    }
}
