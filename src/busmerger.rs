use std::collections::HashMap;
use crate::{io::{BusWriter, BusHeader}, bus_multi::CellUmiIteratorMulti};

pub fn merge_busfiles_on_overlap(busfile1: &str, busfile2: &str, outfile1: &str, outfile2: &str){
    // will extract all busrecords that appear in both inputs and write them to the respective outputs
    // there'll be two output files, each contining the shared reads from the respective input file
    let h: HashMap<String, String> = HashMap::from([
        ("f1".to_string(), busfile1.to_string()), 
        ("f2".to_string(), busfile2.to_string())
    ]);

    let mut writers: HashMap<String, BusWriter> = HashMap::from([
        ("f1".to_string(), BusWriter::new(outfile1, BusHeader::new(16, 12, 20))), 
        ("f2".to_string(), BusWriter::new(outfile2, BusHeader::new(16, 12, 20))), 
    ]);

    let cbumi_merge_iter = CellUmiIteratorMulti::new(&h);
    
    for (_cbumi, record_map) in cbumi_merge_iter{
        // if the CB/UMI is present in both files, write
        if record_map.len() == 2{
            for (name, records) in record_map{
                let w1 = writers.get_mut(&name).unwrap();
                w1.write_records(&records)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use crate::io::{BusRecord, BusIteratorBuffered, setup_busfile};
    use super::*;

    fn get_records(fname: &str) -> Vec<BusRecord>{
        let reader = BusIteratorBuffered::new(fname);
        let records: Vec<BusRecord> = reader.into_iter().collect();
        records
    }

    #[test]
    fn test_merge(){
        let r1 =BusRecord{CB: 0, UMI: 21, EC: 0, COUNT: 2, FLAG: 0};
        let r2 =BusRecord{CB: 1, UMI: 2, EC: 0, COUNT: 12, FLAG: 0}; 
        let r3 =BusRecord{CB: 1, UMI: 3, EC: 0, COUNT:  2, FLAG: 0}; 
        let r4 =BusRecord{CB: 3, UMI: 0, EC: 0, COUNT:  2, FLAG: 0}; 
        let r5 =BusRecord{CB: 3, UMI: 0, EC: 1, COUNT:  2, FLAG: 0}; 

        let v1 = vec![r1.clone(),r2.clone(),r3.clone(), r4.clone(), r5.clone()];

        let s2 = BusRecord{CB: 1, UMI: 2, EC: 1, COUNT: 12, FLAG: 0};
        let s3 = BusRecord{CB: 2, UMI: 3, EC: 1, COUNT:  2, FLAG: 0}; 
        let s4 = BusRecord{CB: 3, UMI: 0, EC: 1, COUNT:  2, FLAG: 0}; 

        let v2 = vec![s2.clone(),s3.clone(), s4.clone()];

        // let input1 = "/tmp/merge1.bus";
        // let input2 = "/tmp/merge2.bus";

        let (input1, _dir1) = setup_busfile(&v1); //input1
        let (input2, _dir2) =  setup_busfile(&v2); // input2

        // TODO dangerous: hardcoded filenames
        let output1 = "/tmp/merge1_out.bus";
        let output2 = "/tmp/merge2_out.bus";

        merge_busfiles_on_overlap(&input1, &input2, output1, output2);

        assert_eq!(get_records(output1), vec![r2, r4,r5]);
        assert_eq!(get_records(output2), vec![s2, s4]);
    }
}