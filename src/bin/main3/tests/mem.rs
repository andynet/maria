use crate::mem::MEMReader;

#[test]
fn mem_iterator() {
    let mems_filename = "data/real/reads_R1.mems";
    let ptrs_filename = "data/real/reads_R1.pointers";

    let mem_reader = MEMReader::new(mems_filename, ptrs_filename);
    for (read_id, mems) in mem_reader {
        println!("{read_id}");
        println!("{mems:?}");
    }
}
