use crate::pf::PFData;

#[test]
fn iter_works() {
    let pfdata = PFData::from_pfgraph("data/pftag/test.gfa");

    for (sa, id, pos) in pfdata.iter() {
        println!("{sa}\t{id}\t{pos}");
    }
}
