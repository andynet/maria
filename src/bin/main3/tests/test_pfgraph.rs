use crate::pf::PFData;

#[test]
fn test() {
    let pfdata = PFData::from_PFGraph("data/pftag/test.gfa");

    for (sa, id, pos) in pfdata {
        println!("{sa}\t{id}\t{pos}");
    }
}
