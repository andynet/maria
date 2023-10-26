use pfg::pf::PFData;

#[test]
fn iter_works() {
    let pfdata = PFData::from_pfgraph("data/pftag/test.gfa");

    for (sa, id, pos) in pfdata.iter() {
        println!("{sa}\t{id}\t{pos}");
    }
}

#[test]
fn can_load_arbitrary_graph() {
    // let pfdata = PFData::from_graph(
    //     "data/pftag/test_arbitrary.gfa",
    //     "data/pftag/triggers.txt"
    // );
}
