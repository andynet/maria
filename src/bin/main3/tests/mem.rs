use crate::mem::MEMReader;

#[test]
fn mem_iterator() {
    let mems_filename = "data/real/reads_R1.mems";
    let ptrs_filename = "data/real/reads_R1.pointers";

    let mem_reader = MEMReader::new(&mems_filename, &ptrs_filename);
    for (read_id, mems) in mem_reader {
        println!("{read_id}");
        println!("{mems:?}");
        for mem in mems {
            // let gp: Vec<GraphPos> = get_graph_positions(&seq, &ref_mem, &stag, &ssa);

            // for gpos in gp {
            //     let (node_id, sign): (String, char) = parse_node_id(gpos.node_id);
            //     println!("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}",
            //         read_id, mem.0, adj, ref_pos, mem.1, node_id, sign, gpos.pos
            //     )
            // }
        }
    }
}
