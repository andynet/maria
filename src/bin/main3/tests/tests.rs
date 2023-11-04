use crate::*;

#[test]
// fn test() {
//     let parser: GFAParser<usize, ()> = GFAParser::new();
//     let gfa = parser.parse_file("data/pftag/test.gfa")
//         .expect("Error parsing GFA file.");
//
//     let overlap = 1;
//     let segments: Vec<Vec<u8>> = parse_segments(&gfa.segments);
//     let paths: Vec<Vec<usize>> = parse_paths(&gfa.paths);
//
//     let segment_join = join(&segments);
//     let path_join = join(&paths);
//
//     let sa  = SuffixArray::create(&*segment_join);
//     let lcp = lcp_array(&segment_join, &sa).decompress();
//     let isa = permutation_invert(&sa);
//
//     let id  = get_node_ids(&segment_join, &isa);
//     let pos = get_node_pos(&segment_join, &isa);
//
//     let mut segment_join_repr = segment_join.clone();
//     for c in segment_join_repr.iter_mut() {
//         match *c {
//             0   => { *c = b'$'; },
//             1   => { *c = b'#'; },
//             2.. => { *c -= 2; }
//         }
//     }
//     for i in 0..sa.len() {
//         println!("{}\t{}\t{}\t{}\t{}\t{}", i, sa[i], lcp[i], id[i], pos[i],
//             str::from_utf8(&segment_join_repr[sa[i]..]).unwrap());
//     }
//     println!();
//
//     let len     = get_lengths(&segments);
//     let mut rc_rank = get_right_context_rank(&path_join, segments.len());
//     let mut seq_pos = get_sequence_position(&path_join, &len, overlap);
//
//     for i in 0..segments.len() {
//         let iperm = argsort(&rc_rank[i]);
//         rc_rank[i] = permutation_apply(&iperm, &rc_rank[i]);
//         seq_pos[i] = permutation_apply(&iperm, &seq_pos[i]);
//     }
//
//     for id in 0..segments.len() {
//         println!("{}\t{}\t{:?}\t{:?}", 
//             id, len[id], seq_pos[id], rc_rank[id]
//         );
//     }
//     println!();
//
//     // print_tag_array(&lcp, &id, &pos, overlap, &len, &seq_pos, &rc_rank);
// }

// #[test]
fn test_arbitrary() {
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let arbitrary_graph = parser.parse_file("data/pftag/test_arbitrary.gfa")
        .expect("Error parsing GFA file.");

    // let (seq_pos, graph_pos) = parse_graph(&arbitrary_graph);

    // for i in 0..seq_pos.len() {
    //     println!("{}\t{}", seq_pos[i], graph_pos[i].node_id);
    //}
}

#[test]
fn test_process_graph() {
    let (ps1, pn1, ns1, nn1) = process_graph("data/real/SARS-CoV2.5.gfa");
    let (ps2, pn2, ns2, nn2) = process_graph2("data/real/SARS-CoV2.5.gfa");

    assert_eq!(ps1, ps2); // path_starts
    assert_eq!(pn1, pn2);
    assert_eq!(ns1, ns2);
    assert_eq!(nn1, nn2);
}

