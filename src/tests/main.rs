use gfa::gfa::GFA;
use crate::*;

fn parse_graph(graph: &GFA<usize, ()>) -> (Vec<usize>, Vec<GraphPos>) {
    let mut len = HashMap::new();
    for seg in &graph.segments { len.insert(seg.name, seg.sequence.len()); }

    let mut result = Vec::new();
    for path in &graph.paths {
        let p = str::from_utf8(&path.segment_names).unwrap();
        let p: Vec<GraphPos> = p.split(',').map(|x| x.parse().unwrap()).collect();
        result.extend_from_slice(&p);
    }

    let mut start = Vec::new();
    let mut s = 0;
    for gp in &result {
        start.push(s);
        let id = gp.id;
        s += *len.get(&id).unwrap();
    }
    start.push(s);

    return (start, result);
}

fn process_graph(filename: &str) -> (Vec<usize>, Vec<String>, Vec<usize>, Vec<GraphPos>) {
    let parser: GFAParser<usize, ()> = GFAParser::new();
    let graph = parser.parse_file(filename)
        .expect("Error parsing GFA file.");

    let (start, graph_pos) = parse_graph(&graph);
    let (path_starts, path_ids) = (
        vec![0, 29850, 59697, 89513, 119327],
        vec![
            "ENA|MW565758|MW565758.1".to_string(),
            "ENA|MW565759|MW565759.1".to_string(),
            "ENA|MW565760|MW565760.1".to_string(),
            "ENA|MW565761|MW565761.1".to_string(),
            "ENA|LR883856|LR883856.1".to_string()
        ]
    );
    return (path_starts, path_ids, start, graph_pos);
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

