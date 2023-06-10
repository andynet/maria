use crate::grammar::Grammar;

#[test]
fn test() {
    let grammar = Grammar::from_file("data/pftag/test_join.txt.plainslp");

    for i in 0..grammar.len() {
        print!("{}", grammar[i] as char);
    }
    println!();
}
