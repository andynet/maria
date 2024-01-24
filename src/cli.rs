use clap::{Parser, Subcommand};

#[rustfmt::skip]
#[test]
fn print_maria_noargs_or_help() {
    let args1 = Args::try_parse_from(["maria"].iter()).err().unwrap();
    let args2 = Args::try_parse_from(["maria", "-h"].iter()).err().unwrap();
    assert_eq!(args1.to_string(), args2.to_string());
    assert_eq!(args1.to_string(), "\
        Usage: \n\
        maria index <graph>.gfa -t <triggers.txt>\n\
        maria align <graph>.gfa <reads>.fastq > <output.gaf> \n\
        \n\
        \n\
        Commands:\n  \
          index  Create a run-length compressed tag array <graph>.tag\n  \
          align  Find all positions of a match in a graph. Matches to the reference can be found by MONI\n  \
          help   Print this message or the help of the given subcommand(s)\n\
        \n\
        Options:\n  \
          -h, --help  Print help\n\
    ");
}

#[rustfmt::skip]
#[test]
fn print_maria_index_noargs_or_help() {
    let args1 = Args::try_parse_from(["maria", "index"].iter()).err().unwrap();
    let args2 = Args::try_parse_from(["maria", "index", "-h"].iter()).err().unwrap();
    assert_eq!(args1.to_string(), args2.to_string());
    assert_eq!(args1.to_string(), "\
        Create a run-length compressed tag array <graph>.tag\n\
        \n\
        Usage: maria index <GFA> -t <TRIGGERS>\n\
        \n\
        Arguments:\n  \
          <GFA>  Graph in GFA format\n\
        \n\
        Options:\n  \
          -t <TRIGGERS>      Trigger file used for prefix-free tag array construction.\n                     \
                             Contains one trigger (e.g. TAA) per line.\n                     \
                             Triggers do not influence the resulting tag array,\n                     \
                             only the time and space complexity of the construction.\n  \
          -h, --help         Print help\n\
    ");
}

#[rustfmt::skip]
#[test]
fn print_maria_align_noargs_or_help() {
    let args1 = Args::try_parse_from(["maria", "align"].iter()).err().unwrap();
    let args2 = Args::try_parse_from(["maria", "align", "-h"].iter()).err().unwrap();
    assert_eq!(args1.to_string(), args2.to_string());
    assert_eq!(args1.to_string(), "\
        Find all positions of a match in a graph. Matches to the reference can be found by MONI\n\
        \n\
        Usage: maria align <GFA> <READS> > output.gaf\n\
        \n\
        Arguments:\n  \
          <GFA>    Graph in GFA format. For <graph>.gfa, tag array <graph>.tag and SLP grammar <graph>.slp need to be present\n  \
          <READS>  File containing reads. For <reads>.fastq, MONI outputs <reads>.mems and <reads>.pointers need to be present\n\
        \n\
        Options:\n  \
          -o <OUTPUT>      Output in GAF format [default: stdout]\n  \
          -h, --help       Print help\n\
    ");
}

#[derive(Parser, Debug)]
#[command(override_usage = "\n\
    maria index <graph>.gfa -t <triggers.txt>\n\
    maria align <graph>.gfa <reads>.fastq > <output.gaf> \n\
")]
pub struct Args {
    #[command(subcommand)]
    pub command: Commands,
}

#[derive(Subcommand, Debug)]
pub enum Commands {
    /// Create a run-length compressed tag array <graph>.tag
    #[command(override_usage = "maria index <GFA> -t <TRIGGERS>", arg_required_else_help = true)]
    Index {
        /// Graph in GFA format
        gfa: String,

        /// Trigger file used for prefix-free tag array construction.
        /// Contains one trigger (e.g. TAA) per line.
        /// Triggers do not influence the resulting tag array,
        /// only the time and space complexity of the construction.
        #[arg(short = 't', verbatim_doc_comment)]
        triggers: String,
    },

    /// Find all positions of a match in a graph.
    /// Matches to the reference can be found by MONI.
    #[command(override_usage = "maria align <GFA> <READS> > output.gaf", arg_required_else_help = true)]
    Align {
        /// Graph in GFA format.
        /// For <graph>.gfa, tag array <graph>.tag and SLP grammar <graph>.slp need to be present.
        gfa: String,

        /// File containing reads.
        /// For <reads>.fastq, MONI outputs <reads>.mems and <reads>.pointers need to be present.
        reads: String,

        /// Output in GAF format [default: stdout]
        #[arg(short = 'o')]
        output: Option<String>,
    },
}
