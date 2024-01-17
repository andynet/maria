use clap::Parser;
use cli::Args;

mod cli;

fn main() {
    let _ = Args::parse();
}
