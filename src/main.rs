mod cli;

use clap::Parser;
use cli::{Cli, META};

fn main() -> std::process::ExitCode {
    let cli = Cli::parse();
    let common = cli.common.clone();
    rsomics_common::run(&common, META, || cli.execute())
}
