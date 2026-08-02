mod cli;

use cli::{Cli, META, build_pool};

fn main() -> std::process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    rsomics_common::run(&output, META, || {
        let pool = build_pool(&cli.threads)?;
        pool.install(|| cli.execute())
    })
}
