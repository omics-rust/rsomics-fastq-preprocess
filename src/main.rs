mod cli;

use cli::{Cli, META};

fn main() -> std::process::ExitCode {
    let cli = rsomics_help::parse::<Cli>();
    let output = cli.output.clone();
    rsomics_common::run(&output, META, || {
        let pool = cli.threads.build()?;
        pool.install(|| cli.execute())
    })
}
