use clap::{Args, Parser, Subcommand};

#[derive(Parser, Debug)]
#[command(version, about, long_about = None)]
#[command(propagate_version = true)]
struct Cli {
    #[command(subcommand)]
    command: Commands,
}

#[derive(Subcommand, Debug)]
enum Commands {
    /// Generate oriented point cloud surfaces from PDB structure files
    Surface(PdbFiles),
    /// Detect interaction sites on protein surface point clouds.
    Sites(SurfaceFiles),
    /// Visualize point cloud surfaces generated from other subcommands.
    Visualize { surface: String },
    /// Export oriented point cloud surfaces.
    Export(SurfaceFiles),
}

#[derive(Args, Debug)]
struct SurfaceFiles {
    /// Oriented point cloud surface files.
    surface: Vec<String>,
}

#[derive(Args, Debug)]
struct PdbFiles {
    /// Protein structure PDB files.
    pdb: Vec<String>,
}

fn main() {
    let _cli = Cli::parse();
}
