#![deny(warnings)]

use clap::Clap;

use anyhow::{bail, Context, Result};

use bincode;

use std::fs::{self, File};
use std::io::prelude::*;

use thermite_aligner::{aligner, index};

fn main() -> Result<()> {
    let opts = ThermiteOpts::parse();

    if opts.verbose {
        println!("{:?}", opts);
    }

    match opts.subcommand {
        SubCommand::Index(index_opts) => {
            let index = index::Index::create_from_file(&index_opts.reference)?;
            let serialized_index = bincode::serialize(&index)?;
            let mut index_file = File::create(&index_opts.index)?;
            index_file.write_all(&serialized_index)?;
        }
        SubCommand::Align(align_opts) => {
            let serialized_index = fs::read(&align_opts.index)?;
            let index = bincode::deserialize(&serialized_index)?;
            drop(serialized_index);

            aligner::align_reads(
                &index,
                &align_opts.query,
                &align_opts.output,
                align_opts.min_mem_len,
            )?;
        }
    }

    Ok(())
}

#[derive(Debug, Clap)]
pub struct Index {
    /// Path to the reference fasta
    pub reference: String,
    /// TAI (Thermite Aligner Index) file to write the index
    pub index: String,
}

#[derive(Debug, Clap)]
pub struct Align {
    /// Path to the reference index
    pub index: String,
    /// Path to the query fastq
    pub query: String,
    /// Path to the output PAF file
    pub output: String,
    /// Minimum length of a Maximal Exact Match
    #[clap(short = 'l', long, default_value = "20")]
    pub min_mem_len: usize,
}

#[derive(Debug, Clap)]
pub enum SubCommand {
    /// Index a reference
    #[clap(name = "index")]
    Index(Index),
    /// Align to an index reference
    #[clap(name = "align")]
    Align(Align),
}

#[derive(Debug, Clap)]
#[clap(name = "thermite", about, version)]
pub struct ThermiteOpts {
    /// Verbose output
    #[clap(long, short)]
    pub verbose: bool,
    #[clap(subcommand)]
    pub subcommand: SubCommand,
}
