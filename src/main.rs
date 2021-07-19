#![deny(warnings)]

use clap::Clap;

use anyhow::Result;

use bincode;

use std::fs::{self, File};
use std::io;
use std::io::prelude::*;
use std::path::Path;
use std::ffi::OsStr;

use thermite_aligner::{aligner, index};
use thermite_aligner::aln_writer::OutputFormat;

fn main() -> Result<()> {
    let opts = ThermiteOpts::parse();

    if opts.verbose {
        println!("{:?}", opts);
    }

    match opts.subcommand {
        SubCommand::Index(index_opts) => {
            let index = index::Index::create_from_file(&index_opts.reference)?;
            let serialized_index = bincode::serialize(&index)?;
            let mut index_file: Box<dyn Write> = match &index_opts.index[..] {
                "-" => Box::new(io::stdout()),
                _ => Box::new(File::create(&index_opts.index)?),
            };
            index_file.write_all(&serialized_index)?;
        }
        SubCommand::Align(align_opts) => {
            let output_fmt = if align_opts.bam {
                let path = Path::new(&align_opts.output);
                let ext = path.extension();
                if ext == Some(OsStr::new("bam")) { OutputFormat::Bam } else { OutputFormat::Sam }
            } else {
                OutputFormat::Paf
            };

            let serialized_index = fs::read(&align_opts.index)?;
            let index = bincode::deserialize(&serialized_index)?;
            drop(serialized_index);

            aligner::align_reads(
                &index,
                &align_opts.queries,
                &align_opts.output,
                output_fmt,
                align_opts.min_seed_len,
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
    #[clap(short = 'o', long = "output", default_value = "-")]
    pub index: String,
}

#[derive(Debug, Clap)]
pub struct Align {
    /// Path to the reference index
    pub index: String,
    /// Paths to the query fastqs
    #[clap(required = true, min_values = 1)]
    pub queries: Vec<String>,
    /// Path to the output file
    #[clap(short = 'o', long = "output", default_value = "-")]
    pub output: String,
    /// Minimum length of an exact seed match
    #[clap(short = 'k', long, default_value = "30")]
    pub min_seed_len: usize,
    /// Output in SAM or BAM format instead of PAF
    #[clap(short = 'a')]
    pub bam: bool,
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
