#![deny(warnings)]

use clap::Clap;

use anyhow::Result;

use bincode;

use std::ffi::OsStr;
use std::fs::File;
use std::io;
use std::io::prelude::*;
use std::path::Path;

use thermite_aligner::aln_writer::OutputFormat;
use thermite_aligner::{aligner, index};

fn main() -> Result<()> {
    let opts = ThermiteOpts::parse();

    if opts.verbose {
        println!("{:?}", opts);
    }

    match opts.subcommand {
        SubCommand::Index(index_opts) => {
            let index =
                index::Index::create_from_files(&index_opts.reference, &index_opts.annotations)?;
            let index_file: Box<dyn Write> = match &index_opts.index[..] {
                "-" => Box::new(io::stdout()),
                _ => Box::new(File::create(&index_opts.index)?),
            };
            bincode::serialize_into(index_file, &index)?;
        }
        SubCommand::Align(align_opts) => {
            let output_fmt = if align_opts.bam {
                let path = Path::new(&align_opts.output);
                let ext = path.extension();
                if ext == Some(OsStr::new("bam")) {
                    OutputFormat::Bam
                } else {
                    OutputFormat::Sam
                }
            } else {
                OutputFormat::Paf
            };

            let index_file = File::open(&align_opts.index)?;
            let index = bincode::deserialize_from(index_file)?;

            aligner::align_reads_from_file(
                &index,
                &align_opts.queries,
                &align_opts.output,
                output_fmt,
                align_opts.min_seed_len,
                align_opts.min_aln_score,
                align_opts.min_total_hit_len,
            )?;
        }
    }

    Ok(())
}

#[derive(Debug, Clap)]
pub struct Index {
    /// Path to the reference fasta. Thermite expects a fasta index with the same
    /// file name and the extension `.fasta.fai` to exist.
    pub reference: String,
    /// Path to the GTF annotations of the reference fasta
    pub annotations: String,
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
    /// Minimum alignment score
    #[clap(short = 's', long, default_value = "40")]
    pub min_aln_score: i32,
    /// Minimum total seed hit length (sum of the lengths of all seed hits)
    #[clap(long, default_value = "40")]
    pub min_total_hit_len: usize,
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
