#![deny(warnings)]

pub mod aligner;
pub mod aln_writer;
pub mod cr_transcriptome;
pub mod index;
pub mod swg;
pub mod txome;
pub mod wrapper;

pub use crate::cr_transcriptome::parse_gtf;
pub use crate::cr_transcriptome::transcript_index::Gene;
pub use crate::cr_transcriptome::transcript_sequence;
pub use crate::cr_transcriptome::transcriptome;
pub use crate::cr_transcriptome::utils;
