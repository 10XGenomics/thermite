use anyhow::Result;

use needletail::*;

use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufWriter};

use crate::aln_writer::*;
use crate::index::*;

pub fn align_reads(
    index: &Index,
    query_paths: &[String],
    output_path: &str,
    min_seed_len: usize,
) -> Result<()> {
    let mut writer: Box<dyn Write> = match output_path {
        "-" => Box::new(BufWriter::new(io::stdout())),
        _ => Box::new(BufWriter::new(File::create(output_path)?)),
    };

    for query_path in query_paths {
        let mut reader = parse_fastx_file(query_path)?;

        while let Some(record) = reader.next() {
            let record = record?;

            if let Some(mem) = index.longest_smem(&record.seq(), min_seed_len) {
                let (mem_ref, ref_idx) = index.idx_to_ref(mem.ref_idx);

                let paf = PafEntry {
                    query_name: &record.id(),
                    query_len: record.seq().len(),
                    query_start: mem.query_idx,
                    query_end: mem.query_idx + mem.len - 1,
                    strand: mem_ref.strand,
                    target_name: &mem_ref.name,
                    target_len: mem_ref.len,
                    target_start: if mem_ref.strand {
                        ref_idx
                    } else {
                        mem_ref.len - ref_idx - mem.len
                    },
                    target_end: if mem_ref.strand {
                        ref_idx + mem.len - 1
                    } else {
                        mem_ref.len - 1 - ref_idx
                    },
                    num_match: mem.len,
                    num_match_gap: mem.len,
                    map_qual: 255,
                };

                write_paf(&mut writer, &paf)?;
            }
        }
    }

    Ok(())
}
