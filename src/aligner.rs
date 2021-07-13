use anyhow::Result;

use needletail::*;

use std::fs::File;
use std::io::BufWriter;

use crate::aln_writer::*;
use crate::index::*;

pub fn align_reads(
    index: &Index,
    query_path: &str,
    output_path: &str,
    min_mem_len: usize,
) -> Result<()> {
    let mut reader = parse_fastx_file(query_path)?;
    let mut writer = BufWriter::new(File::create(output_path)?);

    while let Some(record) = reader.next() {
        let record = record?;
        if let Some(mem) = index.longest_smem(&record.seq(), min_mem_len) {
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

    Ok(())
}
