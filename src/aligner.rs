use anyhow::Result;

use needletail::*;

use noodles::{bam, sam};

use std::convert::TryFrom;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufWriter};
use std::str;

use crate::aln_writer::*;
use crate::index::*;

pub fn align_reads(
    index: &Index,
    query_paths: &[String],
    output_path: &str,
    output_fmt: OutputFormat,
    min_seed_len: usize,
) -> Result<()> {
    let output_writer: Box<dyn Write> = match output_path {
        "-" => Box::new(BufWriter::new(io::stdout())),
        _ => Box::new(BufWriter::new(File::create(output_path)?)),
    };

    let mut writer = match output_fmt {
        OutputFormat::Bam | OutputFormat::Sam => {
            let sam_refs = index
                .refs()
                .into_iter()
                .map(|r| {
                    (
                        r.name.to_owned(),
                        sam::header::ReferenceSequence::new(r.name.to_owned(), r.len as i32)
                            .expect("Error in creating a SAM header reference sequence."),
                    )
                })
                .collect();
            let sam_header = sam::Header::builder()
                .set_reference_sequences(sam_refs)
                .add_program(sam::header::Program::new("thermite"))
                .build();

            match output_fmt {
                OutputFormat::Bam => {
                    let mut w = bam::Writer::new(output_writer);
                    w.write_header(&sam_header)?;
                    w.write_reference_sequences(sam_header.reference_sequences())?;
                    OutputWriter::Bam(w, sam_header)
                }
                OutputFormat::Sam => {
                    let mut w = sam::Writer::new(output_writer);
                    w.write_header(&sam_header)?;
                    OutputWriter::Sam(w)
                }
                _ => unreachable!(),
            }
        }
        OutputFormat::Paf => OutputWriter::Paf(output_writer),
    };

    for query_path in query_paths {
        let mut reader = parse_fastx_file(query_path)?;

        while let Some(record) = reader.next() {
            let record = record?;

            if let Some(mem) = index.longest_smem(&record.seq(), min_seed_len) {
                let (mem_ref, ref_idx) = index.idx_to_ref(mem.ref_idx);
                let query_start = mem.query_idx;
                let query_end = mem.query_idx + mem.len;
                // need to convert interval to always be [left, right) regardless of strand
                let target_start = if mem_ref.strand {
                    ref_idx
                } else {
                    mem_ref.len - ref_idx - mem.len
                };
                let target_end = if mem_ref.strand {
                    ref_idx + mem.len
                } else {
                    mem_ref.len - ref_idx
                };

                match writer {
                    OutputWriter::Bam(_, _) | OutputWriter::Sam(_) => {
                        let flags = if mem_ref.strand {
                            sam::record::Flags::empty()
                        } else {
                            sam::record::Flags::REVERSE_COMPLEMENTED
                        };
                        let read_name = str::from_utf8(record.id())?;
                        let record = sam::Record::builder()
                            .set_read_name(read_name.parse()?)
                            .set_flags(flags)
                            .set_reference_sequence_name(mem_ref.name.parse()?)
                            // 1-based position
                            .set_position(sam::record::Position::try_from(
                                (target_start + 1) as i32,
                            )?)
                            .set_mapping_quality(sam::record::MappingQuality::from(255))
                            .set_cigar(format!("{}M", mem.len).parse()?)
                            .set_template_length(mem.len as i32)
                            .build()?;

                        match writer {
                            OutputWriter::Bam(ref mut w, ref header) => {
                                w.write_sam_record(header.reference_sequences(), &record)?
                            }
                            OutputWriter::Sam(ref mut w) => w.write_record(&record)?,
                            _ => unreachable!(),
                        }
                    }
                    OutputWriter::Paf(ref mut w) => {
                        let paf = PafEntry {
                            query_name: &record.id(),
                            query_len: record.seq().len(),
                            query_start,
                            query_end,
                            strand: mem_ref.strand,
                            target_name: mem_ref.name.as_bytes(),
                            target_len: mem_ref.len,
                            target_start,
                            target_end,
                            num_match: mem.len,
                            num_match_gap: mem.len,
                            map_qual: 255,
                        };

                        write_paf(w, &paf)?;
                    }
                }
            }
        }
    }

    Ok(())
}
