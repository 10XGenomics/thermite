use anyhow::Result;

use needletail::*;

use noodles::{bam, sam};

use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use bio::alignment::Alignment;

use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufWriter};

use crate::aln_writer::*;
use crate::index::*;
use crate::txome::*;

type Kmers<'a> = HashMapFx<&'a [u8], Vec<u32>>;

pub fn align_reads_from_file(
    index: &Index,
    query_paths: &[String],
    output_path: &str,
    output_fmt: OutputFormat,
    min_seed_len: usize,
    min_aln_score: i32,
    min_seed_hit_len: usize,
) -> Result<()> {
    let output_writer: Box<dyn Write> = match output_path {
        "-" => Box::new(BufWriter::new(io::stdout())),
        _ => Box::new(BufWriter::new(File::create(output_path)?)),
    };

    let mut writer = match output_fmt {
        OutputFormat::Sam => {
            let sam_header = build_sam_header(index)?;
            let mut w = sam::Writer::new(output_writer);
            w.write_header(&sam_header)?;
            OutputWriter::Sam(w)
        }
        OutputFormat::Bam => {
            let sam_header = build_sam_header(index)?;
            let mut w = bam::Writer::new(output_writer);
            w.write_header(&sam_header)?;
            w.write_reference_sequences(sam_header.reference_sequences())?;
            OutputWriter::Bam(w, sam_header)
        }
        OutputFormat::Paf => OutputWriter::Paf(output_writer),
    };

    let mut tx_kmer_cache = vec![None; index.txome().txs.len()];

    for query_path in query_paths {
        let mut reader = parse_fastx_file(query_path)?;

        while let Some(record) = reader.next() {
            let record = record?;

            if let Some(aln) = align_read(
                index,
                &mut tx_kmer_cache,
                &record.seq(),
                min_seed_len,
                min_aln_score,
                min_seed_hit_len,
            ) {
                match writer {
                    OutputWriter::Sam(ref mut w) => {
                        let record = aln_to_sam_record(
                            index,
                            &record.id(),
                            &record.seq(),
                            &record.qual().unwrap(),
                            &aln,
                        )?;
                        w.write_record(&record)?;
                    }
                    OutputWriter::Bam(ref mut w, ref header) => {
                        let record = aln_to_sam_record(
                            index,
                            &record.id(),
                            &record.seq(),
                            &record.qual().unwrap(),
                            &aln,
                        )?;
                        w.write_sam_record(header.reference_sequences(), &record)?;
                    }
                    OutputWriter::Paf(ref mut w) => {
                        let paf = PafEntry::new(&record.id(), &record.seq(), &aln)?;
                        write_paf(w, &paf)?;
                    }
                }
            } else {
                // unmapped read
                match writer {
                    OutputWriter::Sam(ref mut w) => {
                        let record = unmapped_sam_record(
                            &record.id(),
                            &record.seq(),
                            &record.qual().unwrap(),
                        )?;
                        w.write_record(&record)?;
                    }
                    OutputWriter::Bam(ref mut w, ref header) => {
                        let record = unmapped_sam_record(
                            &record.id(),
                            &record.seq(),
                            &record.qual().unwrap(),
                        )?;
                        w.write_sam_record(header.reference_sequences(), &record)?;
                    }
                    _ => (),
                }
            }
        }
    }

    Ok(())
}

pub fn align_read<'a>(
    index: &'a Index,
    tx_kmer_cache: &mut [Option<Kmers<'a>>],
    read: &[u8],
    min_seed_len: usize,
    min_aln_score: i32,
    min_total_hit_len: usize,
) -> Option<GenomeAlignment> {
    // TODO: tighten band width when good alignments are found
    let mut aligner = Aligner::new(0, -1, |a, b| if a == b { 1 } else { -1 }, min_seed_len, 30);
    let tx_hits_map = index.intersect_transcripts(read, min_seed_len);
    let mut tx_hits = tx_hits_map.values().collect::<Vec<_>>();
    tx_hits.sort_unstable_by_key(|k| k.total_len);
    let mut best_aln: Option<GenomeAlignment> = None;

    for tx_hit in tx_hits.into_iter().rev() {
        if tx_hit.total_len < min_total_hit_len {
            break;
        }

        let tx = &index.txome().txs[tx_hit.tx_idx];
        if tx_kmer_cache[tx_hit.tx_idx].is_none() {
            tx_kmer_cache[tx_hit.tx_idx] = Some(hash_kmers(&tx.seq, min_seed_len));
        }

        let tx_aln = aligner.semiglobal_with_prehash(
            read,
            &tx.seq,
            tx_kmer_cache[tx_hit.tx_idx].as_ref().unwrap(),
        );
        let gx_aln = tx_to_gx_aln(index, tx_aln, tx_hit.tx_idx);

        if gx_aln.gx_aln.score < min_aln_score {
            continue;
        }

        if let Some(ref mut a) = best_aln {
            if gx_aln.gx_aln.score > a.gx_aln.score {
                *a = gx_aln;
            }
        } else {
            best_aln = Some(gx_aln);
        }
    }

    best_aln
}

fn tx_to_gx_aln(index: &Index, tx_aln: Alignment, tx_idx: usize) -> GenomeAlignment {
    // lift to concatenated reference coordinates
    let lifted_aln = lift_tx_to_gx(&tx_aln, &index.txome().txs[tx_idx]);
    // convert concatenated genome coordinates to coordinates within some genome
    lifted_aln_to_gx_aln(index, lifted_aln, tx_idx, tx_aln)
}

fn lifted_aln_to_gx_aln(
    index: &Index,
    lifted_aln: Alignment,
    tx_idx: usize,
    tx_aln: Alignment,
) -> GenomeAlignment {
    // TODO: could just look up chromosome by using transcript
    let (aln_ref, _ref_idx) = index.idx_to_ref(lifted_aln.ystart);
    // convert alignments from concatenated reference coords to genome coords
    let gx_aln = if aln_ref.strand {
        Alignment {
            ystart: lifted_aln.ystart - aln_ref.start_idx,
            yend: lifted_aln.yend - aln_ref.start_idx,
            ylen: aln_ref.len,
            ..lifted_aln
        }
    } else {
        // need to convert interval to always be [left, right) regardless of strand
        Alignment {
            ystart: aln_ref.len - (lifted_aln.yend - aln_ref.start_idx),
            yend: aln_ref.len - (lifted_aln.ystart - aln_ref.start_idx),
            ylen: aln_ref.len,
            operations: lifted_aln.operations.into_iter().rev().collect(),
            ..lifted_aln
        }
    };
    GenomeAlignment {
        gx_aln,
        tx_aln,
        tx_idx,
        ref_name: aln_ref.name.to_owned(),
        strand: aln_ref.strand,
    }
}
