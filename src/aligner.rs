use anyhow::Result;

use needletail::*;

use noodles::{bam, sam};

use bio::alignment::pairwise::{banded::Aligner, Scoring};
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use bio::alignment::{Alignment, AlignmentOperation};

use std::cmp;
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufWriter};

use crate::aln_writer::*;
use crate::index::*;
use crate::txome::*;

/// Kmers of a transcript sequence.
pub type Kmers<'a> = HashMapFx<&'a [u8], Vec<u32>>;

/// Align reads from a fastq file and output the alignments in paf, sam, or bam format.
pub fn align_reads_from_file(
    index: &Index,
    query_paths: &[String],
    output_path: &str,
    output_fmt: OutputFormat,
    align_opts: &AlignOpts,
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

    let mut tx_kmer_cache = HashMapFx::default();

    for query_path in query_paths {
        let mut reader = parse_fastx_file(query_path)?;

        while let Some(record) = reader.next() {
            let record = record?;
            let alns = align_read(index, &mut tx_kmer_cache, &record.seq(), align_opts);

            if alns.is_empty() {
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

                continue;
            }

            for aln in &alns {
                match writer {
                    OutputWriter::Sam(ref mut w) => {
                        let record = aln_to_sam_record(
                            index,
                            &record.id(),
                            &record.seq(),
                            &record.qual().unwrap(),
                            aln,
                            alns.len(),
                        )?;
                        w.write_record(&record)?;
                    }
                    OutputWriter::Bam(ref mut w, ref header) => {
                        let record = aln_to_sam_record(
                            index,
                            &record.id(),
                            &record.seq(),
                            &record.qual().unwrap(),
                            aln,
                            alns.len(),
                        )?;
                        w.write_sam_record(header.reference_sequences(), &record)?;
                    }
                    OutputWriter::Paf(ref mut w) => {
                        let paf = PafEntry::new(&record.id(), &record.seq(), aln, alns.len())?;
                        write_paf(w, &paf)?;
                    }
                }
            }
        }
    }

    Ok(())
}

/// Attempt to align a single read and return the alignment if it is found.
pub fn align_read<'a>(
    index: &'a Index,
    tx_kmer_cache: &mut HashMapFx<usize, Kmers<'a>>,
    read: &[u8],
    align_opts: &AlignOpts,
) -> Vec<GenomeAlignment> {
    // always make sure reads are uppercase
    let read = read.to_ascii_uppercase();

    let tx_hits_map = index.intersect_transcripts(&read, align_opts.min_seed_len);
    let mut tx_hits = tx_hits_map.values().collect::<Vec<_>>();
    tx_hits.sort_unstable_by_key(|k| k.total_len);

    let mut gx_alns = Vec::with_capacity(8);
    let min_aln_score = (align_opts.min_aln_score_percent * (read.len() as f32)) as i32;
    let mut max_aln_score = min_aln_score;
    let mut band_width = read.len() - (min_aln_score as usize);
    let mut coord_score: HashMapFx<(String, bool, usize), i32> = HashMapFx::default();

    // longest to shortest total seed hit length
    for tx_hit in tx_hits.into_iter().rev() {
        if tx_hit.total_len < align_opts.min_total_hit_len {
            break;
        }

        let tx = &index.txome().txs[tx_hit.tx_idx];
        let kmer_cache = tx_kmer_cache
            .entry(tx_hit.tx_idx)
            .or_insert_with(|| hash_kmers(&tx.seq, align_opts.min_seed_len));

        let tx_aln = {
            // TODO: reuse aligner, just change bandwidth
            // align locally in the transcriptome and allow suffix of the read to be clipped
            // note that the tx sequence can be forwards or revcomp
            let scoring = Scoring::from_scores(-1, -1, 1, -1).yclip(0).xclip_suffix(0);
            let mut aligner = Aligner::with_scoring(scoring, align_opts.min_seed_len, band_width);
            let mut aln = aligner.custom_with_prehash(&read, &tx.seq, kmer_cache);
            aln.operations.retain(|op| match op {
                AlignmentOperation::Yclip(_) => false,
                _ => true,
            });
            aln
        };

        // use the running max alignment score to discard low scoring alignments early
        if tx_aln.score < align_opts.min_aln_score
            || tx_aln.score < max_aln_score - (align_opts.multimap_score_range as i32)
        {
            continue;
        }

        let gx_aln = tx_to_gx_aln(index, tx_aln, tx_hit.tx_idx);

        // ensure that only the max scoring alignment is kept when there are duplicates
        match coord_score.entry((gx_aln.ref_name.clone(), gx_aln.strand, gx_aln.gx_aln.ystart)) {
            Entry::Occupied(mut o) => {
                if gx_aln.gx_aln.score > *o.get() {
                    o.insert(gx_aln.gx_aln.score);
                } else {
                    continue;
                }
            }
            Entry::Vacant(v) => {
                v.insert(gx_aln.gx_aln.score);
            }
        }

        // narrow band when better alignments are found
        band_width = cmp::min(
            band_width,
            read.len() + align_opts.multimap_score_range - (gx_aln.gx_aln.score as usize),
        );
        max_aln_score = cmp::max(max_aln_score, gx_aln.gx_aln.score);
        gx_alns.push(gx_aln);
    }

    gx_alns.retain(|gx_aln| {
        gx_aln.gx_aln.score >= max_aln_score - (align_opts.multimap_score_range as i32)
    });
    // pick an arbitrary max scoring alignment as the primary alignment
    if let Some(gx_aln) = gx_alns
        .iter_mut()
        .find(|gx_aln| gx_aln.gx_aln.score == max_aln_score)
    {
        gx_aln.primary = true;
    }
    gx_alns
}

/// Lift a transcriptome alignment to a genome alignment within a specific chromosome.
fn tx_to_gx_aln(index: &Index, tx_aln: Alignment, tx_idx: usize) -> GenomeAlignment {
    // lift to concatenated reference coordinates
    // tx alignment could be either forwards or reversed
    let lifted_aln = lift_tx_to_gx(&tx_aln, &index.txome().txs[tx_idx]);
    // convert concatenated genome coordinates to coordinates within some chromosome
    // make sure gx alignment is always relative to forwards strand
    lifted_aln_to_gx_aln(index, lifted_aln, tx_idx, tx_aln)
}

/// Convert an alignment that is lifted to the concatenated reference to be within
/// a specific chromosome.
fn lifted_aln_to_gx_aln(
    index: &Index,
    lifted_aln: Alignment,
    tx_idx: usize,
    tx_aln: Alignment,
) -> GenomeAlignment {
    let (aln_ref, _ref_idx) = index.idx_to_ref(lifted_aln.ystart);
    let strand = index.txome().txs[tx_idx].strand;
    // convert alignments from concatenated reference coords to chromosome coords
    let gx_aln = if strand {
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
        strand,
        primary: false,
    }
}

/// Struct to conveniently pass around many alignment parameters.
#[derive(Clone, Debug)]
pub struct AlignOpts {
    /// Min length of a SMEM seed.
    pub min_seed_len: usize,
    /// Min alignment score as a percentage of read length
    pub min_aln_score_percent: f32,
    /// Min alignment score
    pub min_aln_score: i32,
    /// Min total length of all seed hits for a read.
    pub min_total_hit_len: usize,
    /// Range of alignment scores below the max score to output.
    pub multimap_score_range: usize,
}
