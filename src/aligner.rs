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

    for query_path in query_paths {
        let mut reader = parse_fastx_file(query_path)?;

        while let Some(record) = reader.next() {
            let record = record?;
            let alns = align_read(index, &record.seq(), align_opts);

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
pub fn align_read(
    index: &Index,
    read: &[u8],
    align_opts: &AlignOpts,
) -> Vec<GenomeAlignment> {
    // always make sure reads are uppercase
    let read = read.to_ascii_uppercase();

    let mems = index.all_smems(&read, align_opts.min_seed_len);

    let mut gx_alns = Vec::with_capacity(8);
    let min_aln_score = cmp::max(
        (align_opts.min_aln_score_percent * (read.len() as f32)) as i32,
        align_opts.min_aln_score,
    );
    let mut max_aln_score = min_aln_score;
    // saturating sub just in case floating point error
    let mut band_width = read.len().saturating_sub(min_aln_score as usize);
    let mut coord_score: HashMapFx<(String, bool, usize), i32> = HashMapFx::default();

    let scoring = Scoring::from_scores(-1, -1, 1, -1).yclip(0).xclip_suffix(0);
    let mut swg = SwgExtend::new(band_width, scoring);

    // longest to shortest total seed hit length
    for hit in &mems {
        let gx_aln = align_seed_hit(index, read, &mut swg, hit,);

        // use the running max alignment score to discard low scoring alignments early
        if tx_aln.score < align_opts.min_aln_score
            || tx_aln.score < min_aln_score
            || tx_aln.score < max_aln_score - (align_opts.multimap_score_range as i32)
        {
            continue;
        }

        // ensure that only the max scoring alignment is kept when there are duplicates
        // TODO: remove overlapping alignments if they are the same type (exonic, intronic,
        // intergenic)
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
            (read.len() + align_opts.multimap_score_range)
                .saturating_sub(gx_aln.gx_aln.score as usize),
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

pub fn align_seed_hit(index: &Index, read: &[u8], swg: &mut SwgExtend, hit: &Mem, min_score: i32, band_width: usize, x_drop: i32) -> GenomeAlignment {
    // FMD index cannot find a shorter MEM, so this has to be the match
    if mem.len == read.len() {
        // TODO: short circuit
        return;
    }

    // extend seed hit in genome
    let aln_ref = index.get_ref(hit.ref_idx).1;
    let ref_seq = &index.seq()[aln_ref.start_idx..aln_ref.end_idx];
    let relative_hit = {
        let mut h = hit.clone();
        h.ref_idx -= aln_ref.start_idx;
    };
    let gx_aln = {
        let mut a = extend_left_right(swg, relative_hit, read, ref_seq, band_width, x_drop);
        a.ystart += aln_ref.start_idx;
        a.yend += aln_ref.start_idx;
        a
    };

    // intersect with exons/transcripts
    let mut best_tx_aln = None;
    let tx_idxs_iter = index
        .txome()
        .exon_to_tx
        .find(Interval::new(hit.ref_idx..hit.ref_idx + hit.len).unwrap())
        .map(|e| *e.data())
        .for_each(|tx_idx| {
            let tx_seed = lift_mem_to_tx(&index.txome().txs[tx_idx], hit);
            let tx_aln = extend_left_right(swg, tx_seed, read, &index.txome().txs[tx_idx].seq, band_width, x_drop);
            if best_tx_aln.is_none() || tx_aln.score > best_tx_aln.unwrap().1.score {
                best_tx_aln = Some((tx_idx, tx_aln));
            }
        };

    if best_tx_aln.is_some() && best_tx_aln.unwrap().1.score >= gx_aln.score {
        // spliced alignment to txome is better
        let (tx_idx, tx_aln) = best_tx_aln.unwrap();
        // lift to specific chromosome alignment
        let gx_aln = tx_to_gx_aln(index, tx_aln, tx_idx);
    } else {
        // unspliced alignment is better
        let gene_idxs = index
            .txome()
            .gene_intervals
            .find(Interval::new(gx_aln.ystart..gx_aln.yend).unwrap()).map(|e| *e.data()).collect::<Vec<_>>();
        let gx_aln = concat_to_chr_aln(index, gx_aln);

        if gene_idxs.is_empty() {
            // intergenic

        } else {
            // intronic

        }
    }
}

pub fn extend_left_right(swg: &mut SwgExtend, hit: Mem, read: &[u8], ref_seq: &[u8], band_width: usize, x_drop: i32) -> Alignment {
    let x = &read[hit.query_idx + hit.len..];
    let y = &ref_seq[..hit.ref_idx + hit.len..];
    let right_aln = swg.extend(x, y, band_width, x_drop);

    let x = &read[..hit.query_idx].to_owned().reverse();
    let y = &ref_seq[hit.ref_idx.saturating_sub(read.len() + band_width)..hit.ref_idx].to_owned().reverse();
    let left_aln = swg.extend(&x, &y, band_width, x_drop);

    let y_idx = (hit.ref_idx - left_aln.yend, hit.ref_idx + hit.len + right_aln.yend);
    let x_idx = (hit.query_idx - left_aln.xend, hit.query_idx + hit.len + right_aln.xend);
    let match_score = 1;
    let score = left_aln.score + (match_score * mem.len) + right_aln.score;
    let ops = left_aln.operations.into_iter().rev()
        .chain((0..hit.len).iter().map(|_| AlignmentOperation::Match))
        .chain(right_aln.into_iter())
        .collect::<Vec<_>>();

    Alignment {
        score,
        ystart: y_idx.0,
        xstart: x_idx.0,
        yend: y_idx.1,
        xend: x_idx.1,
        ylen: ref_seq.len(),
        xlen: read.len(),
        operations: ops,
        mode: AlignmentMode::Custom,
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
