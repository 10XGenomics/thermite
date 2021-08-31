use anyhow::Result;

use needletail::*;

use noodles::{bam, sam};

use bio::alignment::pairwise::Scoring;
use bio::alignment::sparse::HashMapFx;
use bio::alignment::{Alignment, AlignmentOperation};

use std::cmp;
use std::collections::hash_map::Entry;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufWriter};

use crate::aln_writer::*;
use crate::index::*;
use crate::txome::*;

/// Align reads from fastq files and output the alignments in paf, sam, or bam format.
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

/// Attempt to align a single read and return the alignments that are found.
pub fn align_read(
    index: &Index,
    read: &[u8],
    align_opts: &AlignOpts,
) -> Vec<GenomeAlignment> {
    // always make sure reads are uppercase
    let read = read.to_ascii_uppercase();

    let mems = index.all_smems(&read, align_opts.min_seed_len);

    let mut gx_alns = Vec::with_capacity(4);
    let min_aln_score = cmp::max(
        (align_opts.min_aln_score_percent * (read.len() as f32)) as i32,
        align_opts.min_aln_score,
    );
    let mut max_aln_score = min_aln_score;
    // saturating sub just in case floating point error
    // assumes unit scores
    let mut band_width = read.len().saturating_sub(min_aln_score as usize);
    let mut x_drop = read.len().saturating_sub(min_aln_score as usize);

    let scoring = Scoring::from_scores(-1, -1, 1, -1);
    let mut swg = SwgExtend::new(band_width, scoring);

    for hit in &mems {
        let gx_aln = align_seed_hit(index, read, hit, &mut swg, band_width, x_drop);

        // use the running max alignment score to discard low scoring alignments early
        if tx_aln.score < align_opts.min_aln_score
            || tx_aln.score < min_aln_score
            || tx_aln.score < max_aln_score - (align_opts.multimap_score_range as i32)
        {
            continue;
        }

        // narrow band when better alignments are found
        band_width = cmp::min(
            band_width,
            (read.len() + align_opts.multimap_score_range)
                .saturating_sub(gx_aln.gx_aln.score as usize),
        );
        x_drop = cmp::min(
            x_drop,
            (read.len() + align_opts.multimap_score_range)
                .saturating_sub(gx_aln.gx_aln.score as usize),
        );
        max_aln_score = cmp::max(max_aln_score, gx_aln.gx_aln.score);

        gx_alns.push(gx_aln);
    }

    gx_alns.retain(|gx_aln| {
        gx_aln.gx_aln.score >= max_aln_score - (align_opts.multimap_score_range as i32)
    });

    let gx_alns = filter_overlapping(gx_alns);

    // pick an arbitrary max scoring alignment as the primary alignment
    if let Some(gx_aln) = gx_alns
        .iter_mut()
        .find(|gx_aln| gx_aln.gx_aln.score == max_aln_score)
    {
        gx_aln.primary = true;
    }

    gx_alns
}

/// Align a single seed hit and return a GenomeAlignment.
pub fn align_seed_hit(index: &Index, read: &[u8], hit: &Mem, swg: &mut SwgExtend, band_width: usize, x_drop: i32) -> GenomeAlignment {
    // FMD index cannot find a shorter MEM, so this has to be the match
    if mem.len == read.len() {
        // TODO: short circuit
        // check if mem is fully contained in some exon
        // check if mem intersects gene
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
        let mut a = extend_left_right(ref_seq, &relative_hit, read, swg, band_width, x_drop);
        a.ystart += aln_ref.start_idx;
        a.yend += aln_ref.start_idx;
        a
    };

    // intersect seed hit with exons/transcripts
    // find the one best alignment because they are all from the same gene
    let mut best_tx_aln = None;
    let tx_idxs_iter = index
        .txome()
        .exon_to_tx
        .find(Interval::new(hit.ref_idx..hit.ref_idx + hit.len).unwrap())
        .map(|e| *e.data());
    for &tx_idx in tx_idxs_iter {
        let tx_seed = lift_mem_to_tx(&hit, &index.txome().txs[tx_idx]);
        let tx_aln = extend_left_right(&index.txome().txs[tx_idx].seq, &tx_seed, read, swg, band_width, x_drop);
        if best_tx_aln.is_none() || tx_aln.score > best_tx_aln.unwrap().1.score {
            best_tx_aln = Some((tx_idx, tx_aln));
        }

        // cannot do better than exact match
        let match_score = 1;
        if tx_aln.score >= read.len() * match_score {
            break;
        }
    }

    let ref_name = index.txome().txs[tx_idx].chrom.clone();

    if best_tx_aln.is_some() && best_tx_aln.unwrap().1.score >= gx_aln.score {
        // spliced alignment to txome is better
        let (tx_idx, tx_aln) = best_tx_aln.unwrap();

        // lift to concatenated reference coordinates
        // tx alignment could be either forwards or reversed
        let lifted_aln = lift_tx_to_gx(&tx_aln, &index.txome().txs[tx_idx]);
        // convert concatenated genome coordinates to coordinates within some chromosome
        // make sure gx alignment is always relative to forwards strand
        let gx_aln = concat_to_chr_aln(index, lifted_aln);
        GenomeAlignment {
            gx_aln,
            aln_type: AlnType::Exonic { tx_aln, tx_idx },
            ref_name,
            strand,
            primary: false,
        }
    } else {
        // unspliced alignment is better
        // check if genome alignment intersects any genes
        let gene_idxs = index
            .txome()
            .gene_intervals
            .find(Interval::new(gx_aln.ystart..gx_aln.yend).unwrap()).map(|e| *e.data()).collect::<Vec<_>>();
        let gx_aln = concat_to_chr_aln(index, gx_aln);

        if gene_idxs.is_empty() {
            // intergenic
            GenomeAlignment {
                gx_aln,
                aln_type: AlnType::Intergenic,
                ref_name,
                strand,
                primary: false,
            }
        } else {
            // intronic
            // only choose one gene idx for simplicity
            GenomeAlignment {
                gx_aln,
                aln_type: AlnType::Intronic { gene_idx: gene_idxs[0] },
                ref_name,
                strand,
                primary: false,
            }
        }
    }
}

/// Only pick the max scoring alignment when there multiple alignments overlapping.
pub fn filter_overlapping(mut alns: Vec<GenomeAlignment>) -> Vec<GenomeAlignment> {
    if alns.is_empty() {
        return alns;
    }

    alns.sort_by_key(|a| (&a.ref_name, a.strand, a.gx_aln.ystart));

    let mut max_end = 0;
    let mut res = Vec::new();

    for (i, aln) in alns.into_iter().enumerate() {
        if aln.gx_aln.ystart >= max_end || &aln.ref_name != res.last().unwrap().ref_name || aln.strand != res.last().unwrap().strand {
            max_end = aln.gx_aln.yend;
            res.push(aln);
        } else {
            let curr = res.last().unwrap();
            if aln.gx_aln.score > curr.gx_aln.score {
                *res.last_mut().unwrap() = aln;
            }
            max_end = max_end.max(curr.gx_aln.yend);
        }
    }

    res
}

/// Extend a seed hit left and right using banded SWG alignment.
pub fn extend_left_right(ref_seq: &[u8], hit: Mem, read: &[u8], swg: &mut SwgExtend, band_width: usize, x_drop: i32) -> Alignment {
    let x = &read[hit.query_idx + hit.len..];
    let y = &ref_seq[..hit.ref_idx + hit.len..];
    let right_aln = swg.extend(x, y, band_width, x_drop);

    let x = &read[..hit.query_idx].to_owned().reverse();
    let y = &ref_seq[hit.ref_idx.saturating_sub(read.len() + band_width)..hit.ref_idx].to_owned().reverse();
    let left_aln = swg.extend(&x, &y, band_width, x_drop);

    let y_idx = (hit.ref_idx - left_aln.yend, hit.ref_idx + hit.len + right_aln.yend);
    let x_idx = (hit.query_idx - left_aln.xend, hit.query_idx + hit.len + right_aln.xend);
    // assuming unit scores
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

/// Convert an alignment on the concatenated reference to be within a specific chromosome.
fn concat_to_chr_aln(
    index: &Index,
    aln: Alignment,
) -> Alignment {
    let (aln_ref, _ref_idx) = index.idx_to_ref(aln.ystart);
    // convert alignments from concatenated reference coords to chromosome coords
    if aln_ref.strand {
        Alignment {
            ystart: aln.ystart - aln_ref.start_idx,
            yend: aln.yend - aln_ref.start_idx,
            ylen: aln_ref.len,
            ..aln
        }
    } else {
        // need to convert interval to always be [left, right) regardless of strand
        Alignment {
            ystart: aln_ref.len - (aln.yend - aln_ref.start_idx),
            yend: aln_ref.len - (aln.ystart - aln_ref.start_idx),
            ylen: aln_ref.len,
            operations: aln.operations.into_iter().rev().collect(),
            ..aln
        }
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
    /// Range of alignment scores below the max score to output.
    pub multimap_score_range: usize,
}
