use anyhow::Result;

use needletail::*;

use noodles::{bam, sam};

use bio::alignment::pairwise::banded::Aligner;
use bio::alignment::sparse::{hash_kmers, HashMapFx};
use bio::alignment::Alignment;
use bio::alignment::AlignmentOperation;

use std::convert::TryFrom;
use std::fs::File;
use std::io::prelude::*;
use std::io::{self, BufWriter};
use std::str;

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
                    OutputWriter::Bam(_, _) | OutputWriter::Sam(_) => {
                        let flags = if aln.strand {
                            sam::record::Flags::empty()
                        } else {
                            sam::record::Flags::REVERSE_COMPLEMENTED
                        };
                        let data = {
                            use sam::record::{
                                data::{
                                    field::{Tag, Value},
                                    Field,
                                },
                                Data,
                            };
                            let tx = &index.txome().txs[aln.tx_idx];
                            let gene = &index.txome().genes[tx.gene_idx];
                            let tx_val = format!(
                                "{},{}{},{}",
                                tx.id,
                                if tx.strand { '+' } else { '-' },
                                aln.tx_aln.ystart,
                                aln.tx_aln.cigar(false)
                            );
                            Data::try_from(vec![
                                Field::new(Tag::Other("TX".to_owned()), Value::String(tx_val)),
                                Field::new(
                                    Tag::Other("GX".to_owned()),
                                    Value::String(gene.id.to_owned()),
                                ),
                                Field::new(
                                    Tag::Other("GN".to_owned()),
                                    Value::String(gene.name.to_owned()),
                                ),
                                Field::new(Tag::Other("RE".to_owned()), Value::Char('E')),
                            ])?
                        };
                        let read_name = str::from_utf8(record.id())?;
                        let record = sam::Record::builder()
                            .set_read_name(read_name.parse()?)
                            .set_flags(flags)
                            .set_data(data)
                            .set_reference_sequence_name(aln.ref_name.parse()?)
                            // 1-based position
                            .set_position(sam::record::Position::try_from(
                                (aln.genome_aln.ystart + 1) as i32,
                            )?)
                            .set_mapping_quality(sam::record::MappingQuality::from(255))
                            .set_cigar(to_noodles_cigar(&aln.genome_aln.operations))
                            .set_template_length(
                                (aln.genome_aln.xend - aln.genome_aln.xstart) as i32,
                            )
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
                        let num_match = aln
                            .genome_aln
                            .operations
                            .iter()
                            .filter(|op| match op {
                                AlignmentOperation::Match => true,
                                _ => false,
                            })
                            .count();
                        let num_match_gap = aln
                            .genome_aln
                            .operations
                            .iter()
                            .filter(|op| match op {
                                AlignmentOperation::Yclip(_) => false,
                                _ => true,
                            })
                            .count();
                        let paf = PafEntry {
                            query_name: &record.id(),
                            query_len: record.seq().len(),
                            query_start: aln.genome_aln.xstart,
                            query_end: aln.genome_aln.xend,
                            strand: aln.strand,
                            target_name: aln.ref_name.as_bytes(),
                            target_len: aln.genome_aln.ylen,
                            target_start: aln.genome_aln.ystart,
                            target_end: aln.genome_aln.yend,
                            num_match,
                            num_match_gap,
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

pub fn align_read<'a>(
    index: &'a Index,
    tx_kmer_cache: &mut [Option<Kmers<'a>>],
    read: &[u8],
    min_seed_len: usize,
    min_aln_score: i32,
    min_total_hit_len: usize,
) -> Option<GenomeAlignment> {
    // TODO: tighten band width when good alignments are found
    let mut aligner = Aligner::new(-2, -1, |a, b| if a == b { 1 } else { -1 }, min_seed_len, 8);
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
        let genome_aln = {
            // lift to concatenated genome coordinates
            let a = lift_tx_to_genome(&tx_aln, &index.txome().txs[tx_hit.tx_idx]);
            // convert concatenated genome coordinates to genome coordinates
            convert_aln_index_to_genome(index, a, tx_hit.tx_idx, tx_aln)
        };

        if let Some(ref mut a) = best_aln {
            if genome_aln.genome_aln.score > min_aln_score
                && genome_aln.genome_aln.score > a.genome_aln.score
            {
                *a = genome_aln;
            }
        } else {
            best_aln = Some(genome_aln);
        }
    }

    best_aln
}

fn convert_aln_index_to_genome(
    index: &Index,
    mut genome_aln: Alignment,
    tx_idx: usize,
    tx_aln: Alignment,
) -> GenomeAlignment {
    // TODO: could just look up chromosome by using transcript
    let (aln_ref, _ref_idx) = index.idx_to_ref(genome_aln.ystart);
    // convert to genome coordinates
    genome_aln.ystart -= aln_ref.start_idx;
    genome_aln.yend -= aln_ref.start_idx;
    genome_aln.ylen = aln_ref.len;
    // need to convert interval to always be [left, right) regardless of strand
    if !aln_ref.strand {
        genome_aln.ystart = aln_ref.len - genome_aln.yend;
        genome_aln.yend = aln_ref.len - genome_aln.ystart;
        genome_aln.operations.reverse();
    }
    GenomeAlignment {
        genome_aln,
        tx_aln,
        tx_idx,
        ref_name: aln_ref.name.to_owned(),
        strand: aln_ref.strand,
    }
}

fn to_noodles_cigar(ops: &[AlignmentOperation]) -> sam::record::Cigar {
    use sam::record::{
        cigar::op::{Kind, Op},
        Cigar,
    };

    fn match_op(op: AlignmentOperation, len: usize) -> Op {
        match op {
            AlignmentOperation::Match => Op::new(Kind::SeqMatch, len as u32),
            AlignmentOperation::Subst => Op::new(Kind::SeqMismatch, len as u32),
            AlignmentOperation::Del => Op::new(Kind::Deletion, len as u32),
            AlignmentOperation::Ins => Op::new(Kind::Insertion, len as u32),
            AlignmentOperation::Xclip(l) => Op::new(Kind::SoftClip, l as u32),
            // Y clip is repurposed to represent introns
            AlignmentOperation::Yclip(l) => Op::new(Kind::Skip, l as u32),
        }
    }

    let mut v = Vec::with_capacity(ops.len() / 16);
    let mut prev = None;
    let mut prev_len = 0;

    for &op in ops {
        if prev.is_none() || prev.unwrap() != op {
            if let Some(p) = prev {
                v.push(match_op(p, prev_len));
            }
            prev = Some(op);
            prev_len = 1;
        } else {
            prev_len += 1;
        }
    }
    if ops.len() > 0 {
        v.push(match_op(prev.unwrap(), prev_len));
    }
    Cigar::from(v)
}
