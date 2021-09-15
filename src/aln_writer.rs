use anyhow::Result;

use noodles::{bam, sam};

use bio::alignment::AlignmentOperation;
use bio::alphabets::dna;

use std::convert::TryFrom;
use std::io::Write;
use std::{fmt, str};

use crate::index::{Index, RefType};
use crate::txome::{AlnType, GenomeAlignment};

/// Supported alignment output formats.
#[derive(Copy, Clone, PartialEq)]
pub enum OutputFormat {
    Bam,
    Sam,
    Paf,
}

/// Specific writers for different output formats.
pub enum OutputWriter {
    Bam(bam::Writer<Box<dyn Write>>, sam::Header),
    Sam(sam::Writer<Box<dyn Write>>),
    Paf(Box<dyn Write>),
}

/// A single paf record that corresponds to a single alignment.
#[derive(Clone, PartialEq)]
pub struct PafEntry<'a> {
    pub query_name: &'a [u8],
    pub query_len: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: bool,
    pub target_name: &'a [u8],
    pub target_len: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub num_match: usize,
    pub num_match_gap: usize,
    pub map_qual: u8,
}

impl<'a> PafEntry<'a> {
    /// Create a new paf entry based on a genome alignment.
    pub fn new(
        query_name: &'a [u8],
        query_seq: &[u8],
        aln: &'a GenomeAlignment,
        multimap: usize,
    ) -> Result<Self> {
        let num_match = aln
            .gx_aln
            .operations
            .iter()
            .filter(|op| match op {
                AlignmentOperation::Match => true,
                _ => false,
            })
            .count();
        let num_match_gap = aln
            .gx_aln
            .operations
            .iter()
            .filter(|op| match op {
                AlignmentOperation::Yclip(_) => false,
                _ => true,
            })
            .count();
        Ok(Self {
            query_name,
            query_len: query_seq.len(),
            query_start: aln.gx_aln.xstart,
            query_end: aln.gx_aln.xend,
            strand: aln.strand,
            target_name: aln.ref_name.as_bytes(),
            target_len: aln.gx_aln.ylen,
            target_start: aln.gx_aln.ystart,
            target_end: aln.gx_aln.yend,
            num_match,
            num_match_gap,
            map_qual: multimapq(multimap),
        })
    }
}

impl<'a> fmt::Display for PafEntry<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            String::from_utf8_lossy(self.query_name),
            self.query_len,
            self.query_start,
            self.query_end,
            if self.strand { "+" } else { "-" },
            String::from_utf8_lossy(self.target_name),
            self.target_len,
            self.target_start,
            self.target_end,
            self.num_match,
            self.num_match_gap,
            self.map_qual,
        )
    }
}

/// Write out a paf entry.
pub fn write_paf(writer: &mut dyn Write, paf_entry: &PafEntry) -> Result<()> {
    writeln!(writer, "{}", paf_entry)?;
    Ok(())
}

/// Create a new sam record based on a genome alignment.
pub fn aln_to_sam_record(
    index: &Index,
    query_name: &[u8],
    query_seq: &[u8],
    query_qual: &[u8],
    aln: &GenomeAlignment,
    multimap: usize,
    hit_index: usize,
) -> Result<sam::Record> {
    use sam::record::{
        data::{
            field::{Tag, Value},
            Field,
        },
        Data, Flags,
    };

    let query_seq = if aln.strand {
        query_seq.to_owned()
    } else {
        dna::revcomp(query_seq)
    };
    let query_qual = {
        let mut q = query_qual.to_owned();
        if !aln.strand {
            q.reverse();
        }
        q
    };

    let flags = {
        let mut f = Flags::empty();
        if !aln.strand {
            f |= Flags::REVERSE_COMPLEMENTED;
        }
        if !aln.primary {
            f |= Flags::SECONDARY;
        }
        f
    };

    let mapq = multimapq(multimap);
    let num_mismatch = aln
        .gx_aln
        .operations
        .iter()
        .filter(|op| match op {
            AlignmentOperation::Subst => true,
            _ => false,
        })
        .count();

    let data = {
        let mut d = vec![
            Field::new(Tag::AlignmentScore, Value::Int(aln.gx_aln.score as i64)),
            Field::new(Tag::AlignmentHitCount, Value::Int(multimap as i64)),
            Field::new(Tag::HitIndex, Value::Int(hit_index as i64)),
            Field::new(Tag::Other("nM".to_owned()), Value::Int(num_mismatch as i64)),
        ];
        match aln.aln_type {
            AlnType::Exonic { tx_idx, ref tx_aln } => {
                let tx = &index.txome().txs[tx_idx];
                let gene = &index.txome().genes[tx.gene_idx];
                let tx_val = format!(
                    "{},{}{},{}",
                    tx.id,
                    '+',
                    // 0-based position
                    tx_aln.ystart,
                    to_noodles_cigar(&tx_aln.operations)
                );

                d.push(Field::new(
                    Tag::Other("TX".to_owned()),
                    Value::String(tx_val),
                ));
                d.push(Field::new(
                    Tag::Other("GX".to_owned()),
                    Value::String(gene.id.to_owned()),
                ));
                d.push(Field::new(
                    Tag::Other("GN".to_owned()),
                    Value::String(gene.name.to_owned()),
                ));
                d.push(Field::new(Tag::Other("RE".to_owned()), Value::Char('E')));
            }
            AlnType::Intronic { gene_idx } => {
                let gene = &index.txome().genes[gene_idx];
                d.push(Field::new(
                    Tag::Other("GX".to_owned()),
                    Value::String(gene.id.to_owned()),
                ));
                d.push(Field::new(
                    Tag::Other("GN".to_owned()),
                    Value::String(gene.name.to_owned()),
                ));
                d.push(Field::new(Tag::Other("RE".to_owned()), Value::Char('N')));
            }
            AlnType::Intergenic => {
                d.push(Field::new(Tag::Other("RE".to_owned()), Value::Char('I')));
            }
        }
        Data::try_from(d)?
    };

    let read_name = format_read_name(query_name);
    Ok(sam::Record::builder()
        .set_read_name(read_name.parse()?)
        .set_sequence(format_maybe_empty(&query_seq).parse()?)
        .set_quality_scores(format_maybe_empty(&query_qual).parse()?)
        .set_flags(flags)
        .set_data(data)
        .set_reference_sequence_name(aln.ref_name.parse()?)
        // 1-based position
        .set_position(sam::record::Position::try_from(
            (aln.gx_aln.ystart + 1) as i32,
        )?)
        .set_mapping_quality(sam::record::MappingQuality::from(mapq))
        .set_cigar(to_noodles_cigar(&aln.gx_aln.operations))
        .build()?)
}

/// Create a new sam record for an unmapped read.
pub fn unmapped_sam_record(
    query_name: &[u8],
    query_seq: &[u8],
    query_qual: &[u8],
) -> Result<sam::Record> {
    let read_name = format_read_name(query_name);
    Ok(sam::Record::builder()
        .set_read_name(read_name.parse()?)
        .set_sequence(format_maybe_empty(query_seq).parse()?)
        .set_quality_scores(format_maybe_empty(query_qual).parse()?)
        .set_flags(sam::record::Flags::UNMAPPED)
        .build()?)
}

/// Create a sam header for a certain reference index.
pub fn build_sam_header(index: &Index) -> Result<sam::Header> {
    let sam_refs = index
        .refs()
        .into_iter()
        .filter_map(|r| {
            match r.ref_type {
                RefType::Chr { .. } => {
                    Some((
                        r.name.to_owned(),
                        sam::header::ReferenceSequence::new(r.name.to_owned(), r.len as i32).expect(
                            &format!(
                                "Error in creating a SAM header reference sequence with name \"{}\".",
                                r.name
                            ),
                        ),
                    ))
                },
                _ => None,
            }
        })
        .collect();
    Ok(sam::Header::builder()
        .set_reference_sequences(sam_refs)
        .add_program(sam::header::Program::new("thermite"))
        .build())
}

/// Convert Rust-bio's CIGAR format to Noodles-sam's run-length encoded format.
fn to_noodles_cigar(ops: &[AlignmentOperation]) -> sam::record::Cigar {
    use sam::record::{
        cigar::op::{Kind, Op},
        Cigar,
    };

    fn match_op(op: AlignmentOperation, len: usize) -> Op {
        match op {
            // output 'M' for both match and mismatch
            AlignmentOperation::Match => Op::new(Kind::Match, len as u32),
            AlignmentOperation::Subst => Op::new(Kind::Match, len as u32),
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
        // output 'M' for both match and mismatch
        let op = match op {
            AlignmentOperation::Subst => AlignmentOperation::Match,
            o => o,
        };

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

/// Compute the MAPQ of multimapping reads.
///
/// n == 1 -> 255
/// n == 2 -> 3
/// n == 3 -> 2
/// n == 4 -> 1
/// n >= 5 -> 0
fn multimapq(n: usize) -> u8 {
    if n <= 1 {
        255
    } else if n >= 5 {
        0
    } else {
        (-10.0 * (1.0 - 1.0 / (n as f32)).log10()).round() as u8
    }
}

/// Format a read name by converting it to a string and splitting at
/// the first space char.
fn format_read_name(r: &[u8]) -> &str {
    match r.iter().position(|&c| c == b' ') {
        Some(i) => str::from_utf8(&r[..i]).unwrap(),
        None => str::from_utf8(r).unwrap(),
    }
}

/// Avoid empty strings.
fn format_maybe_empty(s: &[u8]) -> &str {
    if s.is_empty() {
        "*"
    } else {
        str::from_utf8(s).unwrap()
    }
}
