use anyhow::Result;

use noodles::{bam, sam};

use bio::alignment::AlignmentOperation;

use std::convert::TryFrom;
use std::io::Write;
use std::{fmt, str};

use crate::index::Index;
use crate::txome::GenomeAlignment;

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
    pub fn new(query_name: &'a [u8], query_seq: &[u8], aln: &'a GenomeAlignment) -> Result<Self> {
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
            map_qual: 255,
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
) -> Result<sam::Record> {
    use sam::record::{
        data::{
            field::{Tag, Value},
            Field,
        },
        Data, Flags,
    };

    let flags = if aln.strand {
        Flags::empty()
    } else {
        Flags::REVERSE_COMPLEMENTED
    };
    let data = {
        let tx = &index.txome().txs[aln.tx_idx];
        let gene = &index.txome().genes[tx.gene_idx];
        let tx_val = format!(
            "{},{}{},{}",
            tx.id,
            if tx.strand { '+' } else { '-' },
            // 1-based position
            aln.tx_aln.ystart + 1,
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
    let read_name = str::from_utf8(query_name)?;
    Ok(sam::Record::builder()
        .set_read_name(read_name.parse()?)
        .set_sequence(str::from_utf8(query_seq)?.parse()?)
        .set_quality_scores(str::from_utf8(query_qual)?.parse()?)
        .set_flags(flags)
        .set_data(data)
        .set_reference_sequence_name(aln.ref_name.parse()?)
        // 1-based position
        .set_position(sam::record::Position::try_from(
            (aln.gx_aln.ystart + 1) as i32,
        )?)
        .set_mapping_quality(sam::record::MappingQuality::from(255))
        .set_cigar(to_noodles_cigar(&aln.gx_aln.operations))
        .set_template_length((aln.gx_aln.yend - aln.gx_aln.ystart) as i32)
        .build()?)
}

/// Create a new sam record for an unmapped read.
pub fn unmapped_sam_record(
    query_name: &[u8],
    query_seq: &[u8],
    query_qual: &[u8],
) -> Result<sam::Record> {
    let read_name = str::from_utf8(query_name)?;
    Ok(sam::Record::builder()
        .set_read_name(read_name.parse()?)
        .set_sequence(str::from_utf8(query_seq)?.parse()?)
        .set_quality_scores(str::from_utf8(query_qual)?.parse()?)
        .set_flags(sam::record::Flags::UNMAPPED)
        .build()?)
}

/// Create a sam header for a certain reference index.
pub fn build_sam_header(index: &Index) -> Result<sam::Header> {
    let sam_refs = index
        .refs()
        .into_iter()
        .map(|r| {
            (
                r.name.to_owned(),
                sam::header::ReferenceSequence::new(r.name.to_owned(), r.len as i32).expect(
                    &format!(
                        "Error in creating a SAM header reference sequence with name \"{}\".",
                        r.name
                    ),
                ),
            )
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
