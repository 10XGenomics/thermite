use serde::{Deserialize, Serialize};

use bio::alignment::{Alignment, AlignmentOperation};
use bio::data_structures::interval_tree::IntervalTree;

use crate::index::Mem;

/// A transcriptome that holds all the genes and transcripts.
#[derive(Serialize, Deserialize)]
pub struct Txome {
    pub genes: Vec<Gene>,
    pub txs: Vec<Tx>,
    pub gene_intervals: IntervalTree<usize, usize>,
}

/// A single transcript and its associated information.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tx {
    pub id: String,
    pub chrom: String,
    pub strand: bool,
    pub exons: Vec<Exon>,
    pub seq: Vec<u8>,
    pub gene_idx: usize,
}

/// A single gene.
#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

/// A single exon.
#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Exon {
    pub start: usize,
    pub end: usize,
    pub tx_idx: usize,
}

impl Exon {
    /// Calculate the length of the exon.
    pub fn len(&self) -> usize {
        self.end - self.start
    }
}

/// Represents an alignment within the genome.
///
/// This also contains the transcriptome alignment from before this alignment
/// is lifted to genome coordinates.
#[derive(Debug, Clone, PartialEq)]
pub struct GenomeAlignment {
    pub gx_aln: Alignment,
    pub aln_type: AlnType,
    pub ref_name: String,
    pub strand: bool,
    pub primary: bool,
}
// TODO: antisense alignments to genes and transcripts

#[derive(Debug, Clone, PartialEq)]
pub enum AlnType {
    Exonic { tx_aln: Alignment, tx_idx: usize },
    Intronic { gene_idx: usize },
    Intergenic,
}

/// Check if larger interval fully contains the smaller one.
pub fn contains(larger: &(usize, usize), smaller: &(usize, usize)) -> bool {
    smaller.0 >= larger.0 && smaller.1 < larger.1
}

/// Check if two intervals intersect.
pub fn intersect(a: &(usize, usize), b: &(usize, usize)) -> bool {
    (a.0 >= b.0 && a.0 < b.1) || (b.0 >= a.0 && b.0 < a.1)
}

/// Lift a MEM from concatenated genome coordinates to a specific exon in a transcript.
pub fn lift_mem_to_tx(mem: &Mem, tx: &Tx) -> Mem {
    let mut exon_sum = 0;

    for exon in &tx.exons {
        if intersect(
            &(mem.ref_idx, mem.ref_idx + mem.len),
            &(exon.start, exon.end),
        ) {
            let start = mem.ref_idx.saturating_sub(exon.start) + exon_sum;
            let start_offset = exon.start.saturating_sub(mem.ref_idx);
            let end = (mem.ref_idx + mem.len).min(exon.end) - exon.start + exon_sum;
            return Mem {
                ref_idx: start,
                query_idx: mem.query_idx + start_offset,
                len: end - start,
            };
        }
        exon_sum += exon.len();
    }

    unreachable!()
}

/// Lift a transcriptome alignment to a concatenated reference alignment.
///
/// Note that the alignment coordinates are for the concatenated reference
/// with all the chromosomes. This will later need to be converted to be relative
/// to a specific reference (chromosome).
pub fn lift_tx_to_gx(tx_aln: &Alignment, tx: &Tx) -> Alignment {
    let mut aln = tx_aln.clone();
    aln.operations.clear();

    let mut i = tx_aln.ystart;
    let mut op_idx = 0;
    // start position of the current exon in the transcript
    // (prefix sum of all previous exon lengths)
    let mut exon_sum = 0;
    // index of the current exon
    let mut exon_idx = 0;

    // find exon where tx alignment starts
    while exon_sum + tx.exons[exon_idx].len() <= i {
        exon_sum += tx.exons[exon_idx].len();
        exon_idx += 1;
    }

    let diff = i - exon_sum;
    aln.ystart = tx.exons[exon_idx].start + diff;

    while op_idx < tx_aln.operations.len() {
        // TODO: edge case: extra exon will be included when last op is insert and it crosses an exon boundary
        if exon_idx + 1 < tx.exons.len() && exon_sum + tx.exons[exon_idx].len() <= i {
            exon_sum += tx.exons[exon_idx].len();
            exon_idx += 1;
            // intron gap
            // use Y clip because intron variant does not exist
            aln.operations.push(AlignmentOperation::Yclip(
                tx.exons[exon_idx].start - tx.exons[exon_idx - 1].end,
            ));
        }

        match tx_aln.operations[op_idx] {
            AlignmentOperation::Match | AlignmentOperation::Subst | AlignmentOperation::Del => {
                i += 1;
            }
            _ => (),
        }

        aln.operations.push(tx_aln.operations[op_idx]);
        op_idx += 1;
    }

    assert_eq!(i, tx_aln.yend);

    let diff = i - exon_sum;
    aln.yend = tx.exons[exon_idx].start + diff;

    aln
}

#[cfg(test)]
mod test {
    use super::*;

    use bio::alignment::AlignmentMode;

    #[test]
    fn test_lift_mem_to_tx() {
        let exons = vec![
            Exon {
                start: 3,
                end: 6,
                tx_idx: 0,
            },
            Exon {
                start: 10,
                end: 13,
                tx_idx: 0,
            },
        ];
        let tx = Tx {
            id: "".to_owned(),
            chrom: "".to_owned(),
            strand: true,
            exons,
            seq: Vec::new(),
            gene_idx: 0,
        };

        let mem = Mem {
            ref_idx: 4,
            query_idx: 3,
            len: 2,
        };
        let correct_mem = Mem {
            ref_idx: 1,
            query_idx: 3,
            len: 2,
        };
        let res_mem = lift_mem_to_tx(&mem, &tx);
        assert_eq!(res_mem, correct_mem);

        let mem = Mem {
            ref_idx: 9,
            query_idx: 3,
            len: 3,
        };
        let correct_mem = Mem {
            ref_idx: 3,
            query_idx: 4,
            len: 2,
        };
        let res_mem = lift_mem_to_tx(&mem, &tx);
        assert_eq!(res_mem, correct_mem);

        let mem = Mem {
            ref_idx: 12,
            query_idx: 3,
            len: 3,
        };
        let correct_mem = Mem {
            ref_idx: 5,
            query_idx: 3,
            len: 1,
        };
        let res_mem = lift_mem_to_tx(&mem, &tx);
        assert_eq!(res_mem, correct_mem);
    }

    #[test]
    fn test_lift_tx_to_gx() {
        let exons = vec![
            Exon {
                start: 3,
                end: 6,
                tx_idx: 0,
            },
            Exon {
                start: 10,
                end: 13,
                tx_idx: 0,
            },
        ];
        let tx = Tx {
            id: "".to_owned(),
            chrom: "".to_owned(),
            strand: true,
            exons,
            seq: Vec::new(),
            gene_idx: 0,
        };
        let ops = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Subst,
            AlignmentOperation::Ins,
            AlignmentOperation::Del,
        ];
        let aln = Alignment {
            score: 0,
            ystart: 1,
            xstart: 0,
            yend: 4,
            xend: 3,
            ylen: 15,
            xlen: 3,
            operations: ops,
            mode: AlignmentMode::Semiglobal,
        };
        let correct_ops = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Subst,
            AlignmentOperation::Yclip(4),
            AlignmentOperation::Ins,
            AlignmentOperation::Del,
        ];
        let correct_aln = Alignment {
            score: 0,
            ystart: 4,
            xstart: 0,
            yend: 11,
            xend: 3,
            ylen: 15,
            xlen: 3,
            operations: correct_ops,
            mode: AlignmentMode::Semiglobal,
        };
        let gx_aln = lift_tx_to_gx(&aln, &tx);
        assert_eq!(gx_aln, correct_aln);
    }

    #[test]
    fn test_lift_tx_to_gx_insert_end() {
        let exons = vec![Exon {
            start: 3,
            end: 6,
            tx_idx: 0,
        }];
        let tx = Tx {
            id: "".to_owned(),
            chrom: "".to_owned(),
            strand: true,
            exons,
            seq: Vec::new(),
            gene_idx: 0,
        };
        let ops = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Subst,
            AlignmentOperation::Ins,
        ];
        let aln = Alignment {
            score: 0,
            ystart: 1,
            xstart: 0,
            yend: 3,
            xend: 2,
            ylen: 15,
            xlen: 2,
            operations: ops,
            mode: AlignmentMode::Semiglobal,
        };
        let correct_ops = vec![
            AlignmentOperation::Match,
            AlignmentOperation::Subst,
            AlignmentOperation::Ins,
        ];
        let correct_aln = Alignment {
            score: 0,
            ystart: 4,
            xstart: 0,
            yend: 6,
            xend: 2,
            ylen: 15,
            xlen: 2,
            operations: correct_ops,
            mode: AlignmentMode::Semiglobal,
        };
        let gx_aln = lift_tx_to_gx(&aln, &tx);
        assert_eq!(gx_aln, correct_aln);
    }
}
