use serde::{Deserialize, Serialize};

use bio::alignment::{Alignment, AlignmentOperation};
use bio::data_structures::interval_tree::IntervalTree;

#[derive(Serialize, Deserialize)]
pub struct Txome {
    pub genes: Vec<Gene>,
    pub txs: Vec<Tx>,
    pub exon_to_tx: IntervalTree<usize, usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Tx {
    pub id: String,
    pub chrom: String,
    pub strand: bool,
    pub exons: Vec<Exon>,
    pub seq: Vec<u8>,
    pub gene_idx: usize,
}

#[derive(Debug, Clone, PartialEq, Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

#[derive(Debug, Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Exon {
    pub start: usize,
    pub end: usize,
    pub tx_idx: usize,
}

impl Exon {
    pub fn len(&self) -> usize {
        self.end - self.start
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct TxHit {
    pub tx_idx: usize,
    pub hits: usize,
    pub total_len: usize,
}

#[derive(Debug, Clone)]
pub struct GenomeAlignment {
    pub gx_aln: Alignment,
    pub tx_aln: Alignment,
    pub tx_idx: usize,
    pub ref_name: String,
    pub strand: bool,
}

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
