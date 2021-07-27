use serde::{Deserialize, Serialize};

use bio::data_structures::interval_tree::IntervalTree;
use bio::utils::Interval;

#[derive(Serialize, Deserialize)]
pub struct Txome {
    pub genes: Vec<Gene>,
    pub txs: Vec<Tx>,
    pub exon_to_tx: IntervalTree<usize, usize>,
}

#[derive(Clone, Serialize, Deserialize)]
pub struct Tx {
    pub id: String,
    pub chrom: String,
    pub strand: bool,
    pub exons: Vec<Exon>,
    pub seq: Vec<u8>,
    pub gene_idx: usize,
}

#[derive(Clone, PartialEq, Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

#[derive(Copy, Clone, PartialEq, Eq, PartialOrd, Ord, Serialize, Deserialize)]
pub struct Exon {
    pub start: usize,
    pub end: usize,
    pub tx_idx: usize,
}

#[derive(Clone, PartialEq)]
pub struct TxHit {
    pub tx_idx: usize,
    pub hits: usize,
    pub interval: Interval,
}
