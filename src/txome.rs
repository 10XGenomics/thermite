use serde::{Serialize, Deserialize};

#[derive(Serialize, Deserialize)]
pub struct Txome {
    pub genes: Vec<Gene>,
    pub txs: Vec<Tx>,
}

#[derive(Serialize, Deserialize)]
pub struct Tx {
    pub id: String,
    pub chrom: String,
    pub strand: bool,
    pub exons: Vec<Exon>,
    pub seq: Vec<u8>,
    pub gene_idx: usize,
}

#[derive(Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

#[derive(Serialize, Deserialize)]
pub struct Exon {
    pub start: usize,
    pub end: usize,
    pub tx_idx: usize,
}
