use anyhow::Result;

use needletail::*;

use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{FMDIndex, FMIndex};
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};
use bio::io::fasta::IndexedReader;
use bio_types::strand::ReqStrand;

use serde::{Deserialize, Serialize};

use transcriptome;

use std::fs::File;
use std::io::BufReader;

use crate::txome::*;

#[derive(Serialize, Deserialize)]
pub struct Index {
    refs: Vec<Ref>,
    seq: Vec<u8>,
    sa: RawSuffixArray,
    fmd: FMDIndex<BWT, Less, Occ>,
    txome: Txome,
}

impl Index {
    pub fn create_from_files(ref_path: &str, annot_path: &str) -> Result<Self> {
        let mut ref_reader = parse_fastx_file(ref_path)?;
        let mut refs = Vec::new();
        let mut seq = Vec::new();

        while let Some(record) = ref_reader.next() {
            let record = record?;

            seq.extend_from_slice(&record.seq());
            seq.push(b'$');
            refs.push(Ref {
                name: String::from_utf8(record.id().to_vec())?,
                strand: true,
                len: record.seq().len(),
                end_idx: seq.len(),
            });

            let revcomp = dna::revcomp(&record.seq()[..]);
            seq.extend_from_slice(&revcomp);
            seq.push(b'$');
            refs.push(Ref {
                name: String::from_utf8(record.id().to_vec())?,
                strand: false,
                len: record.seq().len(),
                end_idx: seq.len(),
            });
        }

        let sa = suffix_array(&seq);
        let bwt = bwt(&seq, &sa);
        let alpha = dna::n_alphabet();
        let less = less(&bwt, &alpha);
        let occ = Occ::new(&bwt, 3, &alpha);
        let fm = FMIndex::new(bwt, less, occ);
        let fmd = FMDIndex::from(fm);

        let mut ref_fai_reader = IndexedReader::from_file(&ref_path)?;
        let txome =
            transcriptome::Transcriptome::from_reader(BufReader::new(File::open(annot_path)?))
                .expect(&format!(
                    "Error in building transcriptome from annotations {}",
                    annot_path
                ));
        let genes = txome
            .genes
            .iter()
            .map(|g| Gene {
                id: g.id.to_owned(),
                name: g.id.to_owned(),
            })
            .collect::<Vec<_>>();
        let txs = txome
            .transcripts
            .iter()
            .map(|tx| {
                let mut tx_seq = Vec::with_capacity(tx.len() as usize);
                tx.get_sequence(&mut ref_fai_reader, &mut tx_seq)
                    .expect(&format!("Error in reading reference file {}", ref_path));
                let exons = tx
                    .exons
                    .iter()
                    .map(|e| Exon {
                        start: e.start as usize,
                        end: e.end as usize,
                        tx_idx: tx.idx.0 as usize,
                    })
                    .collect::<Vec<_>>();

                Tx {
                    id: tx.id.to_owned(),
                    chrom: tx.chrom.to_owned(),
                    strand: tx.strand == ReqStrand::Forward,
                    exons,
                    seq: tx_seq,
                    gene_idx: tx.gene_idx.0 as usize,
                }
            })
            .collect::<Vec<_>>();
        let txome = Txome { genes, txs };

        Ok(Index {
            refs,
            seq,
            sa,
            fmd,
            txome,
        })
    }

    pub fn longest_smem(&self, query: &[u8], min_mem_len: usize) -> Option<Mem> {
        let mut max_smem = None;
        let intervals = self.fmd.all_smems(query, min_mem_len);

        for interval in intervals {
            let forwards_idxs = interval.0.forward().occ(&self.sa);
            let curr = Mem {
                ref_idx: forwards_idxs[0],
                query_idx: interval.1,
                len: interval.2,
            };
            max_smem = match max_smem {
                None => Some(curr),
                Some(max) if curr.len > max.len => Some(curr),
                _ => max_smem,
            };
        }

        max_smem
    }

    pub fn idx_to_ref(&self, idx: usize) -> (&Ref, usize) {
        let ref_idx = self.refs.partition_point(|x| x.end_idx <= idx);
        (
            &self.refs[ref_idx],
            idx - if ref_idx == 0 {
                0
            } else {
                self.refs[ref_idx - 1].end_idx
            },
        )
    }

    pub fn refs(&self) -> &[Ref] {
        &self.refs
    }
}

#[derive(Copy, Clone, PartialEq)]
pub struct Mem {
    pub ref_idx: usize,
    pub query_idx: usize,
    pub len: usize,
}

#[derive(Clone, PartialEq, Serialize, Deserialize)]
pub struct Ref {
    pub name: String,
    pub strand: bool,
    pub len: usize,
    pub end_idx: usize,
}
