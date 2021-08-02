use anyhow::Result;

use needletail::*;

use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{FMDIndex, FMIndex};
use bio::data_structures::interval_tree::IntervalTree;
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};
use bio::io::fasta::IndexedReader;
use bio::utils::Interval;

use bio_types::strand::ReqStrand;

use serde::{Deserialize, Serialize};

use transcriptome;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::str;

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
        let mut refs = Vec::with_capacity(8);
        let mut seq = Vec::with_capacity(1024);

        #[derive(Clone, Hash, PartialEq, Eq)]
        struct NameStrand(String, bool);
        let mut name_to_ref: HashMap<NameStrand, usize> = HashMap::new();

        // concatenate all reference sequences
        while let Some(record) = ref_reader.next() {
            let record = record?;
            let name = str::from_utf8(record.id())?;

            let start_idx = seq.len();
            seq.extend_from_slice(&record.seq());
            seq.push(b'$');
            name_to_ref.insert(NameStrand(name.to_owned(), true), refs.len());
            refs.push(Ref {
                name: name.to_owned(),
                strand: true,
                len: record.seq().len(),
                start_idx,
                end_idx: seq.len(),
            });

            let start_idx = seq.len();
            let revcomp = dna::revcomp(&record.seq()[..]);
            seq.extend_from_slice(&revcomp);
            seq.push(b'$');
            name_to_ref.insert(NameStrand(name.to_owned(), false), refs.len());
            refs.push(Ref {
                name: name.to_owned(),
                strand: false,
                len: record.seq().len(),
                start_idx,
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
        let transcriptome::Transcriptome {
            genes: txome_genes,
            transcripts: txome_txs,
            ..
        } = transcriptome::Transcriptome::from_reader(BufReader::new(File::open(annot_path)?))
            .expect(&format!(
                "Error in building transcriptome from annotations {}",
                annot_path
            ));

        let genes = txome_genes
            .into_iter()
            .map(|g| Gene {
                id: g.id,
                name: g.name,
            })
            .collect::<Vec<_>>();

        let mut exon_to_tx = IntervalTree::new();
        let txs = txome_txs
            .into_iter()
            .map(|tx| {
                let mut tx_seq = Vec::with_capacity(tx.len() as usize);
                tx.get_sequence(&mut ref_fai_reader, &mut tx_seq)
                    .expect(&format!("Error in reading reference file {}", ref_path));

                let strand = tx.strand == ReqStrand::Forward;
                let tx_ref = &refs[name_to_ref[&NameStrand(tx.chrom.clone(), strand)]];

                let exons = tx
                    .exons
                    .iter()
                    .map(|e| {
                        // transform from coordinates within a single reference sequence to
                        // the coordinates in the entire concatenated reference sequence
                        let exon_start = if strand {
                            (e.start as usize) + tx_ref.start_idx
                        } else {
                            tx_ref.end_idx - 1 - (e.end as usize)
                        };
                        let exon_end = if strand {
                            (e.end as usize) + tx_ref.start_idx
                        } else {
                            tx_ref.end_idx - 1 - (e.start as usize)
                        };
                        let exon_tx_idx = tx.idx.0 as usize;

                        exon_to_tx
                            .insert(Interval::new(exon_start..exon_end).unwrap(), exon_tx_idx);

                        Exon {
                            start: e.start as usize,
                            end: e.end as usize,
                            tx_idx: exon_tx_idx,
                        }
                    })
                    .collect::<Vec<_>>();

                Tx {
                    id: tx.id,
                    chrom: tx.chrom,
                    strand,
                    exons,
                    seq: tx_seq,
                    gene_idx: tx.gene_idx.0 as usize,
                }
            })
            .collect::<Vec<_>>();

        let txome = Txome {
            genes,
            txs,
            exon_to_tx,
        };

        Ok(Index {
            refs,
            seq,
            sa,
            fmd,
            txome,
        })
    }

    pub fn intersect_transcripts(&self, query: &[u8], min_seed_len: usize) -> HashMap<usize, TxHit> {
        let mut tx_hits: HashMap<usize, TxHit> = HashMap::with_capacity(8);
        let intervals = self.fmd.all_smems(query, min_seed_len);

        for interval in intervals {
            let forwards_idxs = interval.0.forward().occ(&self.sa);
            let mem_len = interval.2;

            for ref_idx in &forwards_idxs {
                let seed = Interval::new(*ref_idx..*ref_idx + mem_len).unwrap();
                let tx_idxs = self.txome.exon_to_tx.find(seed);

                for tx_idx in tx_idxs {
                    let mut tx_hit = tx_hits
                        .entry(*tx_idx.data())
                        .or_insert_with(|| TxHit { tx_idx: *tx_idx.data(), hits: 0, total_len: 0 });
                    tx_hit.hits += 1;
                    tx_hit.total_len += mem_len;
                }
            }
        }

        tx_hits
    }

    pub fn longest_smem(&self, query: &[u8], min_seed_len: usize) -> Option<Mem> {
        let mut max_smem = None;
        let intervals = self.fmd.all_smems(query, min_seed_len);

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
        (&self.refs[ref_idx], idx - self.refs[ref_idx].start_idx)
    }

    pub fn refs(&self) -> &[Ref] {
        &self.refs
    }

    pub fn txome(&self) -> &Txome {
        &self.txome
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
    pub start_idx: usize,
    pub end_idx: usize,
}
