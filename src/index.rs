use anyhow::Result;

use needletail::*;

use bio::alphabets::dna;
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{FMDIndex, FMIndex};
use bio::data_structures::suffix_array::{suffix_array, RawSuffixArray};

use serde::{Deserialize, Serialize};

#[derive(Serialize, Deserialize)]
pub struct Index {
    refs: Vec<Ref>,
    seq: Vec<u8>,
    sa: RawSuffixArray,
    fmd: FMDIndex<BWT, Less, Occ>,
}

impl Index {
    pub fn create_from_file(path: &str) -> Result<Self> {
        let mut reader = parse_fastx_file(path)?;
        let mut refs = Vec::new();
        let mut seq = Vec::new();

        while let Some(record) = reader.next() {
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

        Ok(Index { refs, seq, sa, fmd })
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
