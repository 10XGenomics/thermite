// allow unused imports since the transcriptome feature may be disabled
#![allow(unused_imports)]

use anyhow::Result;

use needletail::*;

use bio::alphabets::dna;
use bio::data_structures::interval_tree::IntervalTree;
use bio::io::fasta::IndexedReader;
use bio::utils::Interval;
use bio_opt::alphabets::Alphabet;
use bio_opt::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio_opt::data_structures::fmindex::{FMDIndex, FMIndex};
use bio_opt::data_structures::suffix_array::{RawSuffixArray, SampledSuffixArray, SuffixArray};

use bio_types::strand::ReqStrand;

use libdivsufsort_rs::divsufsort64;

use serde::{Deserialize, Serialize};

#[cfg(feature = "transcriptome")]
use transcriptome;

use std::borrow::Cow;
use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::iter::FromIterator;
use std::{mem, str};

use crate::txome::*;

/// An index of a set of reference sequences, along with the transcriptome.
///
/// Multiple different reference sequences (chromosomes) are concatenated internally,
/// so it is important to convert the concatenated coordinates back to coordinates within
/// a specific chromosome when outputting alignments.
#[derive(Serialize, Deserialize)]
pub struct Index {
    refs: Vec<Ref>,
    sa: SampledSuffixArray<BWT, Less, Occ>,
    txome: Txome,
}

impl Index {
    /// Create an index from a reference fasta file and its annotations file.
    ///
    /// The fasta file is expected to be already indexed, so a .fasta.fai file exists
    /// with the same file name.
    #[cfg(feature = "transcriptome")]
    pub fn create_from_files(
        ref_path: &str,
        annot_path: &str,
        sa_sampling_rate: usize,
        occ_sampling_rate: usize,
    ) -> Result<Self> {
        let mut ref_reader = parse_fastx_file(ref_path)?;
        let mut refs = Vec::with_capacity(8);
        let mut seq = Vec::with_capacity(1024);

        let mut name_to_ref: HashMap<NameStrand, usize> = HashMap::new();

        // concatenate all chromosome sequences
        while let Some(record) = ref_reader.next() {
            let record = record?;
            let name = str::from_utf8(record.id())?.split(' ').next().unwrap();

            // end index includes '$', length does not
            let start_idx = seq.len();
            let mut curr_seq = record.seq().to_vec();
            curr_seq.make_ascii_uppercase();
            seq.extend_from_slice(&curr_seq);
            seq.push(b'$');
            name_to_ref.insert(NameStrand(name.to_owned(), true), refs.len());
            refs.push(Ref {
                name: name.to_owned(),
                ref_type: RefType::Chr { seq: Some(curr_seq) },
                strand: true,
                len: record.seq().len(),
                start_idx,
                end_idx: seq.len(),
            });

            let start_idx = seq.len();
            let mut revcomp = dna::revcomp(&record.seq()[..]);
            revcomp.make_ascii_uppercase();
            seq.extend_from_slice(&revcomp);
            seq.push(b'$');
            name_to_ref.insert(NameStrand(name.to_owned(), false), refs.len());
            refs.push(Ref {
                name: name.to_owned(),
                ref_type: RefType::Chr { seq: None },
                strand: false,
                len: record.seq().len(),
                start_idx,
                end_idx: seq.len(),
            });
        }

        let txome = Self::create_txome(ref_path, annot_path, &name_to_ref, &refs)?;

        // concatenate transcript sequences after chromosome sequences
        // TODO: revcomp transcript sequences for antisense alignments?
        for (i, tx) in txome.txs.iter().enumerate() {
            // end index includes '$', length does not
            let start_idx = seq.len();
            seq.extend_from_slice(&tx.seq);
            seq.push(b'$');
            refs.push(Ref {
                name: "".to_owned(),
                ref_type: RefType::Tx { tx_idx: i },
                strand: true,
                len: tx.seq.len(),
                start_idx,
                end_idx: seq.len(),
            });
        }

        // construct FMD index
        let sa =
            cast_vec_i64_to_usize(divsufsort64(&seq).expect("Suffix array construction failed!"))
                as RawSuffixArray;
        let bwt = bwt(&seq, &sa);
        // use a subset of the required alphabet for FMD index
        let alpha = Alphabet::new(b"ACGNT");
        let less = less(&bwt, &alpha);
        let occ = Occ::new(&bwt, occ_sampling_rate as u32);
        let sa = sa.sample(&seq, bwt, less, occ, sa_sampling_rate);

        Ok(Index { refs, sa, txome })
    }

    fn create_txome(ref_path: &str, annot_path: &str, name_to_ref: &HashMap<NameStrand, usize>, refs: &[Ref]) -> Result<Txome> {
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

        let mut gene_intervals = vec![(usize::MAX, 0); genes.len()];
        let err_msg = format!("Error in reading reference file {}", ref_path);
        let txs = txome_txs
            .into_iter()
            .map(|tx| {
                let gene_idx = tx.gene_idx.0 as usize;
                let mut tx_seq = Vec::with_capacity(tx.len() as usize);
                tx.get_sequence(&mut ref_fai_reader, &mut tx_seq)
                    .expect(&err_msg);
                tx_seq.make_ascii_uppercase();

                let strand = tx.strand == ReqStrand::Forward;
                let tx_ref = &refs[name_to_ref[&NameStrand(tx.chrom.clone(), strand)]];

                let tx_start = if strand {
                    (tx.start() as usize) + tx_ref.start_idx
                } else {
                    tx_ref.end_idx - 1 - (tx.end() as usize)
                };
                let tx_end = if strand {
                    (tx.end() as usize) + tx_ref.start_idx
                } else {
                    tx_ref.end_idx - 1 - (tx.start() as usize)
                };
                gene_intervals[gene_idx] = (
                    gene_intervals[gene_idx].0.min(tx_start),
                    gene_intervals[gene_idx].1.max(tx_end),
                );

                let mut exons = tx
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

                        Exon {
                            start: exon_start as usize,
                            end: exon_end as usize,
                            tx_idx: exon_tx_idx,
                        }
                    })
                    .collect::<Vec<_>>();
                if !strand {
                    // reverse exons because tx sequence will be reversed
                    exons.reverse();
                }

                Tx {
                    id: tx.id,
                    chrom: tx.chrom,
                    strand,
                    exons,
                    seq: tx_seq,
                    gene_idx,
                }
            })
            .collect::<Vec<_>>();

        let gene_intervals = IntervalTree::from_iter(
            gene_intervals
                .into_iter()
                .enumerate()
                .map(|(i, (start, end))| (Interval::new(start..end).unwrap(), i)),
        );

        Ok(Txome {
            genes,
            txs,
            gene_intervals,
        })
    }

    /// Find all SMEMs for a query sequence.
    ///
    /// The MEMs use concatenated reference coordinates.
    pub fn all_smems(&self, query: &[u8], min_seed_len: usize) -> Vec<Mem> {
        let mut mems = Vec::new();
        // creating the fmd index on the fly here is fast since the structs are just wrappers
        let fm = FMIndex::new(self.sa.bwt(), self.sa.less(), self.sa.occ());
        // safe because we only align nucleotides
        let fmd = unsafe { FMDIndex::from_fmindex_unchecked(fm) };
        let intervals = fmd.all_smems(query, min_seed_len);

        for interval in intervals {
            let forwards_idxs = interval.0.forward().occ(&self.sa);
            let query_idx = interval.1;
            let mem_len = interval.2;

            for ref_idx in &forwards_idxs {
                mems.push(Mem {
                    query_idx,
                    ref_idx: *ref_idx,
                    len: mem_len,
                });
            }
        }

        // sorting makes later operations faster
        mems.sort_by_key(|mem| mem.len);
        // longest MEMs first
        mems.reverse();
        mems
    }

    /// Find a single longest supermaximal exact match (SMEM) for a query sequence.
    ///
    /// The MEMs use concatenated reference coordinates.
    pub fn longest_smem(&self, query: &[u8], min_seed_len: usize) -> Option<Mem> {
        let mut max_smem = None;
        // creating the fmd index on the fly here is fast since the structs are just wrappers
        let fm = FMIndex::new(self.sa.bwt(), self.sa.less(), self.sa.occ());
        // safe because we only align nucleotides
        let fmd = unsafe { FMDIndex::from_fmindex_unchecked(fm) };
        let intervals = fmd.all_smems(query, min_seed_len);

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

    /// Get the reference (chromosome) that contains a concatenated coordinate (index)
    /// and convert it to an index within the reference.
    pub fn idx_to_ref(&self, idx: usize) -> (&Ref, usize) {
        let ref_idx = self.refs.partition_point(|x| x.end_idx <= idx);
        (&self.refs[ref_idx], idx - self.refs[ref_idx].start_idx)
    }

    /// Get the references (chromosomes).
    pub fn refs(&self) -> &[Ref] {
        &self.refs
    }

    /// Get the transcriptome.
    pub fn txome(&self) -> &Txome {
        &self.txome
    }

    /// Get a subsequence of a reference chromosome based on an interval
    /// specified with concatenated coordinates.
    pub fn seq_slice(&self, start: usize, end: usize) -> Cow<[u8]> {
        let ref_idx = self.refs.partition_point(|x| x.end_idx <= start);
        let curr_ref = &self.refs[ref_idx];
        // don't worry about slicing a transcript sequence
        let seq = match &curr_ref.ref_type {
            RefType::Chr { seq } => seq,
            _ => unreachable!(),
        };

        match seq {
            Some(s) => {
                let chrom_start = start - curr_ref.start_idx;
                let chrom_end = end - curr_ref.start_idx;
                Cow::Borrowed(&s[chrom_start..chrom_end])
            }
            None => {
                // the reference sequence before a revcomp sequence must be the corresponding forwards sequence
                let prev_ref = &self.refs[ref_idx - 1];
                let chrom_start = curr_ref.end_idx - 1 - end;
                let chrom_end = curr_ref.end_idx - 1 - start;
                let prev_seq = match &prev_ref.ref_type {
                    RefType::Chr { seq } => seq,
                    _ => unreachable!(),
                };
                Cow::Owned(dna::revcomp(
                    &prev_seq.as_ref().unwrap()[chrom_start..chrom_end],
                ))
            }
        }
    }

    /// Print out stats about the index.
    pub fn print_stats(&self) {
        println!("Number of chromosomes\t{}", self.refs.len());
        println!(
            "Length of suffix array\t{}",
            self.sa.len() / self.sa.sampling_rate()
        );
        println!("Length of BWT\t{}", self.sa.bwt().len());
        println!("Length of Less\t{}", self.sa.less().len());
        println!("Number of genes\t{}", self.txome.genes.len());
        println!("Number of transcripts\t{}", self.txome.txs.len());

        println!();

        use bincode::serialized_size;
        println!(
            "Serialized size of index\t{}",
            serialized_size(&self).unwrap()
        );
        println!("\tChromosomes\t{}", serialized_size(&self.refs).unwrap());
        println!(
            "\tSuffix array + FM index\t{}",
            serialized_size(&self.sa).unwrap()
        );
        println!("\t\tBWT\t{}", serialized_size(self.sa.bwt()).unwrap());
        println!("\t\tLess\t{}", serialized_size(self.sa.less()).unwrap());
        println!("\t\tOcc\t{}", serialized_size(self.sa.occ()).unwrap());
        println!("\tTranscriptome\t{}", serialized_size(&self.txome).unwrap());
        println!(
            "\t\tGene interval tree\t{}",
            serialized_size(&self.txome.gene_intervals).unwrap()
        );
    }
}

/// Use unsafe magic to convert i64 vector to usize vector in place.
#[allow(dead_code)]
fn cast_vec_i64_to_usize(v: Vec<i64>) -> Vec<usize> {
    // double check necessary conditions before doing a pointer cast
    // size and alignment must match so the allocator can safely free the memory
    assert_eq!(mem::size_of::<i64>(), mem::size_of::<usize>());
    assert_eq!(mem::align_of::<i64>(), mem::align_of::<usize>());

    unsafe {
        let mut v = mem::ManuallyDrop::new(v);
        let ptr = v.as_mut_ptr();
        let len = v.len();
        let cap = v.capacity();

        Vec::from_raw_parts(ptr as *mut usize, len, cap)
    }
}

/// A single maximal exact match.
#[derive(Debug, Copy, Clone, PartialEq)]
pub struct Mem {
    pub ref_idx: usize,
    pub query_idx: usize,
    pub len: usize,
}

/// A reference type.
#[derive(Debug, Serialize, Deserialize)]
pub enum RefType {
    Chr { seq: Option<Vec<u8>> },
    Tx { tx_idx: usize },
}

/// A single reference (chromosome or transcript).
#[derive(Debug, Serialize, Deserialize)]
pub struct Ref {
    pub name: String,
    pub ref_type: RefType,
    pub strand: bool,
    pub len: usize,
    pub start_idx: usize,
    pub end_idx: usize,
}

#[derive(Debug, Clone, Hash, PartialEq, Eq)]
struct NameStrand(String, bool);
