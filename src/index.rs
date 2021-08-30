use anyhow::Result;

use needletail::*;

use bio::alphabets::{dna, Alphabet};
use bio::data_structures::bwt::{bwt, less, Less, Occ, BWT};
use bio::data_structures::fmindex::{FMDIndex, FMIndex};
use bio::data_structures::interval_tree::IntervalTree;
use bio::data_structures::suffix_array::{RawSuffixArray, SampledSuffixArray, SuffixArray};
use bio::io::fasta::IndexedReader;
use bio::utils::Interval;

use bio_types::strand::ReqStrand;

use libdivsufsort_rs::divsufsort64;

use serde::{Deserialize, Serialize};

use transcriptome;

use std::collections::HashMap;
use std::fs::File;
use std::io::BufReader;
use std::str;

use crate::txome::*;

/// An index of a set of reference sequences, along with the transcriptome.
///
/// Multiple different reference sequences (chromosomes) are concatenated internally,
/// so it is important to convert the concatenated coordinates back to coordinates within
/// a specific chromosome when outputting alignments.
#[derive(Serialize, Deserialize)]
pub struct Index {
    refs: Vec<Ref>,
    seq: Vec<u8>,
    sa: SampledSuffixArray<BWT, Less, Occ>,
    txome: Txome,
}

impl Index {
    /// Create an index from a reference fasta file and its annotations file.
    ///
    /// The fasta file is expected to be already indexed, so a .fasta.fai file exists
    /// with the same file name.
    pub fn create_from_files(
        ref_path: &str,
        annot_path: &str,
        sa_sampling_rate: usize,
        occ_sampling_rate: usize,
    ) -> Result<Self> {
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

        seq.make_ascii_uppercase();

        let sa = divsufsort64(&seq)
            .expect("Suffix array construction failed!")
            .into_iter()
            .map(|x| x as usize)
            .collect::<Vec<_>>() as RawSuffixArray;
        let bwt = bwt(&seq, &sa);
        // use a subset of the required alphabet for FMD index
        let alpha = Alphabet::new(b"ACGNT");
        let less = less(&bwt, &alpha);
        let occ = Occ::new(&bwt, occ_sampling_rate as u32, &alpha);
        let sa = sa.sample(&seq, bwt, less, occ, sa_sampling_rate);

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

        let mut gene_intervals = vec![(seq.len(), 0); genes.len()];
        let mut exon_to_tx = IntervalTree::new();
        let txs = txome_txs
            .into_iter()
            .map(|tx| {
                let gene_idx = tx.gene_idx.0 as usize;
                let mut tx_seq = Vec::with_capacity(tx.len() as usize);
                tx.get_sequence(&mut ref_fai_reader, &mut tx_seq)
                    .expect(&format!("Error in reading reference file {}", ref_path));

                let strand = tx.strand == ReqStrand::Forward;
                let tx_ref = &refs[name_to_ref[&NameStrand(tx.chrom.clone(), strand)]];

                let tx_start = if strand {
                    tx.start() + tx_ref.start_idx
                } else {
                    tx_ref.end_idx - 1 - tx.end()
                };
                let tx_end = if strand {
                    tx.end() + tx_ref.start_idx
                } else {
                    tx_ref.end_idx - 1 - tx.start()
                };
                gene_intervals[gene_idx] = (gene_intervals[gene_idx].0.min(tx_start), gene_intervals[gene_idx].1.max(tx_end));

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

                        exon_to_tx
                            .insert(Interval::new(exon_start..exon_end).unwrap(), exon_tx_idx);

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

        let mut gene_intervals = IntervalTree::from_iter(gene_intervals
                                                     .into_iter()
                                                     .enumerate()
                                                     .map(|(i, (start, end))| (Interval::new(start..end).unwrap(), i)));

        let txome = Txome {
            genes,
            txs,
            exon_to_tx,
            gene_intervals,
        };

        Ok(Index {
            refs,
            seq,
            sa,
            txome,
        })
    }

    /// Find all SMEMs for a query sequence.
    ///
    /// The MEMs use concatenated reference coordinates.
    pub fn all_smems(
        &self,
        query: &[u8],
        min_seed_len: usize,
    ) -> Vec<Mem> {
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
                let seed = Interval::new(*ref_idx..*ref_idx + mem_len).unwrap();
                mems.push(Mem {
                    query_idx,
                    ref_idx: *ref_idx,
                    len: mem_len,
                });
            }
        }

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

    /// Get the concatenated reference sequence.
    pub fn seq(&self) -> &[u8] {
        &self.seq
    }

    /// Print out stats about the index.
    pub fn print_stats(&self) {
        println!("Number of chromosomes\t{}", self.refs.len());
        println!("Length of concatenated sequence\t{}", self.seq.len());
        println!(
            "Length of suffix array\t{}",
            self.sa.len() / self.sa.sampling_rate()
        );
        println!("Length of BWT\t{}", self.sa.bwt().len());
        println!("Length of Less\t{}", self.sa.less().len());
        println!("Number of genes\t{}", self.txome.genes.len());
        println!("Number of transcripts\t{}", self.txome.txs.len());

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
            "\t\tExon -> transcript interval tree\t{}",
            serialized_size(&self.txome.exon_to_tx).unwrap()
        );
    }
}

/// A single maximal exact match.
#[derive(Copy, Clone, PartialEq)]
pub struct Mem {
    pub ref_idx: usize,
    pub query_idx: usize,
    pub len: usize,
}

/// A single reference (chromosome).
#[derive(Clone, PartialEq, Serialize, Deserialize)]
pub struct Ref {
    pub name: String,
    pub strand: bool,
    pub len: usize,
    pub start_idx: usize,
    pub end_idx: usize,
}
