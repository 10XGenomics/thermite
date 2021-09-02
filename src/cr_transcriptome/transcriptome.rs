//
// Copyright (c) 2017 10x Genomics, Inc. All rights reserved.
//

use bio::io::fasta::IndexedReader;
use bio_types::strand::ReqStrand;
use failure::{format_err, Error};

use std::collections::{hash_map::Entry, BTreeMap, HashMap, HashSet};
use std::io::{Read, Seek};
use std::path::{Path, PathBuf};

use crate::Gene;

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Copy)]
pub struct GeneIdx(pub u32);

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Copy)]
pub struct TranscriptIdx(pub u32);

#[derive(Eq, PartialEq, Debug, Clone, Hash)]
pub struct TranscriptomeGene {
    pub idx: GeneIdx,
    pub id: String,
    pub name: String,
    pub properties: Vec<(String, String)>,
}

impl TranscriptomeGene {
    pub fn to_gene(&self) -> Gene {
        Gene {
            id: self.id.clone(),
            name: self.name.clone(),
        }
    }
}

#[derive(Debug)]
pub struct Transcript {
    pub idx: TranscriptIdx,
    pub gene_idx: GeneIdx,
    pub id: String,
    pub chrom: String,
    pub strand: ReqStrand,
    pub exons: Vec<Exon>,
    pub properties: Vec<(String, String)>,
}

impl Transcript {
    pub fn get_sequence<R: Read + Seek>(
        &self,
        fasta_reader: &mut IndexedReader<R>,
        result: &mut Vec<u8>,
    ) -> Result<(), Error> {
        let res = crate::transcript_sequence::get_transcript_sequence(
            fasta_reader,
            &self.chrom,
            &self.exons,
            self.strand,
            result,
        );

        assert_eq!(
            result.len() as u64,
            self.len(),
            "Expected a transcript length of {}, but got {} stitching the exons for {:#?}",
            self.len(),
            result.len(),
            self
        );
        res
    }

    pub fn len(&self) -> u64 {
        self.exons.iter().map(|e| e.len()).sum()
    }

    pub fn gc_content<R: Read + Seek>(
        &self,
        fasta_reader: &mut IndexedReader<R>,
    ) -> Result<f64, Error> {
        let tx_len = self.len();
        let mut buf = Vec::with_capacity(tx_len as usize);
        self.get_sequence(fasta_reader, &mut buf)?;

        let mut gc = 0;
        for c in &buf {
            if *c == b'c' || *c == b'C' || *c == b'g' || *c == b'G' {
                gc += 1;
            }
        }

        if tx_len > 0 {
            Ok((gc as f64) / (tx_len as f64))
        } else {
            Ok(0.0)
        }
    }

    pub fn is_empty(&self) -> bool {
        !self.exons.iter().any(|e| !e.is_empty())
    }

    pub fn start(&self) -> u64 {
        if self.exons.len() == 0 {
            return 0;
        }
        self.exons[0].start
    }

    pub fn end(&self) -> u64 {
        if self.exons.len() == 0 {
            return 0;
        }
        self.exons[self.exons.len() - 1].end
    }
}

pub struct Transcriptome {
    pub genes: Vec<TranscriptomeGene>,
    pub gene_id_to_idx: HashMap<String, GeneIdx>,
    pub transcripts: Vec<Transcript>,
    pub transcript_id_to_idx: HashMap<String, TranscriptIdx>,
    pub gene_to_transcripts: BTreeMap<GeneIdx, Vec<TranscriptIdx>>, // map from gene_id to list of transcript_ids
}

impl Transcriptome {
    pub fn from_reference_path(path: impl AsRef<Path>) -> Result<Transcriptome, Error> {
        let mut path_buf: PathBuf = path.as_ref().into();
        path_buf.push("genes");
        path_buf.push("genes.gtf");

        let rdr = crate::utils::open_maybe_gz(path_buf)?;
        load_from_gtf_reader(std::io::BufReader::new(rdr))
    }

    pub fn from_reader<R: Read + std::io::BufRead>(reader: R) -> Result<Transcriptome, Error> {
        load_from_gtf_reader(reader)
    }

    pub fn dummy() -> Transcriptome {
        Transcriptome {
            genes: Vec::new(),
            gene_id_to_idx: HashMap::new(),
            transcripts: Vec::new(),
            transcript_id_to_idx: HashMap::new(),
            gene_to_transcripts: BTreeMap::new(),
        }
    }
}

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd)]
pub struct Exon {
    pub start: u64,
    pub end: u64,
}

impl Exon {
    pub fn len(&self) -> u64 {
        self.end - self.start
    }
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }
}

/// Read a `Transcriptome` from a GTF file. Handles `gene`, `transcript` and `exon` GTF entries.
fn load_from_gtf_reader<R>(in_gtf: R) -> Result<Transcriptome, Error>
where
    R: Read + std::io::BufRead,
{
    let mut transcripts = Vec::new();
    let mut transcript_id_to_idx = HashMap::new();

    let mut genes = Vec::new();
    let mut gene_id_to_idx: HashMap<String, GeneIdx> = HashMap::new();
    // To keep track of genes that are added by parsing exons/transcripts but
    // later updated by parsing a gene
    let mut genes_not_from_file = HashSet::new();

    let mut gene_idx = 0;
    let mut transcript_idx = 0;

    for (line_num, line) in in_gtf.lines().enumerate() {
        let line = line?;
        if line.starts_with('#') {
            continue;
        }
        let make_err = |msg| {
            format_err!(
                "Error parsing GTF on line {}\nLine = '{}'\n{}",
                line_num,
                line,
                msg
            )
        };

        let rec = match crate::parse_gtf::parse_gtf_line(&line.as_bytes()) {
            Ok((_, rec)) => rec,
            Err(e) => {
                let msg = make_err(format!("{:?}", e));
                return Err(msg);
            }
        };

        // start is 1-based, convert to 0-based
        let start = rec.start - 1;
        // end is 1-based, inclusive. same value as 0-based exclusive.
        let end = rec.end;

        let strand = match rec.strand {
            b"+" => ReqStrand::Forward,
            b"-" => ReqStrand::Reverse,
            _ => return Err(make_err("unknown strand".to_string())),
        };

        let chrom = std::str::from_utf8(rec.seqname)?;

        if rec.feature_type == b"gene" {
            let id = rec.get_attr("gene_id")?;
            let name = rec.get_attr("gene_name").unwrap_or_else(|_| id.clone());
            let idx = GeneIdx(gene_idx);

            let gene = TranscriptomeGene {
                idx: idx.clone(),
                id: id.clone(),
                name,
                properties: rec.all_attributes(),
            };
            // Make sure we don't let duplicates through
            match gene_id_to_idx.entry(id) {
                Entry::Occupied(entry) => {
                    if !genes_not_from_file.remove(entry.key()) {
                        return Err(format_err!(
                            "Duplicate Gene ID found in GTF: {}",
                            entry.key()
                        ));
                    }
                    // If a `transcript` appears before a `gene`, we generate an entry for that gene
                    // but we should override that with the actual gene once we observe it.
                    let idx_int = entry.get().0 as usize;
                    genes[idx_int] = gene;
                }
                Entry::Vacant(entry) => {
                    entry.insert(idx);
                    genes.push(gene);
                    gene_idx += 1;
                }
            }
        } else if rec.feature_type == b"transcript" {
            let id = rec.get_attr("transcript_id")?;
            let gene_id = rec.get_attr("gene_id")?;
            let idx = TranscriptIdx(transcript_idx);

            // handle a missing gene for this tx
            let gene_idx = gene_id_to_idx.entry(gene_id.clone()).or_insert_with(|| {
                let gene_name = rec
                    .get_attr("gene_name")
                    .unwrap_or_else(|_| gene_id.clone());
                let new_gene_idx = GeneIdx(gene_idx);

                let gene = TranscriptomeGene {
                    idx: new_gene_idx,
                    id: gene_id.clone(),
                    name: gene_name,
                    properties: vec![],
                };
                genes_not_from_file.insert(gene_id.clone());
                genes.push(gene);
                gene_idx += 1;
                new_gene_idx
            });

            let transcript = Transcript {
                idx: idx.clone(),
                id: id.clone(),
                gene_idx: gene_idx.clone(),
                chrom: chrom.to_string(),
                strand,
                exons: vec![],
                properties: rec.all_attributes(),
            };

            transcript_id_to_idx.insert(id, idx);
            transcripts.push(transcript);
            transcript_idx += 1;
            continue;
        } else if rec.feature_type == b"exon" {
            let exon = Exon { start, end };
            let transcript_id = rec.get_attr("transcript_id")?;

            // the transcript hasn't been declared -- make a  dummy
            if transcript_id_to_idx.get(&transcript_id).is_none() {
                let gene_id = rec.get_attr("gene_id")?;
                let idx = TranscriptIdx(transcript_idx);

                // handle a missing gene for this tx
                if gene_id_to_idx.get(&gene_id).is_none() {
                    let gene_name = rec
                        .get_attr("gene_name")
                        .unwrap_or_else(|_| gene_id.clone());
                    let new_gene_idx = GeneIdx(gene_idx);

                    let gene = TranscriptomeGene {
                        idx: new_gene_idx,
                        id: gene_id.clone(),
                        name: gene_name,
                        properties: vec![],
                    };
                    genes_not_from_file.insert(gene_id.clone());
                    gene_id_to_idx.insert(gene_id.clone(), new_gene_idx);
                    genes.push(gene);
                    gene_idx += 1;
                }

                let gene_idx = gene_id_to_idx.get(&gene_id).expect("missing gene");

                let transcript = Transcript {
                    idx: idx.clone(),
                    id: transcript_id.clone(),
                    gene_idx: gene_idx.clone(),
                    chrom: chrom.to_string(),
                    strand,
                    exons: vec![],
                    properties: vec![],
                };

                transcript_id_to_idx.insert(transcript_id.clone(), idx);
                transcripts.push(transcript);
                transcript_idx += 1;
            }

            let tx_idx = transcript_id_to_idx.get(&transcript_id).ok_or_else(|| {
                make_err(format!(
                    "this row references transcript_id={} but this \
                    transcript has no preceding 'transcript' row in the GTF",
                    transcript_id
                ))
            })?;

            transcripts[tx_idx.0 as usize].exons.push(exon);
        }
    }

    let mut gene_to_transcripts = BTreeMap::new();
    for tx in &mut transcripts {
        // Sort the exons so they're in coordinate order on the genome
        tx.exons.sort();

        // Tabulate the set of transcripts for each gene
        gene_to_transcripts
            .entry(tx.gene_idx)
            .or_insert_with(|| vec![])
            .push(tx.idx);
    }

    Ok(Transcriptome {
        genes,
        gene_id_to_idx,
        transcripts,
        transcript_id_to_idx,
        gene_to_transcripts,
    })
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::transcript_index::TranscriptIndex;
    use failure::Error;

    // manual test that the GTF path matches the gene_index_tab created in the
    // legacy python code.
    #[test]
    fn test_tx_index() -> Result<(), Error> {
        let f = std::fs::File::open("test/gtf/GRCh38_small.gtf.gz")
            .expect("Failed to open file for reading");
        let rdr = flate2::read::GzDecoder::new(f);
        let txome = Transcriptome::from_reader(std::io::BufReader::new(rdr))?;

        let gene_index_table = "test/gtf/GRCh38_small_gene_index.tab";

        let old = TranscriptIndex::new(gene_index_table)?;
        let new = TranscriptIndex::from_transcriptome(&txome);

        assert_eq!(old, new);
        Ok(())
    }

    #[test]
    fn test_no_gene_name() -> Result<(), Error> {
        //let f = std::fs::File::open("/mnt/opt/refdata_cellranger/GRCh38_and_Sscrofa11.1/genes/genes.gtf")
        let f = std::fs::File::open("test/gtf/no_gene_names.gtf")
            .expect("Failed to open file for reading");
        let _txome = Transcriptome::from_reader(std::io::BufReader::new(f))?;
        Ok(())
    }

    // test duplicates can't be parsed.
    #[test]
    fn test_no_dups() -> Result<(), Error> {
        let gtf = "\
NC_000086.8	BestRefSeq	gene	168758039	168761913	.	-	.	gene_id \"G530011O06Rik\"
NC_000087.8	BestRefSeq	gene	90762409	90766319	.	-	.	gene_id \"G530011O06Rik\"; ";
        let txome = load_from_gtf_reader(std::io::BufReader::new(gtf.as_bytes()));
        let expected = "Duplicate Gene ID found in GTF: G530011O06Rik";
        assert_eq!(txome.is_err(), true);
        if let Err(e) = txome {
            assert_eq!(format!("{}", e), expected)
        }
        // This should succeed
        let alt_gtf = "\
NC_000087.8	BestRefSeq	exon	90765690	90766319	.	-	.	gene_id \"G530011O06Rik\"; transcript_id \"NR_137283.1\";
NC_000086.8	BestRefSeq	gene	168758039	168761913	.	-	.	gene_id \"G530011O06Rik\";";
        let txome2 = load_from_gtf_reader(std::io::BufReader::new(alt_gtf.as_bytes()));
        assert_eq!(txome2.is_ok(), true);
        Ok(())
    }
}
