use failure::{Error, ResultExt};

use crate::transcriptome::Transcriptome;
use serde::{de::DeserializeOwned, Deserialize, Serialize};
use std::collections::HashMap;
use std::path::Path;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TranscriptIndex {
    pub transcript_genes: HashMap<String, Gene>,
    pub transcript_lengths: HashMap<String, i64>,
}

#[derive(Hash, Eq, PartialEq, Debug, Clone, Ord, PartialOrd, Serialize, Deserialize)]
pub struct Gene {
    pub id: String,
    pub name: String,
}

pub fn read_delimited_table<T: DeserializeOwned>(
    file: impl AsRef<Path>,
    delimiter: u8,
    has_headers: bool,
    skip_row: bool,
) -> Result<Vec<T>, Error> {
    let mut rdr = csv::ReaderBuilder::new()
        .delimiter(delimiter)
        .has_headers(has_headers)
        .flexible(true)
        .from_path(file.as_ref())
        .context(format!("Error: {}", file.as_ref().display()))?;
    let mut iter = rdr.deserialize::<T>();
    if skip_row {
        let _ = iter.next();
    }
    let rows = iter.collect::<csv::Result<Vec<T>>>()?;
    Ok(rows)
}

pub fn read_tabular<T: DeserializeOwned>(
    file: impl AsRef<Path>,
    skip_row: bool,
) -> Result<Vec<T>, Error> {
    read_delimited_table(file, b'\t', false, skip_row)
}

// A summary representation of the Gene/Transcript table. Used in previous versions of CR.
// retained here for testing purposes
impl TranscriptIndex {
    pub fn new(tab_file: impl AsRef<Path>) -> Result<TranscriptIndex, Error> {
        let mut transcript_lengths = HashMap::new();
        let mut transcript_genes = HashMap::new();
        let tx_index: Vec<GeneTranscriptRow> = read_tabular(tab_file, true)?;
        for tx in tx_index {
            transcript_lengths.insert(tx.transcript_id.clone(), tx.transcript_len);
            transcript_genes.insert(
                tx.transcript_id.clone(),
                Gene {
                    id: tx.gene_id,
                    name: tx.gene_name,
                },
            );
        }

        Ok(TranscriptIndex {
            transcript_genes,
            transcript_lengths,
        })
    }

    pub fn get_gene_from_transcript(&self, tx: &str) -> &Gene {
        self.transcript_genes.get(tx).expect(&format!(
            "Error: No corresponding gene for transcript: {}",
            tx
        ))
    }

    pub fn get_transcript_length(&self, tx: &str) -> &i64 {
        self.transcript_lengths.get(tx).unwrap()
    }

    pub fn from_transcriptome(txome: &Transcriptome) -> TranscriptIndex {
        let mut transcript_lengths = HashMap::new();
        let mut transcript_genes = HashMap::new();

        for tx in &txome.transcripts {
            transcript_lengths.insert(tx.id.clone(), tx.len() as i64);

            let g = &txome.genes[tx.gene_idx.0 as usize];
            let g = Gene {
                id: g.id.clone(),
                name: g.name.clone(),
            };

            transcript_genes.insert(tx.id.clone(), g);
        }

        TranscriptIndex {
            transcript_genes,
            transcript_lengths,
        }
    }
}

#[derive(Deserialize)]
struct GeneTranscriptRow {
    transcript_id: String,
    gene_id: String,
    gene_name: String,
    transcript_len: i64,
}
