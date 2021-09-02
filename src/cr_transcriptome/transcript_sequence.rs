use crate::transcriptome::Exon;
use bio::io::fasta::IndexedReader;
use bio_types::strand::ReqStrand;
use failure::Error;

/// Construct a transcript sequence from a FASTA file, a strand, and a set of intervals
/// The intervals must be ordered according to the order of the intervals on the genome
pub(crate) fn get_transcript_sequence<R>(
    reader: &mut IndexedReader<R>,
    contig: &str,
    exons: &[Exon],
    strand: ReqStrand,
    result: &mut Vec<u8>,
) -> Result<(), Error>
where
    R: std::io::Read + std::io::Seek,
{
    let tx_len = exons.iter().map(|i| i.len() as usize).sum();

    result.clear();
    result.reserve(tx_len);

    let mut buf = Vec::new();

    for e in exons {
        reader.fetch(contig, e.start.into(), e.end.into())?;

        buf.clear();
        reader.read(&mut buf)?;
        result.extend(&buf);
    }

    if strand == ReqStrand::Reverse {
        revcomp_slice(&mut result[..])
    }

    Ok(())
}

fn revcomp_slice(dna: &mut [u8]) {
    for i in 0..(dna.len() + 1) / 2 {
        let j = dna.len() - 1 - i;
        let fst = bio::alphabets::dna::complement(dna[i]);
        let lst = bio::alphabets::dna::complement(dna[j]);

        dna[i] = lst;
        dna[j] = fst;
    }
}

#[cfg(test)]
mod test {
    use super::*;
    use crate::Transcriptome;
    use std::collections::HashMap;

    fn check_rc(fwd: &[u8], rev: &[u8]) {
        let mut fwd_copy = Vec::from(fwd);
        revcomp_slice(&mut fwd_copy[..]);
        assert_eq!(&fwd_copy[..], rev);
    }

    #[test]
    fn test_rc_slice() {
        check_rc(b"ACGT", b"ACGT");
        check_rc(b"A", b"T");
        check_rc(b"GC", b"GC");
        check_rc(b"GCA", b"TGC");
        check_rc(b"TTT", b"AAA");
        check_rc(b"GTGTT", b"AACAC");
    }

    /// Run this test manually to check that transcript sequences produced by Rust from GTF + genome fasta
    /// match what gencode produces.
    #[ignore]
    #[test]
    fn loader() -> Result<(), Error> {
        //let rdr = utils::open_gzip("test/gtf/GRCh38_small.gtf.gz");
        let rdr =
            std::fs::File::open("/Users/patrick/refdata_cellranger/GRCh38-3.0.0/genes/genes.gtf")?;
        let txome = Transcriptome::from_reader(std::io::BufReader::new(rdr))?;

        let mut fa = bio::io::fasta::IndexedReader::from_file(
            &"/Users/patrick/refdata_cellranger/GRCh38-3.0.0/fasta/genome.fa",
        )
        .unwrap();

        let txs = "/Users/patrick/code/rust-pseudoaligner/test/gencode.v28.transcripts.fa";
        let gencode_txs = bio::io::fasta::Reader::from_file(txs).unwrap();

        let mut tx_seqs = HashMap::new();

        for r in gencode_txs.records() {
            let r = r?;

            let id = r.id();
            let tx_id = id.split('.').next().unwrap();
            tx_seqs.insert(tx_id.to_string(), r.seq().to_vec());
        }

        let mut loaded_tx_seq = Vec::new();
        let mut not_found = 0;

        for tx in &txome.transcripts {
            let tx_id = &tx.id;

            get_transcript_sequence(&mut fa, &tx.chrom, &tx.exons, tx.strand, &mut loaded_tx_seq)?;

            if let Some(real_seq) = tx_seqs.get(tx_id) {
                //println!("found tx_id in gencode: {}", tx_id);
                // compare gencode to loaded
                if &loaded_tx_seq != real_seq {
                    println!("tx_id: {}", tx_id);
                    println!("gencode: {}", std::str::from_utf8(&real_seq).unwrap());
                    println!("me: {}", std::str::from_utf8(&loaded_tx_seq).unwrap());
                    assert!(false);
                }
            } else {
                //println!("couldn't find tx_id in gencode: {}", tx_id);
                not_found += 1;
            }
        }

        println!(
            "tried to validate {}, not found: {}",
            &txome.transcripts.len(),
            not_found
        );

        Ok(())
    }
}
