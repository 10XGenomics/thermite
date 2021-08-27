use rust_htslib::bam::{record::Record, HeaderView};

use bio::alignment::sparse::HashMapFx;

use noodles::sam;

use bincode::serialized_size;

use std::sync::Arc;

use crate::aligner::*;
use crate::aln_writer;
use crate::index::*;

/// Wrapper around core thermite aligner code to mimic Orbit (STAR).
#[derive(Clone)]
pub struct ThermiteAligner {
    index: Arc<Index>,
    pub align_opts: AlignOpts,
    header_view: HeaderView,
}

impl ThermiteAligner {
    /// Create a new thermite aligner instance.
    pub fn new(ref_path: &str, annot_path: &str) -> Self {
        let index = Arc::new(Index::create_from_files(ref_path, annot_path, 16, 128).unwrap());
        let align_opts = AlignOpts {
            min_seed_len: 20,
            min_aln_score_percent: 0.66,
            min_aln_score: 30,
            min_total_hit_len: 40,
            multimap_score_range: 1,
        };
        let header_view = {
            let mut writer = sam::Writer::new(Vec::with_capacity(64));
            writer
                .write_header(&aln_writer::build_sam_header(&index).unwrap())
                .unwrap();
            HeaderView::from_bytes(writer.get_ref())
        };
        Self {
            index,
            align_opts,
            header_view,
        }
    }

    /// Align a single read and return a rust_htslib Record
    pub fn align_read(&self, name: &[u8], read: &[u8], qual: &[u8]) -> Vec<Record> {
        let mut tx_kmer_cache = HashMapFx::default();
        let alns = align_read(&self.index, &mut tx_kmer_cache, read, &self.align_opts);

        if alns.is_empty() {
            return vec![sam_noodles_to_htslib(
                &aln_writer::unmapped_sam_record(name, read, qual).unwrap(),
                &self.header_view,
            )];
        }

        alns.iter()
            .map(|aln| {
                sam_noodles_to_htslib(
                    &aln_writer::aln_to_sam_record(&self.index, name, read, qual, aln, alns.len())
                        .unwrap(),
                    &self.header_view,
                )
            })
            .collect::<Vec<_>>()
    }

    /// Estimate the amount of memory the index takes up.
    pub fn est_mem(&self) -> usize {
        serialized_size(self.index.as_ref()).unwrap() as usize
    }
}

/// Convert a Noodles SAM record to a rust_htslib Record.
fn sam_noodles_to_htslib(noodles_sam: &sam::Record, header_view: &HeaderView) -> Record {
    // TODO: directly create rust_htslib Record
    let mut writer = sam::Writer::new(Vec::with_capacity(64));
    writer.write_record(noodles_sam).unwrap();
    Record::from_sam(header_view, writer.get_ref()).unwrap()
}
