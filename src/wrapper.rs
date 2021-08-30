use rust_htslib::bam::{record::Record, HeaderView};

use bio::alignment::sparse::HashMapFx;

use noodles::sam;

use bincode::deserialize_from;

use std::fs::{self, File};
use std::path::Path;
use std::sync::Arc;

use crate::aligner::*;
use crate::aln_writer;
use crate::index::*;

/// Wrapper around core thermite aligner code to mimic Orbit (STAR).
#[derive(Clone)]
pub struct ThermiteAligner {
    index: Arc<Index>,
    align_opts: AlignOpts,
    header_view: HeaderView,
}

unsafe impl Send for ThermiteAligner {}

impl ThermiteAligner {
    /// Create a new thermite aligner instance from an existing Thermite Aligner Index file.
    pub fn new(index_path: &Path) -> Self {
        let index = Arc::new(
            deserialize_from(
                File::open(index_path).expect(&format!("Failed to open {}", index_path.display())),
            )
            .unwrap(),
        );
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
        let alns = align_read(&self.index, read, &self.align_opts);

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
    pub fn est_mem(index_path: &Path) -> usize {
        fs::metadata(index_path)
            .expect(&format!("Failed to open {}", index_path.display()))
            .len() as usize
    }

    /// Get a reference to the alignment options.
    pub fn opts(&self) -> &AlignOpts {
        &self.align_opts
    }

    /// Get a mutable reference to the alignment options.
    pub fn opts_mut(&mut self) -> &mut AlignOpts {
        &mut self.align_opts
    }

    /// Get a reference to the rust_htslib HeaderView.
    pub fn header_view(&self) -> &HeaderView {
        &self.header_view
    }
}

/// Convert a Noodles SAM record to a rust_htslib Record.
fn sam_noodles_to_htslib(noodles_sam: &sam::Record, header_view: &HeaderView) -> Record {
    // TODO: directly create rust_htslib Record
    let mut writer = sam::Writer::new(Vec::with_capacity(64));
    writer.write_record(noodles_sam).unwrap();
    Record::from_sam(header_view, writer.get_ref()).unwrap()
}
