use anyhow::Result;

use std::fmt;
use std::io::Write;

pub struct PafEntry<'a> {
    pub query_name: &'a [u8],
    pub query_len: usize,
    pub query_start: usize,
    pub query_end: usize,
    pub strand: bool,
    pub target_name: &'a [u8],
    pub target_len: usize,
    pub target_start: usize,
    pub target_end: usize,
    pub num_match: usize,
    pub num_match_gap: usize,
    pub map_qual: u8,
}

impl<'a> fmt::Display for PafEntry<'a> {
    fn fmt(&self, f: &mut fmt::Formatter) -> fmt::Result {
        write!(
            f,
            "{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t",
            String::from_utf8_lossy(self.query_name),
            self.query_len,
            self.query_start,
            self.query_end,
            if self.strand { "+" } else { "-" },
            String::from_utf8_lossy(self.target_name),
            self.target_len,
            self.target_start,
            self.target_end,
            self.num_match,
            self.num_match_gap,
            self.map_qual,
        )
    }
}

pub fn write_paf(writer: &mut dyn Write, paf_entry: &PafEntry) -> Result<()> {
    writeln!(writer, "{}", paf_entry)?;
    Ok(())
}
