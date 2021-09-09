use bio::alignment::pairwise::{MatchFunc, Scoring, MIN_SCORE};
use bio::alignment::{Alignment, AlignmentMode, AlignmentOperation};

use block_aligner::cigar::*;
use block_aligner::scan_block::*;
use block_aligner::scores::*;

/// Smith-Waterman-Gotoh banded extension alignment.
#[allow(non_snake_case)]
pub struct SwgExtend<F: MatchFunc> {
    D: Vec<i32>,
    C: Vec<i32>,
    R: Vec<i32>,
    trace: Vec<AlignmentOperation>,
    scoring: Scoring<F>,
    max_band_width: usize,
}

impl<F: MatchFunc> SwgExtend<F> {
    /// Allocate space for alignments up to a certain band size.
    pub fn new(max_band_width: usize, scoring: Scoring<F>) -> Self {
        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                return Self {
                    D: Vec::new(),
                    C: Vec::new(),
                    R: Vec::new(),
                    trace: Vec::new(),
                    scoring,
                    max_band_width,
                };
            }
        }

        Self {
            D: vec![0; max_band_width * 2 + 1],
            C: vec![0; max_band_width * 2 + 1],
            R: vec![0; max_band_width * 2 + 1],
            trace: Vec::new(),
            scoring,
            max_band_width,
        }
    }

    /// Extension alignment using SIMD-accelerated block aligner.
    fn extend_simd(&mut self, x: &[u8], y: &[u8], band_width: usize, x_drop: i32) -> Alignment {
        let block_size = if band_width <= 16 { 16 } else { 32 };
        let gaps = Gaps {
            open: (self.scoring.gap_open + self.scoring.gap_extend) as i8,
            extend: self.scoring.gap_extend as i8,
        };
        let q = PaddedBytes::from_bytes::<NucMatrix>(x, block_size);
        let r = PaddedBytes::from_bytes::<NucMatrix>(y, block_size);
        let matrix = NucMatrix::new_simple(
            self.scoring.match_fn.score(b'A', b'A') as i8,
            self.scoring.match_fn.score(b'A', b'C') as i8,
        );

        // align with X-drop threshold and store trace
        let a =
            Block::<_, true, true>::align(&q, &r, &matrix, gaps, block_size..=block_size, x_drop);
        let res = a.res();
        let cigar = a.trace().cigar(res.query_idx, res.reference_idx);
        let operations = {
            let mut o = Vec::with_capacity(x.len() + y.len() + 4);
            let mut i = 0;
            let mut j = 0;

            for k in 0..cigar.len() {
                let op_len = cigar.get(k);

                match op_len.op {
                    Operation::M => {
                        for _ in 0..op_len.len {
                            o.push(if x[i] == y[j] {
                                AlignmentOperation::Match
                            } else {
                                AlignmentOperation::Subst
                            });
                            i += 1;
                            j += 1;
                        }
                    }
                    Operation::I => {
                        o.resize(o.len() + op_len.len, AlignmentOperation::Ins);
                        i += op_len.len;
                    }
                    Operation::D => {
                        o.resize(o.len() + op_len.len, AlignmentOperation::Del);
                        j += op_len.len;
                    }
                    _ => unreachable!(),
                };
            }

            if res.query_idx < x.len() {
                o.push(AlignmentOperation::Xclip(x.len() - res.query_idx));
            }

            o
        };

        Alignment {
            score: res.score,
            ystart: 0,
            xstart: 0,
            yend: res.reference_idx,
            xend: res.query_idx,
            ylen: y.len(),
            xlen: x.len(),
            operations,
            mode: AlignmentMode::Custom,
        }
    }

    /// Align with a certain band width and X-drop threshold, ending the alignment
    /// at the highest scoring index.
    #[allow(non_snake_case)]
    pub fn extend(&mut self, x: &[u8], y: &[u8], band_width: usize, x_drop: i32) -> Alignment {
        assert!(
            band_width <= self.max_band_width,
            "Band width of {} must be less than the max band width of {}!",
            band_width,
            self.max_band_width
        );

        if x.len() == 0 || y.len() == 0 {
            return Alignment {
                score: 0,
                ystart: 0,
                xstart: 0,
                yend: 0,
                xend: 0,
                ylen: y.len(),
                xlen: x.len(),
                operations: if x.len() > 0 {
                    vec![AlignmentOperation::Xclip(x.len())]
                } else {
                    Vec::new()
                },
                mode: AlignmentMode::Custom,
            };
        }

        #[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
        {
            if is_x86_feature_detected!("avx2") {
                return self.extend_simd(x, y, band_width, x_drop);
            }
        }

        let w = band_width * 2 + 1;
        let mut max_score = 0i32;
        let mut max_idx = (0usize, 0usize);

        // initialize leftmost column
        self.D[0] = 0;
        self.C[0] = 0;
        self.R[0] = 0;
        self.set_trace(0, 0, AlignmentOperation::Ins);
        for i in 1..w {
            self.C[i] = MIN_SCORE;
            self.R[i] = (i as i32) * self.scoring.gap_extend + self.scoring.gap_open;
            self.D[i] = self.R[i];
            self.set_trace(0, i, AlignmentOperation::Ins);
        }

        // handle first couple of columns where the band always starts at the
        // first row of the full DP matrix
        for j in 1..=band_width.min(y.len()) {
            let mut band_max = MIN_SCORE;
            let mut prev_D = MIN_SCORE;

            // compute each cell in the column
            for i in 0..w.min(x.len() + 1) {
                self.C[i] = (self.C[i] + self.scoring.gap_extend)
                    .max(self.D[i] + self.scoring.gap_extend + self.scoring.gap_open);
                self.R[i] = if i == 0 {
                    MIN_SCORE
                } else {
                    (self.R[i - 1] + self.scoring.gap_extend)
                        .max(self.D[i - 1] + self.scoring.gap_extend + self.scoring.gap_open)
                };
                let d = if i == 0 {
                    MIN_SCORE
                } else {
                    prev_D + self.scoring.match_fn.score(x[i - 1], y[j - 1])
                };
                prev_D = self.D[i];

                let (curr_D, dir) =
                    self.triple_max(d, self.C[i], self.R[i], i > 0 && x[i - 1] == y[j - 1]);
                self.D[i] = curr_D;
                self.set_trace(j, i, dir);

                if self.D[i] > max_score {
                    max_score = self.D[i];
                    max_idx = (i, j);
                }

                band_max = band_max.max(self.D[i]);
            }

            // x drop
            if band_max < max_score - x_drop {
                break;
            }
        }

        // handle the rest of the columns that shift down in each iteration
        for j in (band_width + 1)..y.len() + 1 {
            let mut band_max = MIN_SCORE;

            for i in (j - band_width)..(j - band_width + w).min(x.len() + 1) {
                // compute the index of the current cell within the column
                let band_idx = i - (j - band_width);

                self.C[band_idx] = if band_idx >= w - 1 {
                    MIN_SCORE
                } else {
                    (self.C[band_idx + 1] + self.scoring.gap_extend)
                        .max(self.D[band_idx + 1] + self.scoring.gap_extend + self.scoring.gap_open)
                };
                self.R[band_idx] = if band_idx == 0 {
                    MIN_SCORE
                } else {
                    (self.R[band_idx - 1] + self.scoring.gap_extend)
                        .max(self.D[band_idx - 1] + self.scoring.gap_extend + self.scoring.gap_open)
                };
                let d = self.D[band_idx] + self.scoring.match_fn.score(x[i - 1], y[j - 1]);

                let (curr_D, dir) =
                    self.triple_max(d, self.C[band_idx], self.R[band_idx], x[i - 1] == y[j - 1]);
                self.D[band_idx] = curr_D;
                self.set_trace(j, band_idx, dir);

                if self.D[band_idx] > max_score {
                    max_score = self.D[band_idx];
                    max_idx = (i, j);
                }

                band_max = band_max.max(self.D[band_idx]);
            }

            // x drop
            if band_max < max_score - x_drop {
                break;
            }
        }

        Alignment {
            score: max_score,
            ystart: 0,
            xstart: 0,
            yend: max_idx.1,
            xend: max_idx.0,
            ylen: y.len(),
            xlen: x.len(),
            operations: self.trace(max_idx.0, max_idx.1, x.len(), band_width),
            mode: AlignmentMode::Custom,
        }
    }

    /// Compute the traceback path ending at a certain position.
    fn trace(
        &self,
        mut i: usize,
        mut j: usize,
        len: usize,
        band_width: usize,
    ) -> Vec<AlignmentOperation> {
        let mut traceback = Vec::with_capacity(i + j + 4);
        if i < len {
            traceback.push(AlignmentOperation::Xclip(len - i));
        }

        while i > 0 || j > 0 {
            let band_idx = i - j.saturating_sub(band_width);
            traceback.push(self.get_trace(j, band_idx));

            match self.get_trace(j, band_idx) {
                AlignmentOperation::Match => {
                    i -= 1;
                    j -= 1;
                }
                AlignmentOperation::Subst => {
                    i -= 1;
                    j -= 1;
                }
                AlignmentOperation::Ins => {
                    i -= 1;
                }
                AlignmentOperation::Del => {
                    j -= 1;
                }
                _ => unreachable!(),
            }
        }

        traceback.reverse();
        traceback
    }

    /// Set a trace direction.
    fn set_trace(&mut self, j: usize, i: usize, op: AlignmentOperation) {
        let w = self.max_band_width * 2 + 1;
        if self.trace.len() <= j * w {
            self.trace
                .resize(self.trace.len() + w, AlignmentOperation::Match);
        }
        self.trace[j * w + i] = op;
    }

    /// Get a trace direction.
    fn get_trace(&self, j: usize, i: usize) -> AlignmentOperation {
        let w = self.max_band_width * 2 + 1;
        self.trace[j * w + i]
    }

    /// Determine the max of 3, recording the trace direction.
    fn triple_max(&self, d: i32, c: i32, r: i32, m: bool) -> (i32, AlignmentOperation) {
        let score = d.max(c).max(r);
        let dir = if score == d {
            if m {
                AlignmentOperation::Match
            } else {
                AlignmentOperation::Subst
            }
        } else if score == c {
            AlignmentOperation::Del
        } else {
            AlignmentOperation::Ins
        };
        (score, dir)
    }
}

#[cfg(test)]
mod test {
    use super::*;

    use bio::alignment::AlignmentOperation::*;

    #[test]
    fn test_swg_extend() {
        let scoring = Scoring::from_scores(-1, -1, 1, -1);
        let mut swg = SwgExtend::new(4, scoring);

        let x = b"AAAAAAAA";
        let y = b"AAAAAAAA";
        let correct_aln = Alignment {
            score: 8,
            ystart: 0,
            xstart: 0,
            yend: 8,
            xend: 8,
            ylen: 8,
            xlen: 8,
            operations: vec![Match, Match, Match, Match, Match, Match, Match, Match],
            mode: AlignmentMode::Custom,
        };
        let res_aln = swg.extend(x, y, 1, 1);
        assert_eq!(res_aln, correct_aln);

        let x = b"AAAAATTT";
        let y = b"AAAAAAAA";
        let correct_aln = Alignment {
            score: 5,
            ystart: 0,
            xstart: 0,
            yend: 5,
            xend: 5,
            ylen: 8,
            xlen: 8,
            operations: vec![Match, Match, Match, Match, Match, Xclip(3)],
            mode: AlignmentMode::Custom,
        };
        let res_aln = swg.extend(x, y, 1, 1);
        assert_eq!(res_aln, correct_aln);

        let x = b"AAATAAAA";
        let y = b"AAAAAAAA";
        let correct_aln = Alignment {
            score: 6,
            ystart: 0,
            xstart: 0,
            yend: 8,
            xend: 8,
            ylen: 8,
            xlen: 8,
            operations: vec![Match, Match, Match, Subst, Match, Match, Match, Match],
            mode: AlignmentMode::Custom,
        };
        let res_aln = swg.extend(x, y, 1, 1);
        assert_eq!(res_aln, correct_aln);

        let x = b"AAATTTT";
        let y = b"AAACCTTTT";
        let correct_aln = Alignment {
            score: 4,
            ystart: 0,
            xstart: 0,
            yend: 9,
            xend: 7,
            ylen: 9,
            xlen: 7,
            operations: vec![Match, Match, Match, Del, Del, Match, Match, Match, Match],
            mode: AlignmentMode::Custom,
        };
        let res_aln = swg.extend(x, y, 2, 3);
        assert_eq!(res_aln, correct_aln);
    }
}
