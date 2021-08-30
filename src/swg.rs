use bio::alignment::{Alignment, AlignmentOperation};
use bio::alignment::pairwise::{Scoring, MIN_SCORE};

/// Smith-Waterman-Gotoh banded extension alignment.
pub struct SwgExtend<F> {
    D: Vec<i32>,
    C: Vec<i32>,
    R: Vec<i32>,
    trace: Vec<Vec<AlignmentOperation>>,
    scoring: Scoring<F>,
    max_band_width: usize,
}

impl<F> SwgExtend<F> {
    pub fn new(max_band_width: usize, scoring: Scoring<F>) -> Self {
        Self {
            D: vec![0; max_band_width * 2 + 1],
            C: vec![0; max_band_width * 2 + 1],
            R: vec![0; max_band_width * 2 + 1],
            trace: Vec::new(),
            scoring,
            max_band_width,
        }
    }

    pub fn extend(&mut self, x: &[u8], y: &[u8], band_width: usize, x_drop: i32) -> Alignment {
        assert!(band_width <= self.max_band_width);

        let mut x_idx_offset = 0;
        let w = band_width * 2 + 1;

        self.D[0] = 0;
        self.C[0] = 0;
        self.R[0] = 0;
        self.set_trace(0, 0, AlignmentOperation::Insert, w);
        for i in 1..w {
            self.C[i] = MIN_SCORE;
            self.R[i] = i * self.scoring.gap_extend + self.scoring.gap_open;
            self.D[i] = self.R[i];
            self.set_trace(0, i, AlignmentOperation::Insert, w);
        }

        let mut max_score = 0i32;
        let mut max_idx = (0usize, 0usize);

        for j in 1..y.len() + 1 {
            let mut band_max = MIN_SCORE;
            let mut prev_D = MIN_SCORE;

            for i in 0..w {
                let x_idx = x_idx_offset + i;
                if x_idx > x.len() {
                    break;
                }

                self.C[i] = (self.C[i] + self.scoring.gap_extend)
                    .max(self.D[i] + self.scoring.gap_extend + self.scoring.gap_open);

                R[i] = if i == 0 { i32::MIN } else {
                    (R[i - 1] + self.scoring.gap_extend).max(D[i - 1] + self.scoring.gap_extend + self.scoring.gap_open)
                };

                let d = if x_idx == 0 {
                    i32::MIN
                } else {
                    let m = self.scoring.match_fn(x[x_idx - 1], y[j - 1]);
                    D[j - 1][i - if j > band_width { 1 } else { 0 }] + m
                };

                D[j][i] = d.max(C[i]).max(R[i]);

                if D[j][i] > max_score {
                    max_score = D[j][i];
                    max_idx = (x_idx, j);
                }

                band_max = band_max.max(D[j][i]);

                if x_idx > x.len() {
                    break;
                }
            }

            if j >= band_width {
                x_idx_offset += 1;
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
            operations: self.trace(max_idx.0, max_idx.1, band_width),
            mode: AlignmentMode::Custom,
        }
    }

    fn trace(&self, mut i: usize, mut j: usize, band_width: usize) -> Vec<AlignmentOperation> {
        let mut traceback = Vec::with_capacity(i + j + 4);

        while i > 0 || j > 0 {
            let band_idx = i - j.saturating_sub(band_width);
            traceback.push(self.trace[j][band_idx]);

            match self.trace[j][band_idx] {
                AlignmentOperation::Match {
                    i -= 1;
                    j -= 1;
                },
                AlignmentOperation::Mismatch {
                    i -= 1;
                    j -= 1;
                },
                AlignmentOperation::Insertion {
                    i -= 1;
                },
                AlignmentOperation::Deletion {
                    j -= 1;
                },
                _ => unreachable!(),
            }
        }

        traceback
    }

    fn set_trace(&mut self, j: usize, i: usize, op: AlignmentOperation, w: usize) {
        if self.trace.len() <= j {
            self.trace.push(vec![AlignmentOperation::Match; w]);
        }
        self.trace[j][i] = op;
    }
}
