# thermite
Spliced short read aligner implemented in Rust.

## Running with CLI
First:
```
cargo build --release
```

Generating PAF, SAM, and BAM alignments with some basic test reference/reads:
```
make -C data
```

Computing metrics:
```
make -C data metrics dataset=chrM
make -C data metrics dataset=chr21
cat data/chr21_comparison_metrics.txt data/chrM_comparison_metrics.txt > data/comparison_metrics.txt
```

## Running with Cellranger
Use the thermite [branch](https://github.com/10XDev/cellranger/tree/dl/thermite) of Cellranger.

First, use thermite to build an index on your reference. Cellranger expects a `reference_path`
directory, and the thermite reference should be located at `thermite/index.tai` in the directory.

Next, build Cellranger using the thermite branch and generate the `.mro` file for a dataset.
Edit the `.mro` file to make sure that the thermite index is in the `reference_path` and
set `aligner = "thermite"`.

Finally, run the entire Cellranger pipestance with the modified `.mro` file.

## Algorithm
Thermite has two main modes: indexing and aligning.

Thermite builds an index on a reference genome and transcriptome, so it needs both the
genome sequence and the transcriptome annotations. During indexing, a new sequence is built
by concatenating chromosomes (forwards and revcomp) and transcripts (forwards and revcomp).
Then, an FMD-index is built on this sequence.

When aligning a read, SuperMaximal Exact Matches between the read and the concatenated
sequence are found in thermite by using the FMD-index. Each of the SMEM seeds are extended
through SIMD DP to get full alignments. Alignments that hit transcript sequences are lifted
over to their locations in the genome, producing spliced genome alignments. Finally, overlapping
alignments from a read are deduplicated, with transcriptomic alignments being favored over
intronic or intergenic alignments.
