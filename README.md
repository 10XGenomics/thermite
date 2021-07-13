# thermite
RNA aligner in Rust.

## Running
```
cd data
cargo run index test_ref.fasta test_index.tai
cargo run align test_index.tai test_query.fastq test_aln.paf
```
