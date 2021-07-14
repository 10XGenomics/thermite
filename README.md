# thermite
RNA aligner in Rust.

## Running
```
cd data
cargo run -- index -o test_ref.fasta.tai test_ref.fasta
cargo run -- align -k 3 -o test_query.paf test_ref.fasta.tai test_query.fastq
```
