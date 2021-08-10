# thermite
RNA aligner in Rust.

## Running
```
cd data

samtools faidx test_ref.fasta

cargo run -- index -o test_ref.fasta.tai test_ref.fasta test_ref.gtf

cargo run -- align -s3 -k3 --min-total-hit-len=3 -o test_query.paf test_ref.fasta.tai test_query.fastq
# or
cargo run -- align -a -s3 -k3 --min-total-hit-len=3 -o test_query.sam test_ref.fasta.tai test_query.fastq
```
Or
```
make -C data
```
