# thermite
RNA aligner in Rust.

## Running
First:
```
cargo build
```

Generating PAF, SAM, and BAM alignments:
```
make -C data
# or
make -C data dataset=chrM
```

Computing metrics:
```
make -C data metrics dataset=chrM
```
