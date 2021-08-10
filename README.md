# thermite
RNA aligner in Rust.

## Running
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
