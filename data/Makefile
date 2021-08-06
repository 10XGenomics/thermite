# Test the Thermite aligner.

# Path to the thermite executable.
export PATH:=../target/debug:${PATH}

# The reference genome.
ref=test_ref

# The minimum seed length.
k=3
# Minimum score.
s=3
# Minimum sum of the lengths of all seed hits.
min_total_hit_len=3

thermite_flags=-s $s -k $k --min-total-hit-len $(min_total_hit_len)

SHELL=bash -e -o pipefail
ifneq ($(shell command -v zsh),)
# Set pipefail to ensure that all commands of a pipe succeed.
SHELL=zsh -e -o pipefail
# Report run time and memory usage with zsh.
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
endif

# Test the Thermite aligner.
all: test_query.paf test_query.sam test_query.bam

# Remove test output.
clean:
	rm -f test_query.paf test_query.sam test_query.bam

.DELETE_ON_ERROR:
.SECONDARY:

# Index the fasta reference for random access.
%.fasta.fai: %.fasta
	samtools faidx $<

# Index the reference genome.
%.fasta.tai: %.fasta %.fasta.fai %.gtf
	thermite index -o $@ $*.fasta $*.gtf

# Align the query to the reference (output paf).
%.paf: $(ref).fasta.tai %.fastq
	thermite align $(thermite_flags) -o $@ $^

# Align the query to the reference (output sam).
%.sam: $(ref).fasta.tai %.fastq
	thermite align -a $(thermite_flags) -o $@ $^

# Align the query to the reference (output bam).
%.bam: $(ref).fasta.tai %.fastq
	thermite align -a $(thermite_flags) -o $@ $^
