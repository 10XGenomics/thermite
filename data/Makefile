# Test the Thermite aligner.

# Path to the thermite executable.
export PATH:=../target/debug:${PATH}

# The reference genome.
ref=test_ref

# The minimum seed length.
k=3

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

# Index the reference genome.
%.fasta.tai: %.fasta
	thermite index -o $@ $<

# Align the query to the reference (output paf).
%.paf: $(ref).fasta.tai %.fastq
	thermite align -k$k -o $@ $^

# Align the query to the reference (output sam).
%.sam: $(ref).fasta.tai %.fastq
	thermite align -a -k$k -o $@ $^

# Align the query to the reference (output bam).
%.bam: $(ref).fasta.tai %.fastq
	thermite align -a -k$k -o $@ $^
