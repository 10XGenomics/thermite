# Test the Thermite aligner.

# Path to the thermite executable.
export PATH:=../target/debug:${PATH}

dataset=test

ifeq ($(dataset),test)
	# Use integration test data.
	# Set the testing reference genome name.
	ref=test_ref
	# Set the testing reads file as the input.
	query=test_query.fastq
	# Set the name of the output files.
	aln=test_query
	# Set the file to compare to for metrics.
	compare=
	# Set flags for indexing.
	index_flags=
	# Set flags for alignment.
	aln_flags=-k3 --min-total-hit-len=3
else ifeq ($(dataset),chrM)
	# Use the human mitochondrial genome data.
	ref=GRCh38-2020-A-chrM
	chr=M
	query=pbmc10k_chrM_500k_reads.fastq.gz
	aln=thermite_pbmc10k_chrM_500k_reads
	compare=cellranger_pbmc10k_chrM_500k_reads.bam
	index_flags=
	aln_flags=-k20 --min-total-hit-len=40
else ifeq ($(dataset), chr21)
	# Use chromosome 21 dataset.
	ref=GRCh38-2020-A-chr21
	chr=21
	query=pbmc10k_chr21_500k_reads.fastq.gz
	aln=thermite_pbmc10k_chr21_500k_reads
	compare=cellranger_pbmc10k_chr21_500k_reads.bam
	index_flags=
	aln_flags=-k20 --min-total-hit-len=40
else
	$(error The dataset name $(dataset) is invalid!)
endif

SHELL=bash -e -o pipefail
ifneq ($(shell command -v zsh),)
# Set pipefail to ensure that all commands of a pipe succeed.
SHELL=zsh -e -o pipefail
# Report run time and memory usage with zsh.
export REPORTTIME=1
export TIMEFMT=time user=%U system=%S elapsed=%E cpu=%P memory=%M job=%J
endif

# Test the Thermite aligner.
all: $(aln).paf $(aln).sam $(aln).bam

# Alignment metrics.
metrics: chr$(chr)_comparison_metrics.txt

# Remove test output.
clean:
	rm -f $(ref).fasta.tai
	rm -f $(aln).paf $(aln).sam $(aln).bam
	rm -f comparison_metrics.txt

.DELETE_ON_ERROR:
.SECONDARY:

# Index the reference genome.
$(ref).fasta.tai: $(ref).fasta $(ref).fasta.fai $(ref).gtf
	thermite index $(index_flags) -o $@ $< $(ref).gtf

# Align the query to the reference (output paf).
$(aln).paf: $(ref).fasta.tai $(query)
	thermite align $(aln_flags) -o $@ $^

# Align the query to the reference (output sam).
$(aln).sam: $(ref).fasta.tai $(query)
	thermite align -a $(aln_flags) -o $@ $^

# Align the query to the reference (output bam).
$(aln).bam: $(ref).fasta.tai $(query)
	thermite align -a $(aln_flags) -o $@ $^

# Get STAR comparison metrics.
chr$(chr)_comparison_metrics.txt: $(aln).bam $(compare)
	python3 ../lib/testing/alignment_metrics.py $(compare) $< > $@
