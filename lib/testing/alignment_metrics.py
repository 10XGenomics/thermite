import csv
import os
import pysam
import argparse
from collections.abc import Iterable
from collections import namedtuple
from dataclasses import dataclass

# TODO: use data class instead of named tuple

PAF_TUPLE = namedtuple(
    "paf_record",
    [
        "query_name",
        "query_len",
        "query_start",
        "query_end",
        "strand",
        "target_name",
        "target_len",
        "target_start",
        "target_end",
        "num_match_residues",
        "alignment_len",
        "mapping_qual",
    ],
)


def main():
    parser = argparse.ArgumentParser(
        description="Get alignment metrics. Assumes reads are in same order in paf file."
    )
    parser.add_argument("in1", help="input sam/bam file 1")
    parser.add_argument("in2", help="input sam/bam file 2")
    args = parser.parse_args()

    reader1, reader1_type = get_alignment_reader(args.in1)
    reader2, reader2_type = get_alignment_reader(args.in2)
    metrics = Metrics()

    for row1, row2 in zip(
        reader1,
        reader2,
    ):
        metrics.n_reads += 1
        # if reader1_type == "paf":
        #     print("Not supporting PAF atm")
        #     exit()
        #     while row1[0].startswith("@"):
        #         row1 = next(reader1)
        # if reader2_type == "paf":
        #     print("Not supporting PAF atm")
        #     exit()
        #     while row2[0].startswith("@"):
        #         row2 = next(reader2)
        while row1.is_secondary:
            row1 = next(reader1)
        while row2.is_secondary:
            row2 = next(reader2)

        if row1.query_name != row2.query_name:
            print(f"query names not matching up:{row1.query_name}, {row2.query_name}")
            exit()

        metrics.n_in1_identical_align += sam_query_identical_to_reference(row1)
        metrics.n_in2_identical_align += sam_query_identical_to_reference(row2)
        metrics.n_in1_unaligned += row1.is_unmapped
        metrics.n_in2_unaligned += row2.is_unmapped
        metrics.n_same_chromosome_align += row1.reference_name == row2.reference_name
        try:
            metrics.n_same_gene_align += row1.get_tag("GX") == row2.get_tag("GX")
            metrics.n_reads_on_genes += 1
        except KeyError:
            pass

        metrics.n_concordant_align += (
            1 if row1.compare(row2) == 0 else 0
        )  # will this work?
        # should just compare chromosome, start position, end position, and strand
    print(f"file1: {args.in1}, file2: {args.in2}")
    print(
        f"file1 identical alignment to ref fraction: {metrics.n_in1_identical_align/metrics.n_reads}"
    )
    print(
        f"file2 identical alignment to ref fraction: {metrics.n_in2_identical_align/metrics.n_reads}"
    )
    print(f"file1 unaligned reads fraction: {metrics.n_in1_unaligned/metrics.n_reads}")
    print(f"file2 unaligned reads fraction: {metrics.n_in2_unaligned/metrics.n_reads}")
    print(
        f"file1 and file2 reads on same chr fraction: {metrics.n_same_chromosome_align/metrics.n_reads}"
    )

    print(
        f"file1 and file2 identical alignments fraction: {metrics.n_concordant_align/metrics.n_reads}"
    )
    print(
        f"file1 and file2 reads on same gene fraction: {metrics.n_same_gene_align/metrics.n_reads_on_genes}"
    )


def get_alignment_reader(path: str) -> Iterable:
    _, ext = os.path.splitext(path)
    if ext == ".bam" or ext == ".sam":
        pysam.sort("-n", "-o", f"namesorted_{path}", path)
        samfile = pysam.AlignmentFile(f"namesorted_{path}", "rb")
        return (samfile.fetch(until_eof=True), "sam")
    elif ext == ".paf":
        return (csv.reader(open(path), delimiter="\t"), "paf")


def paf_query_identical_to_reference(alignment: list) -> int:
    return alignment.query_len == alignment.alignment_len


def sam_query_identical_to_reference(
    alignment: pysam.AlignedSegment,
) -> int:
    return alignment.query_alignment_length == alignment.reference_length


def sam_query_unmapped(alignment: pysam.AlignedSegment) -> int:
    if alignment.is_unmapped:
        return


def parse_paf_alignment(alignment: list) -> namedtuple:
    return PAF_TUPLE(*alignment)


def convert_sam_to_paf_alignment(alignment: pysam.AlignedSegment) -> namedtuple:
    strand = "-" if alignment.is_reverse else "+"
    return PAF_TUPLE(
        alignment.query_name,
        str(alignment.infer_query_length()),
        str(alignment.query_alignment_start),
        str(alignment.query_alignment_start + alignment.infer_query_length() - 1),
        strand,
        alignment.reference_name,
        alignment.reference_length,
        str(alignment.reference_start),
        str(alignment.reference_end - 1),
        str(len(alignment.get_aligned_pairs())),
        str(alignment.query_alignment_length),
        str(alignment.mapping_quality),
    )


def compare_sam_to_paf(
    sam_alignment: pysam.AlignedSegment, paf_alignment: namedtuple
) -> int:
    strand = "-" if sam_alignment.is_reverse else "+"
    comparison = (
        (sam_alignment.query_name == paf_alignment.query_name)
        and (str(sam_alignment.infer_query_length()) == paf_alignment.query_len)
        and (str(sam_alignment.query_alignment_start) == paf_alignment.query_start)
        and (
            str(
                sam_alignment.query_alignment_start
                + sam_alignment.infer_query_length()
                - 1
            )
            == paf_alignment.query_end
        )
        and (strand == paf_alignment.strand)
        and (sam_alignment.reference_name == paf_alignment.target_name)
        and (str(sam_alignment.reference_start) == paf_alignment.target_start)
        and (str(sam_alignment.reference_end - 1) == paf_alignment.target_end)
        and (
            str(len(sam_alignment.get_aligned_pairs()))
            == paf_alignment.num_match_residues
        )
        and (str(sam_alignment.mapping_quality) == paf_alignment.mapping_qual)
    )
    return comparison


@dataclass
class Metrics:
    n_reads: int = 0
    n_reads_on_genes: int = 0
    n_in1_identical_align: int = 0
    n_in2_identical_align: int = 0
    n_concordant_align: int = 0
    n_in1_unaligned: int = 0
    n_in2_unaligned: int = 0
    n_overlapping_align: int = 0
    n_same_gene_align: int = 0
    n_same_chromosome_align: int = 0
    n_identical_align: int = 0


if __name__ == "__main__":
    main()
