import csv
import os
import pysam
import argparse
from collections.abc import Iterable
from collections import namedtuple

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

ALIGN_TUPLE = namedtuple(
    "alignment",
    [
        "read_name",
        "query_len",
        "query_start",
        "query_end",
        "target_name",
        "target_start",
        "target_end",
        "num_match_residue",
    ],
)


def main():
    parser = argparse.ArgumentParser(
        description="Get alignment metrics. Assumes reads are in same order in paf file."
    )
    parser.add_argument("in1", help="input sam/bam/paf file 1")
    parser.add_argument("in2", help="input sam/bam/paf file 2")
    args = parser.parse_args()

    reader1, reader1_type = get_alignment_reader(args.in1)
    reader2, reader2_type = get_alignment_reader(args.in2)

    n_reads = 0
    n_in1_identical_align = 0
    n_in2_identical_align = 0
    n_concordant_align = 0

    for row1, row2 in zip(reader1, reader2):
        n_reads += 1
        if reader1_type == "paf":
            row1 = parse_paf_alignment(row1[:-1])
            n_in1_identical_align += paf_query_identical_to_reference(row1)
            print(row1)
        else:
            n_in1_identical_align += sam_query_identical_to_reference(row1)
            row1 = convert_sam_to_paf_alignment(row1)
        if reader2_type == "paf":
            row2 = parse_paf_alignment(row2[:-1])
            n_in2_identical_align += paf_query_identical_to_reference(row2)
        else:
            n_in2_identical_align += sam_query_identical_to_reference(row2)
            row2 = convert_sam_to_paf_alignment(row2)
            print(row2)

        n_concordant_align += row1 == row2

    print(f"paf1 identical alignment to ref fraction: {n_in1_identical_align/n_reads}")
    print(f"paf2 identical alignment to ref fraction: {n_in2_identical_align/n_reads}")
    print(f"paf1 and paf2 identical alignments fraction: {n_concordant_align/n_reads}")


def get_alignment_reader(path: str) -> Iterable:
    _, ext = os.path.splitext(path)
    if ext == ".bam" or ext == ".sam":
        samfile = pysam.AlignmentFile(path)
        return (samfile.fetch(), "sam")  # will this work?
    elif ext == ".paf":
        return (csv.reader(open(path), delimiter="\t"), "paf")


def paf_query_identical_to_reference(alignment: list) -> int:
    if alignment.query_len == alignment.alignment_len:
        return 1
    else:
        return 0


def sam_query_identical_to_reference(
    alignment: pysam.AlignedSegment,
) -> int:
    if alignment.query_alignment_length == alignment.reference_length:
        return 1
    else:
        return 0


def parse_paf_alignment(alignment: list) -> namedtuple:
    return PAF_TUPLE(*alignment)


def convert_sam_to_paf_alignment(alignment: pysam.AlignedSegment) -> namedtuple:
    strand = "-" if alignment.is_reverse else "+"
    return PAF_TUPLE(
        alignment.query_name,
        str(alignment.infer_query_length()),
        str(alignment.query_alignment_start),
        str(alignment.query_alignment_start + alignment.query_length),
        strand,
        alignment.reference_name,
        str(alignment.reference_end),
        str(alignment.reference_start),
        str(alignment.reference_end - 1),
        len(alignment.get_aligned_pairs()),
        str(alignment.query_alignment_length),
        str(alignment.mapping_quality),
    )


if __name__ == "__main__":
    main()
