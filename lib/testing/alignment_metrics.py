import csv
import argparse
from collections import namedtuple


def main():
    parser = argparse.ArgumentParser(
        description="Get alignment metrics. Assumes reads are in same order in paf file."
    )
    parser.add_argument("inpaf1", help="input paf file 1")
    parser.add_argument("inpaf2", help="input paf file 2")
    args = parser.parse_args()

    paf_file_path1 = args.inpaf1
    paf_file_path2 = args.inpaf2

    paf_file1_reader = csv.reader(open(paf_file_path1), delimiter="\t")
    paf_file2_reader = csv.reader(open(paf_file_path2), delimiter="\t")

    n_reads = 0
    n_paf1_identical_align = 0
    n_paf2_identical_align = 0
    n_concordant_align = 0

    for row1, row2 in zip(paf_file1_reader, paf_file2_reader):
        n_reads += 1
        n_paf1_identical_align += query_identical_to_reference(row1[:-1])
        n_paf2_identical_align += query_identical_to_reference(row2[:-1])
        print(row1, row2)
        if row1 == row2:
            n_concordant_align += 1
    print(f"paf1 identical alignment to ref fraction: {n_paf1_identical_align/n_reads}")
    print(f"paf2 identical alignment to ref fraction: {n_paf2_identical_align/n_reads}")
    print(f"paf1 and paf2 identical alignments fraction: {n_concordant_align/n_reads}")


def query_identical_to_reference(alignment: list) -> int:
    alignment = parse_alignment(alignment)
    query_length = alignment.query_len
    align_block_length = alignment.alignment_len
    if query_length == align_block_length:
        return 1
    else:
        return 0


def parse_alignment(alignment: list) -> namedtuple:
    align_tuple = namedtuple(
        "alignment",
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
    return align_tuple(*alignment)


if __name__ == "__main__":
    main()