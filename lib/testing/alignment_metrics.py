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
        description="Get alignment metrics. Assumes reads are in same order in sam/bam file."
    )
    parser.add_argument("in1", help="input sam/bam file 1 (cellranger)")
    parser.add_argument("in2", help="input sam/bam file 2 (thermite)")
    parser.add_argument("-r", "--ref", help="reference fasta file")
    parser.add_argument(
        "-p",
        "--print",
        nargs="+",
        choices=["all", "identical", "overlap", "gene", "unmapped"],
        help="when should mismatching alignments be printed",
    )
    args = parser.parse_args()

    when_print = set(args.print) if args.print else set()
    reader1, reader1_type = get_alignment_reader(args.in1)
    reader2, reader2_type = get_alignment_reader(args.in2)
    ref_file = pysam.FastaFile(args.ref) if args.ref else None
    metrics = Metrics()

    for row1 in reader1:
        print_aln = "all" in when_print

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
        # while row1.is_secondary:
        #     row1 = next(reader1)
        # while row2.is_secondary:
        #     row2 = next(reader2)

        row2s = []
        while True:
            row2 = next(reader2)
            row2s.append(row2)
            if row2.is_unmapped or row2.get_tag("HI") == row2.get_tag("NH"):
                break

        if row1.query_name != row2s[0].query_name:
            print(
                f"query names not matching up: {row1.query_name}, {row2s[0].query_name}"
            )
            exit()

        metrics.n_in1_identical_align += sam_query_identical_to_reference(row1)
        metrics.n_in2_identical_align += sam_query_identical_to_reference(row2s[0])

        metrics.n_in1_unaligned += row1.is_unmapped
        metrics.n_in2_unaligned += row2s[0].is_unmapped
        if "unmapped" in when_print and row1.is_unmapped != row2s[0].is_unmapped:
            print_aln = True

        metrics.n_same_chromosome_align += row1.reference_name in (
            r.reference_name for r in row2s
        )

        aln_overlap = any((queries_overlap(row1, r) for r in row2s))
        metrics.n_overlapping_align += aln_overlap
        if "overlap" in when_print and not aln_overlap:
            print_aln = True

        aln_identical = any((queries_identical(row1, r) for r in row2s))
        metrics.n_identical_align += aln_identical
        if "identical" in when_print and not aln_identical:
            print_aln = True

        if row1.has_tag("GX"):
            same_gene = len(get_gx_tags([row1]) & get_gx_tags(row2s)) > 0
            metrics.n_same_gene_align += same_gene
            metrics.n_reads_on_genes += 1
            if "gene" in when_print and not same_gene:
                print_aln = True

        if print_aln:
            if ref_file and not row1.is_unmapped:
                s = ref_file.fetch(
                    row1.reference_name, row1.reference_start, row1.reference_end
                )
                print("ref:", s)

            print(row1.tostring(reader1))
            print()

            for row2 in row2s:
                if ref_file and not row2.is_unmapped:
                    s = ref_file.fetch(
                        row2.reference_name, row2.reference_start, row2.reference_end
                    )
                    print("ref:", s)

                print(row2.tostring(reader2))

            print()
            print()

    print(f"file1: {args.in1}, file2: {args.in2}")
    metrics_to_markdown(metrics)


def get_gx_tags(rows: list) -> set:
    gxs = set()
    for r in rows:
        if r.has_tag("GX"):
            gxs.update(r.get_tag("GX").split(";"))
    return gxs


def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGNTacgnt", "TGCNAtgcna"))[::-1]


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


def queries_identical(
    alignment1: pysam.AlignedSegment, alignment2: pysam.AlignedSegment
) -> int:
    return (
        alignment1.reference_name == alignment2.reference_name
        and alignment1.reference_start == alignment2.reference_start
        and alignment1.reference_end == alignment2.reference_end
        and alignment1.is_reverse == alignment2.is_reverse
    )


def queries_overlap(alignment1: pysam.AlignedSegment, alignment2: pysam.AlignedSegment):
    return (
        alignment1.reference_name == alignment2.reference_name
        and alignment1.is_reverse == alignment2.is_reverse
        and (
            (
                alignment1.reference_end > alignment2.reference_start
                and alignment1.reference_start < alignment2.reference_end
            )
            or (
                alignment1.reference_start < alignment2.reference_end
                and alignment1.reference_end > alignment2.reference_start
            )
        )
    )


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


def metrics_to_markdown(metrics: dataclass):

    print(
        f"""
|metric|fraction|
|------|--------|
|file1 identical alignment to ref | {round(metrics.n_in1_identical_align / metrics.n_reads, 3)}|
|file2 identical alignment to ref | {round(metrics.n_in2_identical_align / metrics.n_reads, 3)}|
|file1 unaligned reads |            {round(metrics.n_in1_unaligned / metrics.n_reads, 3)}|
|file2 unaligned reads |            {round(metrics.n_in2_unaligned / metrics.n_reads, 3)}|
|file1 and file2 reads on same chr |{round(metrics.n_same_chromosome_align / metrics.n_reads, 3)}|
|file1 and file2 identical alignments |{round(metrics.n_identical_align / metrics.n_reads, 3)}|
|file1 and file2 overlapping align |{round(metrics.n_overlapping_align / metrics.n_reads, 2)}|
|file1 and file2 reads on same gene |{round(metrics.n_same_gene_align / metrics.n_reads_on_genes, 3)}|
    """
    )


@dataclass
class Metrics:
    n_reads: int = 0
    n_reads_on_genes: int = 0
    n_in1_identical_align: int = 0
    n_in2_identical_align: int = 0
    n_in1_unaligned: int = 0
    n_in2_unaligned: int = 0
    n_overlapping_align: int = 0
    n_same_gene_align: int = 0
    n_same_chromosome_align: int = 0
    n_identical_align: int = 0


if __name__ == "__main__":
    main()
