import sys
from helpers import parse_fasta
from seqs import Seqs


def main():
    # read in the sequences from the fasta file into a Seqs class
    seqs = Seqs(*parse_fasta(sys.argv[1]))
    # 2 approx
    # seqs.global_align("approx")
    # prim mst
    # seqs.global_align("prim_mst", "list")
    seqs.global_align("prim_mst", "fheap")
    seqs.output_alignments()


if __name__ == "__main__":
    main()
