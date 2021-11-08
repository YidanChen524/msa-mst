import sys
import time
from helpers import parse_fasta
from seqs import Seqs


def main():
    # 2 approx
    print("##########2 Approx###########")
    seqs_approx = Seqs(*parse_fasta(sys.argv[1]))
    start_approx = time.time()
    seqs_approx.global_align("approx")
    end_approx = time.time()
    print(";Score: ", seqs_approx.sp_score())
    print(";Time: ", end_approx - start_approx)
    seqs_approx.output_alignments()
    print("\n")

    # prim mst list
    print("#####Prim MST with List######")
    seqs_prim_mst_list = Seqs(*parse_fasta(sys.argv[1]))
    start_prim_mst_list = time.time()
    seqs_prim_mst_list.global_align("prim_mst", "list")
    end_prim_mst_list = time.time()
    print(";Score: ", seqs_prim_mst_list.sp_score())
    print(";Time: ", end_prim_mst_list - start_prim_mst_list)
    seqs_prim_mst_list.output_alignments()
    print("\n")

    # prim mist fheap
    print("#####Prim MST with Fheap#####")
    seqs_prim_mst_fheap = Seqs(*parse_fasta(sys.argv[1]))
    start_prim_mst_fheap = time.time()
    seqs_prim_mst_fheap.global_align("prim_mst", "fheap")
    end_prim_mst_fheap = time.time()
    print(";Score: ", seqs_prim_mst_fheap.sp_score())
    print(";Time: ", end_prim_mst_fheap - start_prim_mst_fheap)
    seqs_prim_mst_fheap.output_alignments()
    print("\n")


if __name__ == "__main__":
    main()
