import os
import sys
import time
from ..helpers import parse_fasta
from ..classes.seqs import Seqs


def msa_timer(seqs, *options):
    start = time.perf_counter()
    seqs.global_align(*options)
    end = time.perf_counter()
    return end - start


def msa_output(fname, *options):
    seqs = Seqs(*parse_fasta(fname))
    t = msa_timer(seqs, *options)
    print(";Score: ", seqs.sp_score())
    print(";Time: ", t)
    seqs.output_alignments()
    print("\n")


def msa_time_score(fname, *options):
    seqs = Seqs(*parse_fasta(fname))
    t = msa_timer(seqs, *options)
    s = seqs.sp_score()
    return t, s


def test():
    # 2 approx
    print("##########2 Approx###########")
    msa_output(sys.argv[1], "approx")

    # prim mst list
    print("#####Prim MST with List######")
    msa_output(sys.argv[1], "prim_mst", "list")

    # prim mist fheap
    print("#####Prim MST with Fheap#####")
    msa_output(sys.argv[1], "prim_mst", "fheap")


def compare(n, m, fname):
    n_test = len(n) * len(m)
    time_approx = [0] * n_test
    score_approx = [0] * n_test
    time_prim_mst_list = [0] * n_test
    score_prim_mst_list = [0] * n_test
    time_prim_mst_fheap = [0] * n_test
    score_prim_mst_fheap = [0] * n_test
    count = 0
    for i in n:
        for j in m:
            test_file = f"{os.path.dirname(os.path.abspath(__file__))}/test_seqs/test_{i}_{j}.fa"
            time_approx[count], score_approx[count] = msa_time_score(test_file, "approx")
            time_prim_mst_list[count], score_prim_mst_list[count] = msa_time_score(test_file, "prim_mst", "list")
            time_prim_mst_fheap[count], score_prim_mst_fheap[count] = msa_time_score(test_file, "prim_mst", "fheap")
            count += 1
    with open(fname, "w") as f:
        f.write("Approx\n")
        f.write(str(time_approx))
        f.write("\n")
        f.write(str(score_approx))
        f.write("\n\n")

        f.write("Prim mst list\n")
        f.write(str(time_prim_mst_list))
        f.write("\n")
        f.write(str(score_prim_mst_list))
        f.write("\n\n")

        f.write("Prim mst fheap\n")
        f.write(str(time_prim_mst_fheap))
        f.write("\n")
        f.write(str(score_prim_mst_fheap))
        f.write("\n\n")


if __name__ == "__main__":
    compare(list(range(10, 301, 10)), [10], "results/result_m_10.txt")
