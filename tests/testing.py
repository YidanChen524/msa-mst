import os
import sys
import time
sys.path.append(os.path.join(os.path.dirname(__file__), '../src'))
from helpers import parse_fasta
from classes.seqs import Seqs


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


def compare():
    time_approx = [0] * 100
    score_approx = [0] * 100
    time_prim_mst_list = [0] * 100
    score_prim_mst_list = [0] * 100
    time_prim_mst_fheap = [0] * 100
    score_prim_mst_fheap = [0] * 100
    count = 0
    for i in range(10, 101, 10):
        for j in range(10, 101, 10):
            test_file = f"test_seqs/test_{i}_{j}.fa"
            time_approx[count], score_approx[count] = msa_time_score(test_file, "approx")
            time_prim_mst_list[count], score_prim_mst_list[count] = msa_time_score(test_file, "prim_mst", "list")
            time_prim_mst_fheap[count], score_prim_mst_fheap[count] = msa_time_score(test_file, "prim_mst", "fheap")
            count += 1
    with open("result.txt", "w") as f:
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
    compare()
    # test()
