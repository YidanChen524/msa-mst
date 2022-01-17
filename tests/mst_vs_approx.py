"""
comparing the running time and sp_score and length of alignments of msa-mst using prim_mst and approx
"""
import sys
import time
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../seqs'))
from seqs import Seqs


num_of_seq = list(range(10, 301, 10))
seq_len = 10
time_mst = []
time_approx = []
score_mst = []
score_approx = []
len_alm_mst = []
len_alm_approx = []

for n in num_of_seq:
    fname = f"related_seqs/test_{n}_{seq_len}.fa"
    seqs = Seqs(fname)
    # approx
    start = time.perf_counter()
    seqs.global_align("approx")
    end = time.perf_counter()
    time_approx.append(end - start)
    score_approx.append(seqs.sp_score())
    len_alm_approx.append(seqs.length_of_alignment())
    # mst list
    start = time.perf_counter()
    seqs.global_align("prim_mst", "list")
    end = time.perf_counter()
    time_mst.append(end - start)
    score_mst.append(seqs.sp_score())
    len_alm_mst.append(seqs.length_of_alignment())

plt.figure(1)
plt.plot(num_of_seq, time_approx, label="approx", alpha=0.7)
plt.plot(num_of_seq, time_mst, label="prim_mst", alpha=0.7)
plt.legend(loc="upper left")
plt.title(f"running time vs number of sequences (seq_len={seq_len})")
plt.xlabel("number of sequences")
plt.ylabel("running time")
plt.savefig("running time vs number of sequences")

plt.figure(2)
plt.plot(num_of_seq, score_approx, label="approx", alpha=0.7)
plt.plot(num_of_seq, score_mst, label="prim_mst", alpha=0.7)
plt.legend(loc="upper left")
plt.title(f"sp score vs number of sequences (seq_len={seq_len})")
plt.xlabel("number of sequences")
plt.ylabel("sp score")
plt.savefig("sp score vs number of sequences")

plt.figure(3)
plt.plot(num_of_seq, len_alm_approx, label="approx", alpha=0.7)
plt.plot(num_of_seq, len_alm_mst, label="prim_mst", alpha=0.7)
plt.legend(loc="upper left")
plt.title(f"alignment length vs number of sequences (seq_len={seq_len})")
plt.xlabel("number of sequences")
plt.ylabel("length of alignment")
plt.savefig("length of alignment vs number of sequences")
