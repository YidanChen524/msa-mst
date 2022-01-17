"""
test if how we pick the first vector affects sp score
"""
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../seqs'))
from seqs import Seqs


num_of_seq = list(range(10, 301, 10))
seq_len = 10
sp_score_first = []
sp_score_random = []

for n in num_of_seq:
    fname = f"random_seqs/test_{n}_{seq_len}.fa"
    # unsorted
    seqs1 = Seqs(fname)
    seqs1.global_align("prim_mst", "list")
    sp_score_first.append(seqs1.sp_score())
    # sorted
    seqs2 = Seqs(fname)
    seqs2.global_align("prim_mst", "list", "sorted")
    sp_score_random.append(seqs2.sp_score())

plt.plot(num_of_seq, sp_score_first, label="first", alpha=0.7)
plt.plot(num_of_seq, sp_score_random, label="random", alpha=0.7)
plt.legend(loc="upper left")
plt.title(f"sp score vs number of sequences (length={seq_len})")
plt.xlabel("number of sequences")
plt.ylabel("sp score")
plt.savefig("sp score vs number of sequences")
