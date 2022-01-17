"""
test if sorting the edges of the guided tree before alignments improves sp score
"""
import sys
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../seqs'))
from seqs import Seqs


num_of_seq = list(range(10, 101, 10))
seq_len = 20
sp_score_unsorted = []
sp_score_sorted = []
sp_score_reverse = []

for n in num_of_seq:
    fname = f"related_seqs/test_{n}_{seq_len}.fa"
    # unsorted
    seqs1 = Seqs(fname)
    seqs1.global_align("prim_mst", "list")
    sp_score_unsorted.append(seqs1.sp_score())
    # sorted
    seqs2 = Seqs(fname)
    seqs2.global_align("prim_mst", "list", "sorted")
    sp_score_sorted.append(seqs2.sp_score())
    # reverse
    seqs3 = Seqs(fname)
    seqs3.global_align("prim_mst", "list", "reverse")
    sp_score_reverse.append(seqs3.sp_score())

plt.plot(num_of_seq, sp_score_unsorted, label="unsorted", alpha=0.7)
plt.plot(num_of_seq, sp_score_sorted, label="sorted", alpha=0.7)
plt.plot(num_of_seq, sp_score_reverse, label="reverse", alpha=0.7)
plt.legend(loc="upper left")
plt.title(f"sp score vs number of sequences (length={seq_len})")
plt.xlabel("number of sequences")
plt.ylabel("sp score")
plt.savefig("sp score vs number of sequences")
