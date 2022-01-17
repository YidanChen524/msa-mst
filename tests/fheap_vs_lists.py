"""
comparing the running time of msa-mst using list vs using fibonacci heap as the number of sequences increases
"""
import sys
import time
import os
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), '../seqs'))
from seqs import Seqs


# time vs number of sequences
num_of_seq = list(range(10, 1001, 10))
seq_len = 10
time_list = []
time_fheap = []

for n in num_of_seq:
    fname = f"random_seqs/test_{n}_{seq_len}.fa"
    seqs = Seqs(fname)
    # list
    start = time.perf_counter()
    seqs.build_mst("list")
    end = time.perf_counter()
    time_list.append(end - start)
    # fheap
    start = time.perf_counter()
    seqs.build_mst("fheap")
    end = time.perf_counter()
    time_fheap.append(end - start)

plt.plot(num_of_seq, time_list, label="list", alpha=0.7)
plt.plot(num_of_seq, time_fheap, label="fheap", alpha=0.7)
plt.legend(loc="upper left")
plt.title(f"running time vs number of sequences (length={seq_len})")
plt.xlabel("number of sequences")
plt.ylabel("running time")
plt.savefig("running time vs number of sequences")
