import time
import random
import matplotlib.pyplot as plt
from seqs import Seqs

# seqs with different length
time_approx = []
time_prim_mst_list = []
time_prim_mst_fheap = []
m = 5
for j in range(10, 101, 10):
    names = [None] * m
    seqs = [None] * m
    for i in range(m):
        names[i] = "seq " + str(i+1)
        seqs[i] = "".join(random.choices(['A', 'T', 'C','G'], k=j))
    # approx
    seqs_approx = Seqs(names, seqs)
    start_approx = time.time()
    seqs_approx.global_align("approx")
    end_approx = time.time()
    time_approx.append(end_approx-start_approx)
    # prim mst list
    seqs_prim_mst_list = Seqs(names, seqs)
    start_prim_mst_list = time.time()
    seqs_prim_mst_list.global_align("prim_mst", "list")
    end_prim_mst_list = time.time()
    time_prim_mst_list.append(end_prim_mst_list-start_prim_mst_list)
    # prim mst fheap
    seqs_prim_mst_fheap = Seqs(names, seqs)
    start_prim_mst_fheap = time.time()
    seqs_prim_mst_fheap.global_align("prim_mst", "fheap")
    end_prim_mst_fheap = time.time()
    time_prim_mst_fheap.append(end_prim_mst_fheap-start_prim_mst_fheap)

print(time_approx)
print(time_prim_mst_list)
print(time_prim_mst_fheap)



