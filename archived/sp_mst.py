import sys
from Bio import SeqIO
from graph import Graph
from config import gap, score, mapping
from utils import pairwise_distance


def sp_mst(seqs, option="list"):
    # calculate distance matrix
    n = len(seqs)
    G = [[0 for _ in range(n)] for _ in range(n)]
    for i in range(n):
        for j in range(i+1, n):
            G[i][j] = G[j][i] = pairwise_distance(seqs[i], seqs[j])
    # find the minimum spanning tree
    g = Graph(G)
    if option == "list":
        tree = g.prim_list()
    else:
        tree = g.prim_fibonacci_heap()
    # multiple sequence alignment based on mst
    tree.sort(key=lambda e: e[2])
    group = [None] * n
    count = 0
    alignments = []
    for edge in tree:
        ind1 = edge[0]
        ind2 = edge[1]
        if not group[ind1] and group[ind2]:
            group[ind1] = group[ind2] = count
            count += 1

    # return the alignment
    return tree


if __name__ == "__main__":
    # read in sequences
    names = []
    seqs = []
    for seq in SeqIO.parse(sys.argv[1], "fasta"):
        names.append(str(seq.id))
        seqs.append(str(seq.seq))
    # run exact algorithm
    optimal_alignment = sp_mst(seqs, "fibonacci_heap")
    print(optimal_alignment)
    ## print results
    #for i in range(len(seqs)):
    #    print(f">{names[i]}")
    #    print(optimal_alignment[i] + "\n")
    ## write results to file
    #with open("sp_approx_output.txt", "w") as f:
    #    for i in range(len(seqs)):
    #        f.write(f">{names[i]}" + "\n")
    #        f.write(optimal_alignment[i] + "\n\n")
