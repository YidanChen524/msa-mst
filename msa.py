import sys
import argparse
from Bio import SeqIO
from graph import Graph


# define default gap and score matrix
gap = 5
score = [[0, 5, 2, 5],
         [5, 0, 5, 2],
         [2, 5, 0, 5],
         [5, 2, 5, 0]]


# define how nucleotides are mapped to indices
mapping = {'A': 0, 'C': 1, 'G': 2, 'T': 3,
           'a': 0, 'c': 1, 'g': 2, 't': 3,
           'N':0, 'R':0, 'S':1}


def pairwise(s1, s2, alignment=False):

    # fill out T table
    m = len(s1) + 1
    n = len(s2) + 1
    T = [[None for j in range(n)] for i in range(m)]

    for i in range(m):
        for j in range(n):
            v0 = v1 = v2 = v3 = None
            if i == 0 and j == 0:
                v0 = 0
            if i > 0 and j >= 0:
                v1 = T[i-1][j] + gap
            if i >= 0 and j > 0:
                v2 = T[i][j-1] + gap
            if i > 0 and j > 0:
                v3 = T[i-1][j-1] + score[mapping[s1[i-1]]][mapping[s2[j-1]]]
            T[i][j] = min([v for v in [v0, v1, v2, v3] if v is not None])

    # pairwise score
    if not alignment:
        return T[m-1][n-1]

    # alignments
    i = m - 1
    j = n - 1
    a1 = a2 = ''
    while i > 0 or j > 0:
        v = T[i][j]
        if i > 0 and j > 0 and v == T[i-1][j-1] + score[mapping[s1[i-1]]][mapping[s2[j-1]]]:
            a1 = s1[i-1] + a1
            a2 = s2[j-1] + a2
            i -= 1
            j -= 1
        elif i > 0 and v == T[i-1][j] + gap:
            a1 = s1[i-1] + a1
            a2 = '-' + a2
            i -= 1
        elif j > 0 and v == T[i][j-1] + gap:
            a2 = s2[j-1] + a2
            a1 = '-' + a1
            j -= 1

    return [a1, a2]


def exact():
    """run the exact algorithm up to 3 sequences"""
    pass


def approx():
    """run the approx algorithm"""
    pass


def mst(seqs, data_structure="list"):
    """run the mst algorithm"""
    # calculate the distance matrix between sequences
    m = len(seqs)
    G = [[0 for i in range(m)] for j in range(m)]
    for i in range(m):
        for j in range(m):
            G[i][j] = pairwise(seqs[i], seqs[j])
    # run mst algorithm and return a mst
    g = Graph(G)
    print(G)
    print(g.prim_list())
    print(g.prim_fibonacci_heap())
    # align the sequences based on mst


if __name__ == "__main__":
    seqs = ['GTTCCGAAAGGCTAGCGCTAGGCGCC',
            'ATGGATTTATCTGCTCTTCG',
            'TGCATGCTGAAACTTCTCAACCA']
    mst(seqs)

