from config import gap, score, mapping
from graph import Graph
from helpers import string_concat


class AlmNode:
    def __init__(self, num_seqs=0):
        self.val = ['-'] * num_seqs
        self.next = None


class Seqs:
    def __init__(self, names, seqs):
        self.names = names
        self.seqs = seqs
        self.len = len(self.seqs)
        self.alignments = AlmNode()
        self.if_aligned = [False] * self.len
        self.D = [[0 for _ in range(self.len)] for _ in range(self.len)]
        self.distance_matrix()

    def __repr__(self):
        d = {}
        for i in range(self.len):
            d[self.names[i]] = self.seqs[i]
        return repr(d)

    def distance_matrix(self):
        for i in range(self.len):
            for j in range(i+1, self.len):
                self.D[i][j] = self.D[j][i] = self.pairwise_score(i, j)

    def pairwise_score(self, ind1: int, ind2: int) -> int:
        """calculate the optimal sp_score between 2 sequences"""
        seq1 = self.seqs[ind1]
        seq2 = self.seqs[ind2]
        # fill out T table
        m = len(seq1) + 1
        n = len(seq2) + 1
        T = [[None for _ in range(n)] for _ in range(m)]
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
                    v3 = T[i-1][j-1] + score[mapping[seq1[i-1]]][mapping[seq2[j-1]]]
                T[i][j] = min([v for v in [v0, v1, v2, v3] if v is not None])
        return T[m-1][n-1]

    def pairwise_alignment(self, ind1, ind2):
        seq1 = self.seqs[ind1]
        seq2 = self.seqs[ind2]
        # fill out T table
        m = len(seq1) + 1
        n = len(seq2) + 1
        T = [[None for _ in range(n)] for _ in range(m)]
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
                    v3 = T[i-1][j-1] + score[mapping[seq1[i-1]]][mapping[seq2[j-1]]]
                T[i][j] = min([v for v in [v0, v1, v2, v3] if v is not None])
        # alignments
        i = m - 1
        j = n - 1
        a1 = a2 = ''
        while i > 0 or j > 0:
            v = T[i][j]
            if i > 0 and j > 0 and v == T[i-1][j-1] + score[mapping[seq1[i-1]]][mapping[seq2[j-1]]]:
                a1 = string_concat(seq1[i-1], a1)
                a2 = string_concat(seq2[j-1], a2)
                i -= 1
                j -= 1
            elif i > 0 and v == T[i-1][j] + gap:
                a1 = string_concat(seq1[i-1], a1)
                a2 = string_concat('-', a2)
                i -= 1
            elif j > 0 and v == T[i][j-1] + gap:
                a2 = string_concat(seq2[j-1], a2)
                a1 = string_concat('-', a1)
                j -= 1
        return a1, a2

    def add_to_alignments(self, src, dest):
        """align seqs[ind] to the alignments so far"""
        a1, a2 = self.pairwise_alignment(src, dest)
        ptr = self.alignments
        if not self.if_aligned[src] and not self.if_aligned[dest]:
            for i in range(len(a1)):
                ptr.next = AlmNode(self.len)
                ptr = ptr.next
                ptr.val[src] = a1[i]
                ptr.val[dest] = a2[i]
        else:
            for i in range(len(a1)):
                if ptr.next.val[src] != a1[i]:
                    node = AlmNode(self.len)
                    node.next = ptr.next
                    ptr.next = node
                ptr = ptr.next
                ptr.val[dest] = a2[i]
        self.if_aligned[src] = self.if_aligned[dest] = True

    def global_align(self, option, suboption=None):
        if option == "prim_mst":
            g = Graph(self.D)
            if suboption == "list":
                mst = g.prim_list_mst()
            elif suboption == "fheap":
                mst = g.prim_fheap_mst()
            for src, dest, _ in mst:
                self.add_to_alignments(src, dest)
        if option == "approx":
            _, center = min((sum(val), ind) for (ind, val) in enumerate(self.D))
            for i in range(self.len):
                if i != center:
                    self.add_to_alignments(center, i)

    def output_alignments(self, save=False):
        for i in range(self.len):
            ptr = self.alignments
            a = ''
            while ptr.next:
                a = string_concat(a, ptr.next.val[i])
                ptr = ptr.next
            print(f">{self.names[i]}")
            print(a, '\n')

