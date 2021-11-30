from ..config import gap, score, mapping
from ..helpers import string_concat
from .graph import Graph


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
        self.guide_tree = None
        self.distance_matrix()

    def __repr__(self):
        d = {}
        for i in range(self.len):
            d[self.names[i]] = self.seqs[i]
        return repr(d)

    def distance_matrix(self):
        for i in range(self.len):
            for j in range(i + 1, self.len):
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
                    v1 = T[i - 1][j] + gap
                if i >= 0 and j > 0:
                    v2 = T[i][j - 1] + gap
                if i > 0 and j > 0:
                    v3 = T[i - 1][j - 1] + score[mapping[seq1[i - 1]]][mapping[seq2[j - 1]]]
                T[i][j] = min([v for v in [v0, v1, v2, v3] if v is not None])
        return T[m - 1][n - 1]

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
                    v1 = T[i - 1][j] + gap
                if i >= 0 and j > 0:
                    v2 = T[i][j - 1] + gap
                if i > 0 and j > 0:
                    v3 = T[i - 1][j - 1] + score[mapping[seq1[i - 1]]][mapping[seq2[j - 1]]]
                T[i][j] = min([v for v in [v0, v1, v2, v3] if v is not None])
        # alignments
        i = m - 1
        j = n - 1
        a1 = a2 = ''
        while i > 0 or j > 0:
            v = T[i][j]
            if i > 0 and j > 0 and v == T[i - 1][j - 1] + score[mapping[seq1[i - 1]]][mapping[seq2[j - 1]]]:
                a1 = string_concat(seq1[i - 1], a1)
                a2 = string_concat(seq2[j - 1], a2)
                i -= 1
                j -= 1
            elif i > 0 and v == T[i - 1][j] + gap:
                a1 = string_concat(seq1[i - 1], a1)
                a2 = string_concat('-', a2)
                i -= 1
            elif j > 0 and v == T[i][j - 1] + gap:
                a2 = string_concat(seq2[j - 1], a2)
                a1 = string_concat('-', a1)
                j -= 1
        return a1, a2

    def add_to_alignments(self, src, dest):
        """align seqs[ind] to the alignments so far"""
        a1, a2 = self.pairwise_alignment(src, dest)
        ptr = self.alignments
        # initialize alignments nodes
        if not self.alignments.next:
            for i in range(len(a1)):
                ptr.next = AlmNode(self.len)
                ptr = ptr.next
                ptr.val[src] = a1[i]
                ptr.val[dest] = a2[i]
        # when msa node has not been aligned yet (for sorted guide tree)
        elif not self.if_aligned[src]:
            pass
        # align dest node based on msa node
        else:
            for i in range(len(a1)):
                while ptr.next and ptr.next.val[src] != a1[i] and ptr.next.val[src] == '-':
                    ptr = ptr.next
                if not ptr.next or ptr.next.val[src] != a1[i]:
                    node = AlmNode(self.len)
                    node.next = ptr.next
                    ptr.next = node
                ptr = ptr.next
                ptr.val[dest] = a2[i]
        self.if_aligned[src] = self.if_aligned[dest] = True

    def global_align(self, option, option2=None, option3=None):
        if option == "prim_mst":
            g = Graph(self.D)
            if option2 == "list":
                self.guide_tree = g.prim_list_mst()
            elif option2 == "fheap":
                self.guide_tree = g.prim_fheap_mst()
            if option3 == "sorted":
                self.guide_tree.sort(key=lambda x: x[2])
            for src, dest, _ in self.guide_tree:
                self.add_to_alignments(src, dest)
        if option == "approx":
            _, center = min((sum(val), ind) for (ind, val) in enumerate(self.D))
            self.guide_tree = [None] * (self.len - 1)
            pos = 0
            for i in range(self.len):
                if i != center:
                    self.add_to_alignments(center, i)
                    self.guide_tree[pos] = (center, i, None)
                    pos += 1

    def output_alignments(self, save=False):
        for i in range(self.len):
            ptr = self.alignments
            a = ''
            while ptr.next:
                a = string_concat(a, ptr.next.val[i])
                ptr = ptr.next
            print(f">{self.names[i]}")
            print(a, '\n')

    def test_alignments(self):
        """
        test if alignments are correct:
        - no positions with all gaps
        - induced pairwise alignments equals pairwise alignments when stripped extra gaps
        """
        # check if there are positions with all gaps
        ptr = self.alignments
        pos = 0
        while ptr.next:
            if "".join(ptr.next.val) == '-' * self.len:
                return f"Test Fail: All Gaps at Position {pos}"
            else:
                ptr = ptr.next
                pos += 1
        # check if induced pairwise alignments equals pairwise alignments when stripped extra gaps
        for i, j, _ in self.guide_tree:
            a1, a2 = self.pairwise_alignment(i, j)
            msa1 = msa2 = ""
            ptr = self.alignments
            while ptr.next:
                if ptr.next.val[i] != '-' or ptr.next.val[j] != '-':
                    msa1 += ptr.next.val[i]
                    msa2 += ptr.next.val[j]
                ptr = ptr.next
            if a1 != msa1 or a2 != msa2:
                return "Test Fail: Induced Pairwise Alignment not Correct"

        return "Test Success"

    def sp_score(self):
        ptr = self.alignments
        cost = 0
        while ptr.next:
            for i in range(self.len):
                for j in range(i + 1, self.len):
                    if ptr.next.val[i] == ptr.next.val[j]:
                        cost += 0
                    elif ptr.next.val[i] == '-' or ptr.next.val[j] == '-':
                        cost += gap
                    else:
                        cost += score[mapping[ptr.next.val[i]]][mapping[ptr.next.val[j]]]
            ptr = ptr.next
        return cost
