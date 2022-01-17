from .config import gap, score, mapping
from .helpers import string_concat, parse_fasta
from .graph import Graph


class AlmNode:
    def __init__(self):
        self.values = {}
        self.next = None


class Seqs:
    def __init__(self, filename):
        names, seqs = parse_fasta(filename)
        self.names = names
        self.seqs = seqs
        self.len = len(self.seqs)
        self.D = self.distance_matrix()
        self.guide_tree = None
        self.alignments = None
        self.if_aligned = None
        self.aligned_group = None

    def distance_matrix(self):
        """returns distance matrix for sequences stored in Seqs object"""
        D = [[0 for _ in range(self.len)] for _ in range(self.len)]
        for i in range(self.len):
            for j in range(i + 1, self.len):
                D[i][j] = D[j][i] = self.pairwise_score(i, j)
        return D

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
        """align src and dest sequences and add them to self.alignments"""
        a1, a2 = self.pairwise_alignment(src, dest)
        if not self.if_aligned[src]:
            ptr_src = AlmNode()
            self.alignments.append(ptr_src)
            self.aligned_group[src] = len(self.alignments) - 1
            self.if_aligned[src] = True
            for i in range(len(a1)):
                ptr_src.next = AlmNode()
                ptr_src = ptr_src.next
                ptr_src.values[src] = a1[i]
        if not self.if_aligned[dest]:
            ptr_dest = AlmNode()
            self.alignments.append(ptr_dest)
            self.aligned_group[dest] = len(self.alignments) - 1
            self.if_aligned[dest] = True
            for i in range(len(a2)):
                ptr_dest.next = AlmNode()
                ptr_dest = ptr_dest.next
                ptr_dest.values[dest] = a2[i]
        # merge 2 alignment groups
        self.merge_alignments(src, dest, a1, a2)

    def merge_alignments(self, src, dest, a1, a2):
        """merge 2 alignment groups based on src and dest sequences' alignments"""
        # extract corresponding alignment groups
        src_ptr, dest_ptr = self.alignments[self.aligned_group[src]], self.alignments[self.aligned_group[dest]]
        src_keys, dest_keys = src_ptr.next.values.keys(), dest_ptr.next.values.keys()
        # put sequences in dest group into src group
        for key in dest_keys:
            self.aligned_group[key] = self.aligned_group[src]
        # update src alignments
        for i in range(len(a1)):
            while src_ptr.next and src_ptr.next.values[src] != a1[i] and src_ptr.next.values[src] == '-' \
                    and dest_ptr.next and dest_ptr.next.values[dest] != a2[i] and dest_ptr.next.values[dest] == '-':
                src_ptr = src_ptr.next
                dest_ptr = dest_ptr.next
                for key in dest_keys:
                    src_ptr.values[key] = dest_ptr.values[key]
            while src_ptr.next and src_ptr.next.values[src] != a1[i] and src_ptr.next.values[src] == '-':
                src_ptr = src_ptr.next
                for key in dest_keys:
                    src_ptr.values[key] = '-'
            while dest_ptr.next and dest_ptr.next.values[dest] != a2[i] and dest_ptr.next.values[dest] == '-':
                dest_ptr = dest_ptr.next
                new_node = AlmNode()
                for key in src_keys:
                    new_node.values[key] = '-'
                for key in dest_keys:
                    new_node.values[key] = dest_ptr.values[key]
                new_node.next = src_ptr.next
                src_ptr.next = new_node
                src_ptr = src_ptr.next
            if src_ptr.next and src_ptr.next.values[src] != a1[i] and a1[i] == '-':
                new_node = AlmNode()
                for key in src_keys:
                    new_node.values[key] = '-'
                for key in dest_keys:
                    new_node.values[key] = dest_ptr.next.values[key]
                new_node.next = src_ptr.next
                src_ptr.next = new_node
                src_ptr = src_ptr.next
                dest_ptr = dest_ptr.next
            elif dest_ptr.next and dest_ptr.next.values[dest] != a2[i] and a2[i] == '-':
                for key in dest_keys:
                    src_ptr.next.values[key] = '-'
                src_ptr = src_ptr.next
            else:
                if src_ptr.next and dest_ptr.next:
                    src_ptr = src_ptr.next
                    dest_ptr = dest_ptr.next
                    for key in dest_keys:
                        src_ptr.values[key] = dest_ptr.values[key]
        else:
            while src_ptr.next and dest_ptr.next:
                src_ptr = src_ptr.next
                dest_ptr = dest_ptr.next
                for key in dest_keys:
                    src_ptr.values[key] = dest_ptr.values[key]
            while src_ptr.next:
                src_ptr = src_ptr.next
                for key in dest_keys:
                    src_ptr.values[key] = '-'
            while dest_ptr.next:
                new_node = AlmNode()
                for key in src_keys:
                    new_node.values[key] = '-'
                src_ptr.next = new_node
                src_ptr = src_ptr.next
                dest_ptr = dest_ptr.next
                for key in dest_keys:
                    src_ptr.values[key] = dest_ptr.values[key]

    def build_mst(self, option):
        g = Graph(self.D)
        if option == "list":
            self.guide_tree = g.prim_list_mst()
        elif option == "fheap":
            self.guide_tree = g.prim_fheap_mst()

    def global_align(self, option, option2=None, option3=None, option4="first"):
        # initialize values
        self.alignments = []
        self.if_aligned = [False] * self.len
        self.aligned_group = [None] * self.len
        if option == "prim_mst":
            g = Graph(self.D)
            if option2 == "list":
                if option4 == "random":
                    self.guide_tree = g.prim_list_mst(initialize="random")
                else:
                    self.guide_tree = g.prim_list_mst(initialize="first")
            elif option2 == "fheap":
                self.guide_tree = g.prim_fheap_mst()
            if option3 == "sorted":
                self.guide_tree.sort(key=lambda x: x[2])
            elif option3 == "reverse":
                self.guide_tree.sort(key=lambda x: x[2], reverse=True)
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
        self.alignments = self.alignments[self.aligned_group[0]]

    def output_alignments(self):
        for i in range(self.len):
            ptr = self.alignments
            a = ''
            while ptr.next:
                a = string_concat(a, ptr.next.values[i])
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
            if "".join(ptr.next.values.values()) == '-' * self.len:
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
                if ptr.next.values[i] != '-' or ptr.next.values[j] != '-':
                    msa1 += ptr.next.values[i]
                    msa2 += ptr.next.values[j]
                ptr = ptr.next
            if a1 != msa1 or a2 != msa2:
                print(a1)
                print(msa1)
                print(a2)
                print(msa2)
                return "Test Fail: Induced Pairwise Alignment not Correct"
        return "Test Success"

    def sp_score(self):
        """return sp score of the alignment"""
        ptr = self.alignments
        cost = 0
        while ptr.next:
            for i in range(self.len):
                for j in range(i + 1, self.len):
                    if ptr.next.values[i] == ptr.next.values[j]:
                        cost += 0
                    elif ptr.next.values[i] == '-' or ptr.next.values[j] == '-':
                        cost += gap
                    else:
                        cost += score[mapping[ptr.next.values[i]]][mapping[ptr.next.values[j]]]
            ptr = ptr.next
        return cost

    def length_of_alignment(self):
        """return the length of alignment"""
        count = 0
        ptr = self.alignments.next
        while ptr:
            count += 1
            ptr = ptr.next
        return count

    def average_sp_score(self):
        """return the sp score per length of the alignment"""
        count = self.length_of_alignment()
        score = self.sp_score()
        return score/count
