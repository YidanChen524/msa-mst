"""
Helper functions
"""

from config import gap, score, mapping


def pairwise_distance(s1, s2):
    """calculate the optimal sp_score between 2 sequences"""
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
    return T[m-1][n-1]


def sp_score(alignments):
    """calculate the sp score of alignments"""
    pass


def align_2_groups(group1, group2):
    """align 2 groups of sequences"""
    pass
