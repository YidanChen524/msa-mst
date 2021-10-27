import sys
from Bio import SeqIO
from config import gap, score, mapping


def sp_exact_3(seqs):
    # a helper function to calculate sp score for a single position
    def sp(i, j, k):
        s1 = s2 = s3 = None
        if i != '-' and j != '-':
            s1 = score[mapping[seqs[0][i - 1]]][mapping[seqs[1][j - 1]]]
        if i != '-' and k != '-':
            s2 = score[mapping[seqs[0][i - 1]]][mapping[seqs[2][k - 1]]]
        if j != '-' and k != '-':
            s3 = score[mapping[seqs[1][j - 1]]][mapping[seqs[2][k - 1]]]
        if i == '-' and j == '-':
            s1 = 0
        if i == '-' and k == '-':
            s2 = 0
        if j == '-' and k == '-':
            s3 = 0
        return sum([gap if s is None else s for s in [s1, s2, s3]])
    # define a dynamic table to store intermediate results
    l = len(seqs[0]) + 1
    m = len(seqs[1]) + 1
    n = len(seqs[2]) + 1
    T = [[[None for k in range(n)] for j in range(m)] for i in range(l)]
    # fill out the dynamic table T
    for i in range(l):
        for j in range(m):
            for k in range(n):
                v0 = v1 = v2 = v3 = v4 = v5 = v6 = v7 = None
                if i == 0 and j == 0 and k == 0:
                    v0 = 0
                if i > 0 and j > 0 and k > 0:
                    v1 = T[i - 1][j - 1][k - 1] + sp(i, j, k)
                if i > 0 and j > 0 and k >= 0:
                    v2 = T[i - 1][j - 1][k] + sp(i, j, '-')
                if i > 0 and j >= 0 and k > 0:
                    v3 = T[i - 1][j][k - 1] + sp(i, '-', k)
                if i >= 0 and j > 0 and k > 0:
                    v4 = T[i][j - 1][k - 1] + sp('-', j, k)
                if i > 0 and j >= 0 and k >= 0:
                    v5 = T[i - 1][j][k] + sp(i, '-', '-')
                if i >= 0 and j > 0 and k >= 0:
                    v6 = T[i][j - 1][k] + sp('-', j, '-')
                if i >= 0 and j >= 0 and k > 0:
                    v7 = T[i][j][k - 1] + sp('-', '-', k)
                T[i][j][k] = min([v for v in [v0, v1, v2, v3, v4, v5, v6, v7] if v is not None])
    # Find an optimal alignment
    a1 = a2 = a3 = ''
    i = l - 1
    j = m - 1
    k = n - 1
    while i != 0 or j != 0 or k != 0:
        v = T[i][j][k]
        if i > 0 and j > 0 and k > 0 and v == T[i - 1][j - 1][k - 1] + sp(i, j, k):
            a1 = seqs[0][i - 1] + a1
            a2 = seqs[1][j - 1] + a2
            a3 = seqs[2][k - 1] + a3
            i -= 1
            j -= 1
            k -= 1
        elif i > 0 and j > 0 and k >= 0 and v == T[i - 1][j - 1][k] + sp(i, j, '-'):
            a1 = seqs[0][i - 1] + a1
            a2 = seqs[1][j - 1] + a2
            a3 = '-' + a3
            i -= 1
            j -= 1
        elif i > 0 and j >= 0 and k > 0 and v == T[i - 1][j][k - 1] + sp(i, '-', k):
            a1 = seqs[0][i - 1] + a1
            a2 = '-' + a2
            a3 = seqs[2][k - 1] + a3
            i -= 1
            k -= 1
        elif i >= 0 and j > 0 and k > 0 and v == T[i][j - 1][k - 1] + sp('-', j, k):
            a1 = '-' + a1
            a2 = seqs[1][j - 1] + a2
            a3 = seqs[2][k - 1] + a3
            j -= 1
            k -= 1
        elif i > 0 and j >= 0 and k >= 0 and v == T[i - 1][j][k] + sp(i, '-', '-'):
            a1 = seqs[0][i - 1] + a1
            a2 = '-' + a2
            a3 = '-' + a3
            i -= 1
        elif i >= 0 and j > 0 and k >= 0 and v == T[i][j - 1][k] + sp('-', j, '-'):
            a1 = '-' + a1
            a2 = seqs[1][j - 1] + a2
            a3 = '-' + a3
            j -= 1
        elif i >= 0 and j >= 0 and k > 0 and v == T[i][j][k - 1] + sp('-', '-', k):
            a1 = '-' + a1
            a2 = '-' + a2
            a3 = seqs[2][k - 1] + a3
            k -= 1
    # return the optimal score and optimal alignment
    return T[l-1][m-1][n-1], [a1, a2, a3]


if __name__ == "__main__":
    # read in sequences
    names = []
    seqs = []
    for seq in SeqIO.parse(sys.argv[1], "fasta"):
        names.append(str(seq.id))
        seqs.append(str(seq.seq))
    # run exact algorithm
    optimal_score, optimal_alignment = sp_exact_3(seqs)
    # print results
    print("; optimal score: " + str(optimal_score) + "\n")
    print(f">{names[0]}\n" + optimal_alignment[0] + "\n")
    print(f">{names[1]}\n" + optimal_alignment[1] + "\n")
    print(f">{names[2]}\n" + optimal_alignment[2] + "\n")
    # write results to file
    with open("../archived/sp_exact_3_output.txt", "w") as f:
        f.write("; optimal score: " + str(optimal_score) + "\n\n")
        f.write(f">{names[0]}\n" + optimal_alignment[0] + "\n\n")
        f.write(f">{names[1]}\n" + optimal_alignment[1] + "\n\n")
        f.write(f">{names[2]}\n" + optimal_alignment[2] + "\n")
