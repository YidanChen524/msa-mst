import sys
from Bio import SeqIO
from config import gap, score, mapping


def sp_approx(seqs):
    # a helper function that calculating pairwise score and alignments
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
            return T[m-1][n-1], []
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
        return T[m-1][n-1], [a1, a2]
    # find the "center" string
    n_seqs = len(seqs)
    min_index = None
    min_sum = None
    for i in range(n_seqs):
        acc = 0
        for j in range(n_seqs):
            if j != i:
                acc += pairwise(seqs[i], seqs[j])[0]
        if min_sum is None or acc < min_sum:
            min_index = i
            min_sum = acc
    s = seqs[min_index]
    # Multiple Alignments
    for index in range(n_seqs):
        if index != min_index:
            m1 = seqs[min_index]
            [a1, a2] = pairwise(s, seqs[index], alignment=True)[1]
            i = j = k = 0
            while i < len(m1) and j < len(a1):
                if m1[i] == a1[j]:
                    i += 1
                    j += 1
                    k += 1
                elif m1[i] == '-':
                    a2 = a2[:k] + '-' + a2[k:]
                    i += 1
                    k += 1
                else:
                    for l in range(index):
                        seqs[l] = seqs[l][:k] + '-' + seqs[l][k:]
                    if min_index > index:
                        seqs[min_index] = seqs[min_index][:k] + '-' + seqs[min_index][k:]
                    j += 1
                    k += 1
            if i < len(m1):
                a2 = a2 + '-' * (len(m1) - i)
            if j < len(a1):
                for l in range(index):
                    seqs[l] = seqs[l] + '-' * (len(a1) - j)
                if min_index > index:
                    seqs[min_index] = seqs[min_index] + '-' * (len(a1) - j)
            seqs[index] = a2
    return min_index, seqs


if __name__ == "__main__":
    # read in sequences
    names = []
    seqs = []
    for seq in SeqIO.parse(sys.argv[1], "fasta"):
        names.append(str(seq.id))
        seqs.append(str(seq.seq))
    # run approx algorithm
    min_index, optimal_alignment = sp_approx(seqs)
    # print results
    print("; center string: " + names[min_index] + "\n")
    for i in range(len(seqs)):
        print(f">{names[i]}")
        print(optimal_alignment[i] + "\n")
    # write results to file
    with open("sp_approx_output.txt", "w") as f:
        f.write("; center string: " + names[min_index] + "\n\n")
        for i in range(len(seqs)):
            f.write(f">{names[i]}" + "\n")
            f.write(optimal_alignment[i] + "\n\n")
